use bstr::BString;
use std::collections::BTreeMap;
use std::error::Error;
use std::fs::File;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process;

use clap::Parser;
use hashbrown::{HashMap, HashSet};

use rayon::prelude::*;

use itertools::Itertools;
use num_iter::range_step_inclusive;

use noodles::bam;
use noodles::bgzf;
use noodles::sam;
// use noodles::bam::io::Writer as BamWriter;
// use noodles::sam::alignment::io::Write as SamWrite;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record_buf::data::field::Value as RecordBufDataValue;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::{map::program::tag as program_tag, map::Program};
use noodles::sam::header::Header as SamHeader;
use noodles::sam::header::{Builder as HeaderBuilder, ReferenceSequences};

// This lets us write `#[derive(Deserialize)]`.
use serde::Deserialize;

#[derive(Parser, Debug)]
#[command(author, version, long_about = None)]
#[command(name = "split_bams_per_cluster")]
#[command(about = "Split BAM files in per cluster BAM files based on cell barcodes.")]
#[command(
    about = "Split (multiple) BAM file(s) in per cluster BAM files based on a list of cell barcodes per cluster."
)]
struct Cli {
    #[arg(
        short = 's',
        long = "sample_bam",
        required = true,
        help = "Sample name to BAM filename mapping TSV file.",
        long_help = "Sample name to BAM filename mapping TSV file consisting of 2 columns:\n\
        \u{20} 1) sample: sample name (same name as used in cluster to cell barcode mapping TSV file)\n\
        \u{20} 2) bam_filename: BAM filename"
    )]
    sample_to_bam_tsv_path: PathBuf,
    #[arg(
        short = 'c',
        long = "cluster_cb_sample",
        required = true,
        help = "Cluster to original cell barcode, new cell barcode and sample name mapping TSV file.",
        long_help = "Cluster too original cell barcode, new cell barcode and sample name mapping TSV file consisting of 4 columns:\n\
        \u{20} 1) cluster: cluster\n\
        \u{20} 2) cell_barcode_input: input cell barcode (as written in input BAM files)\n\
        \u{20} 3) cell_barcode_output: output cell barcode (as to be written to output cluster BAM file)\n\
        \u{20} 4) sample: sample name (same name as used in sample to bam mapping TSV file)"
    )]
    cluster_to_cb_and_sample_tsv_path: PathBuf,
    #[arg(
        short = 'o',
        long = "output_prefix",
        required = true,
        help = "Output prefix.",
        long_help = "Output prefix used to create output cluster BAM files for each cluster."
    )]
    output_prefix: PathBuf,
    #[arg(
        short = 'f',
        long = "fragment_reads_only",
        required = false,
        help = "Only keep reads that will be used to create scATAC-seq fragments.",
        long_help = "Only keep reads that will be used to create scATAC-seq fragments:\n\
        \u{20} - read is properly paired.\n\
        \u{20} - read and its pair are located on the same chromosome.\n\
        \u{20} - read and its pair have a mapping quality of 30 or higher.\n\
        \u{20} - insert size is at least 10 in absolute value.\n\
        \u{20} - read is primary alignment."
    )]
    fragment_reads_only: bool,
}

// Sample name to BAM filename TSV record.
#[derive(Debug, Deserialize)]
struct SampleToBamRecord {
    sample: String,
    bam_filename: String,
}

// Cluster to orginal cell barcode, new cell barcode and sample name TSV file record.
#[derive(Debug, Deserialize)]
struct ClusterToCbSampleRecord {
    cluster: String,
    cell_barcode_input: String,
    cell_barcode_output: String,
    sample: String,
}

type BamToSampleMapping = HashMap<String, String>;
type BamToSampleBTreeMapping = BTreeMap<String, String>;

type SampleSet = HashSet<String>;
type ClusterToSamplesMapping = HashMap<String, SampleSet>;

#[derive(Debug)]
struct CellBarcodeOutputAndCluster {
    cell_barcode_output: String,
    cluster: String,
}

type CellBarcodeInputToCellBarcodeOutputAndClusterMapping =
    HashMap<String, CellBarcodeOutputAndCluster>;
type SampleToCellBarcodeInputToCellBarcodeOutputAndClusterMapping =
    HashMap<String, CellBarcodeInputToCellBarcodeOutputAndClusterMapping>;

type BamFileToBamReaderAndHeaderMapping =
    HashMap<String, (bam::io::Reader<bgzf::MultithreadedReader<File>>, SamHeader)>;

type ClusterToBamWriterAndHeaderMapping =
    HashMap<String, (bam::io::Writer<bgzf::MultithreadedWriter>, SamHeader)>;

fn read_sample_to_bam_tsv_file(
    sample_to_bam_tsv_path: &Path,
) -> Result<BamToSampleBTreeMapping, Box<dyn Error>> {
    let mut bam_to_sample_mapping: BamToSampleMapping = BamToSampleMapping::new();

    // Build a CSV reader for a plain TSV file.
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .escape(None)
        .double_quote(false)
        .quoting(false)
        .comment(Some(b'#'))
        .from_path(sample_to_bam_tsv_path)?;

    for result in rdr.deserialize() {
        let sample_to_bam_record: SampleToBamRecord = result?;

        // Add sample to BAM filename mapping to the hashmap.
        bam_to_sample_mapping
            .entry(sample_to_bam_record.bam_filename.clone())
            .or_insert(sample_to_bam_record.sample.clone());
    }

    // Create a BTreeMap from the BAM to sample HashMap to have a deterministic read
    // order in the output cluster BAM files for reads with the same start position.
    let bam_to_sample_mapping: BamToSampleBTreeMapping = bam_to_sample_mapping
        .iter()
        .sorted_by_key(|(k, _v)| *k)
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect::<_>();

    Ok(bam_to_sample_mapping)
}

fn read_cluster_to_cb_and_sample_tsv_file(
    cluster_to_cb_and_sample_tsv_path: &Path,
) -> Result<
    (
        ClusterToSamplesMapping,
        SampleToCellBarcodeInputToCellBarcodeOutputAndClusterMapping,
    ),
    Box<dyn Error>,
> {
    let mut cluster_to_samples_mapping: ClusterToSamplesMapping = ClusterToSamplesMapping::new();

    let mut sample_to_cb_input_to_cb_output_and_cluster_mapping: SampleToCellBarcodeInputToCellBarcodeOutputAndClusterMapping = SampleToCellBarcodeInputToCellBarcodeOutputAndClusterMapping::new();

    // Build a CSV reader for a plain TSV file.
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .escape(None)
        .double_quote(false)
        .quoting(false)
        .comment(Some(b'#'))
        .from_path(cluster_to_cb_and_sample_tsv_path)?;

    for result in rdr.deserialize() {
        let cluster_to_cb_sample_record: ClusterToCbSampleRecord = result?;

        // Create cluster to samples mapping.
        cluster_to_samples_mapping
            .entry(cluster_to_cb_sample_record.cluster.clone())
            .or_default()
            .insert(cluster_to_cb_sample_record.sample.clone());

        // Create nested hashmap with multiple levels:
        //   - level 1: Sample name to cell barcode input mapping
        //   - level 2: Cell barcode input to cell barcode output mapping
        //   - level 3: Cell barcode output to cluster mapping
        sample_to_cb_input_to_cb_output_and_cluster_mapping
            .entry(cluster_to_cb_sample_record.sample.clone())
            .or_default()
            .entry(cluster_to_cb_sample_record.cell_barcode_input.clone())
            .or_insert(CellBarcodeOutputAndCluster {
                cell_barcode_output: cluster_to_cb_sample_record.cell_barcode_output.clone(),
                cluster: cluster_to_cb_sample_record.cluster.clone(),
            });
    }

    Ok((
        cluster_to_samples_mapping,
        sample_to_cb_input_to_cb_output_and_cluster_mapping,
    ))
}

fn fix_pg_bam_header_lines<'a>(
    mut merged_header_builder: HeaderBuilder,
    header: &'a SamHeader,
    sample_id: &'a str,
) -> HeaderBuilder {
    // Add sample name to "ID" and "PP" fields to make them unique over all sample BAM files.
    for (program_id, mut program_map) in header.programs().clone() {
        // Add sample name to "ID" field in "@PG" line.
        let mut program_id = program_id.clone();

        program_id.extend(format!(".{}", sample_id).into_bytes());

        // Add sample name to "PP" field in "@PG" line.
        if let Some(pp_value) = program_map.other_fields_mut().get_mut(b"PP") {
            let mut pp_field = pp_value.clone();
            pp_field.extend(format!(".{}", sample_id).into_bytes());
            *pp_value = pp_field;
        }

        merged_header_builder = merged_header_builder.add_program(program_id, program_map.clone());
    }

    merged_header_builder
}

fn split_bams_per_cluster(
    bam_to_sample_mapping: &BamToSampleBTreeMapping,
    cluster_to_samples_mapping: &ClusterToSamplesMapping,
    sample_to_cb_input_to_cb_output_and_cluster_mapping: &SampleToCellBarcodeInputToCellBarcodeOutputAndClusterMapping,
    output_prefix: &Path,
    fragment_reads_only: bool,
    cmd_line_str: &str,
) -> Result<(), Box<dyn Error>> {
    let mut bam_file_to_bam_reader_and_header_mapping = BamFileToBamReaderAndHeaderMapping::new();
    let mut cluster_to_bam_writer_and_header_mapping = ClusterToBamWriterAndHeaderMapping::new();

    // Get all cluster names and sort them.
    let mut clusters: Vec<&String> = cluster_to_samples_mapping.keys().collect();
    clusters.sort_unstable();

    let mut cluster_to_bam_records: HashMap<String, Vec<RecordBuf>> = HashMap::new();
    let mut cluster_to_out_of_range_bam_records: HashMap<String, Vec<RecordBuf>> = HashMap::new();

    let mut reference_sequences: Option<ReferenceSequences> = None;

    // Construct BAM header for each per cluster BAM file by combining headers
    // of each per sample BAM file.
    for cluster in clusters.iter() {
        // Get all bam filenames per cluster and sort them.
        let mut bam_filenames = cluster_to_samples_mapping
            .get(cluster.as_str())
            .unwrap()
            .iter()
            .flat_map(|sample| {
                bam_to_sample_mapping
                    .iter()
                    .filter(|(_bam_filename, s)| *s == sample)
                    .map(|(bam_filename, _sample)| bam_filename)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        bam_filenames.sort_unstable();

        // Create merged BAM header per cluster BAM file.
        let mut cluster_bam_header_builder = sam::Header::builder();

        for (i, bam_filename) in bam_filenames.iter().enumerate() {
            let sample_id = bam_to_sample_mapping.get(*bam_filename).unwrap();
            println!(
                "cluster: {}, bam_filename: {}, sample_id: {}",
                cluster, bam_filename, sample_id
            );

            let (ref mut bam_reader, header) = bam_file_to_bam_reader_and_header_mapping
                .entry(bam_filename.to_string())
                .or_insert({
                    // Read input BAM file with 4 decompressing worker threads.
                    let bam_decompressing_worker_threads = NonZeroUsize::new(4).unwrap();
                    let bam_file = File::open(bam_filename)?;

                    let bam_bgzf_reader = bgzf::MultithreadedReader::with_worker_count(
                        bam_decompressing_worker_threads,
                        bam_file,
                    );

                    // After decompressing BGZF blocks of BAM file in parallel,
                    // decode the actual BAM records.
                    let mut bam_reader = bam::io::Reader::from(bam_bgzf_reader);

                    // Read header.
                    let header = bam_reader.read_header()?;

                    (bam_reader, header)
                });

            header
                .reference_sequences()
                .iter()
                .for_each(|reference_sequence| {
                    let reference_sequence_name =
                        String::from_utf8(reference_sequence.0.to_vec()).unwrap();
                    let reference_sequence_length = reference_sequence.1.length();
                    // println!("{}: {}", reference_sequence_name, reference_sequence_length);
                });
            // println!("\nheader\n{:?}\n", header.header());
            // println!("\nread_groups:\n{:?}\n", header.read_groups());
            // println!("\nprograms:\n{:?}\n", header.programs());
            // println!("\ncomments:\n{:?}\n", header.comments());

            // Add "@HD" and "@SQ" lines from first BAM file to merged header.
            if i == 0 {
                cluster_bam_header_builder = cluster_bam_header_builder
                    .set_header(header.header().unwrap().clone())
                    .set_reference_sequences(header.reference_sequences().clone());

                // Store reference sequences for later use to get chromosome names and lengths.
                reference_sequences = Some(header.reference_sequences().clone());
            }

            // Add "@RG" lines from each BAM file to merged header.
            if !header.read_groups().is_empty() {
                for (read_group_id, read_group_map) in header.read_groups() {
                    cluster_bam_header_builder = cluster_bam_header_builder
                        .add_read_group(read_group_id.clone(), read_group_map.clone());
                }
            }

            // Add "@PG"" lines from each BAM file to merged header.
            if !header.programs().is_empty() {
                cluster_bam_header_builder =
                    fix_pg_bam_header_lines(cluster_bam_header_builder, &header, sample_id);
            }

            // Add "@CO" lines from each BAM file to merged header.
            if !header.comments().is_empty() {
                for comment in header.comments() {
                    cluster_bam_header_builder =
                        cluster_bam_header_builder.add_comment(comment.clone());
                }
            }
        }

        // Add "@PG" header line for "split_bams_per_cluster".
        cluster_bam_header_builder =
            cluster_bam_header_builder.add_program("split_bams_per_cluster", {
                let mut program_map = Map::<Program>::default();
                program_map
                    .other_fields_mut()
                    .insert(program_tag::NAME, BString::from("split_bams_per_cluster"));
                program_map.other_fields_mut().insert(
                    program_tag::VERSION,
                    BString::from(env!("CARGO_PKG_VERSION")),
                );
                program_map
                    .other_fields_mut()
                    .insert(program_tag::COMMAND_LINE, BString::from(cmd_line_str));
                program_map
            });

        let cluster_bam_header = cluster_bam_header_builder.build();

        // Create per cluster BAM file writer.

        let mut cluster_bam_path = PathBuf::from(output_prefix);
        cluster_bam_path
            .as_mut_os_string()
            .push(format!("{}.bam", cluster));

        // Write records to BGZF-compressed BAM file file with 3 worker threads.
        // Set the compression level to 7, to match the default compression level (6)
        // of bgzip of HTSlib (which has libdeflate level 7 mapped to level 6:
        // https://github.com/samtools/htslib/pull/1488).
        let cluster_bam_file = File::create(cluster_bam_path)?;
        let compression_level = bgzf::writer::CompressionLevel::try_from(7)?;
        let cluster_bam_compressing_worker_threads = NonZeroUsize::new(3).unwrap();

        let mut cluster_bam_bgzf_writer = bgzf::multithreaded_writer::Builder::default()
            .set_compression_level(compression_level)
            .set_worker_count(cluster_bam_compressing_worker_threads)
            .build_from_writer(cluster_bam_file);

        // After decompressing BGZF blocks of BAM file in parallel,
        // decode the actual BAM records.
        let mut cluster_bam_writer = bam::io::Writer::from(cluster_bam_bgzf_writer);

        cluster_bam_writer.write_header(&cluster_bam_header)?;

        cluster_to_bam_writer_and_header_mapping
            .entry(cluster.to_string())
            .or_insert((cluster_bam_writer, cluster_bam_header));

        // Create empty vector to store BAM records for current cluster.
        cluster_to_bam_records
            .entry(cluster.to_string())
            .or_default();
        // Create empty vector to store BAM records for current cluster that do not fit in the current region.
        cluster_to_out_of_range_bam_records
            .entry(cluster.to_string())
            .or_default();
    }

    // println!("refence_sequences: {:?}", reference_sequences);
    let reference_sequences = reference_sequences.unwrap();

    // Filter out BAM files that are not in any requested cluster.
    let bam_to_sample_mapping: BamToSampleBTreeMapping = bam_to_sample_mapping
        .iter()
        .filter(|(k, _v)| bam_file_to_bam_reader_and_header_mapping.get(*k).is_some())
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();

    let chunk_size: i64 = 10_000_000;

    // Loop over each chromosome and fetch reads in chunks from each BAM file and sort
    // them by position before writing them to the per cluster BAM file.
    for (reference_sequence_id, reference_sequence) in reference_sequences.iter().enumerate() {
        let chrom_end = reference_sequence.1.length().get() as i64;

        // Fetch reads from each BAM file for each cluster in chunks of 10_000_000 bp and sort them by position.
        for (start, end) in
            range_step_inclusive(1, chrom_end + chunk_size - 1, chunk_size).tuple_windows()
        {
            // Make sure that the end of the chunk is not larger than the end of the chromosome.
            let end = if end > chrom_end { chrom_end } else { end };

            println!(
                "reference_sequence_id: {}, start: {}, end: {}, chrom_end: {}",
                reference_sequence_id, start, end, chrom_end
            );

            println!(
                "cluster_to_out_of_range_bam_records = {:?}",
                cluster_to_out_of_range_bam_records
            );

            // Check which out of range records from previous iteration fit in the
            // current interval and add them to the correct cluster BAM record vector.
            cluster_to_out_of_range_bam_records.iter_mut().for_each(
                |(cluster, cluster_out_of_region_bam_records)| {
                    cluster_out_of_region_bam_records.retain(|record| {
                        println!("recovered record: {:?}", record);
                        if let Some(cluster_bam_records) =
                            cluster_to_bam_records.get_mut(cluster.as_str())
                        {
                            if record.reference_sequence_id().unwrap() == reference_sequence_id
                                && record.alignment_start().unwrap().get() <= end as usize
                            {
                                println!("placed.");
                                // Add BAM record to correct per cluster BAM record vector.
                                cluster_bam_records.push(record.clone());
                                false
                            } else {
                                println!("not placed yet.");
                                // Keep BAM record in out of range cluster BAM record vector for next iteration.
                                true
                            }
                        } else {
                            // This should never happen as both cluster_to_bam_records and cluster_to_out_of_range_bam_records
                            // contain the same cluster names as keys.
                            true
                        }
                    })
                },
            );

            bam_to_sample_mapping
                .iter()
                .try_for_each(|(bam_filename, sample)| {
                    println!("bam_filename: {}", bam_filename);
                    let (bam_reader, header) = bam_file_to_bam_reader_and_header_mapping
                        .get_mut(bam_filename)
                        .unwrap();

                    let cb_input_to_cb_output_and_cluster_mapping =
                        sample_to_cb_input_to_cb_output_and_cluster_mapping
                            .get(sample)
                            .unwrap();

                    // Filter reads of current chunk and write them to a per cluster vector.
                    for (i, result) in bam_reader.records().enumerate() {
                        let record = result?;

                        // Keep only read that will be used to create scATAC-seq fragment or all
                        // reads, depending on the settings.
                        let keep_read = match fragment_reads_only {
                            true => {
                                let mut keep_read = false;

                                // Only keep reads that will be used to create scATAC-seq fragments.
                                //   - read is properly paired.
                                //   - read and its pair are located on the same chromosome.
                                //   - read and its pair have a mapping quality of 30 or higher.
                                //   - insert size is at least 10 in absolute value.
                                //   - read is primary alignment.
                                if record.flags().is_properly_segmented()
                                    && record.reference_sequence_id().unwrap().unwrap()
                                        == record.mate_reference_sequence_id().unwrap().unwrap()
                                    && record.mapping_quality() >= MappingQuality::new(30)
                                    && record.template_length().abs() >= 10
                                    && !record.flags().is_secondary()
                                    && !record.flags().is_supplementary()
                                {
                                    if let Some(mate_mapq) =
                                        record.data().get(&Tag::MATE_MAPPING_QUALITY).transpose()?
                                    {
                                        if mate_mapq.as_int() >= Some(30) {
                                            // Keep current BAM record.
                                            keep_read = true;
                                        }
                                    }
                                }

                                keep_read
                            }
                            false => true,
                        };

                        if !keep_read {
                            // Skip unwanted reads.
                            continue;
                        }

                        if record.reference_sequence_id().is_none() {
                            continue;
                        }
                        if record.alignment_start().is_none() {
                            continue;
                        }

                        _ = if let Some(Value::String(cb)) =
                            record.data().get(&Tag::CELL_BARCODE_ID).transpose()?
                        {
                            // Add BAM record with updated full barcode name to per
                            // cluster vector, if the barcode was in the list of
                            // filtered barcodes.
                            if let Some(cb_output_and_cluster) =
                                cb_input_to_cb_output_and_cluster_mapping.get(&cb.to_string())
                            {
                                let CellBarcodeOutputAndCluster {
                                    cell_barcode_output: cb_output,
                                    cluster,
                                } = cb_output_and_cluster;

                                let mut record_buf = RecordBuf::try_from_alignment_record(
                                    &header, &record,
                                )?;
                                let data = record_buf.data_mut();

                                data.insert(
                                    Tag::CELL_BARCODE_ID,
                                    RecordBufDataValue::String(BString::from(
                                        cb_output.to_string(),
                                    )),
                                );

                                // Add current BAM record to correct per cluster vector.
                                if record.reference_sequence_id().unwrap()?
                                    != reference_sequence_id
                                    || record.alignment_start().unwrap()?.get() > end as usize
                                {
                                    println!("check: record {}: {:?}, reference_sequence_id: {} != {}, end: {}", i, record, record.reference_sequence_id().unwrap()?, reference_sequence_id, record.alignment_start().unwrap()?.get());
                                    if let Some(cluster_out_of_region_bam_records) =
                                        cluster_to_out_of_range_bam_records
                                            .get_mut(cluster.as_str())
                                    {
                                        cluster_out_of_region_bam_records.push(record_buf);
                                        // Skip reads that are not in the current region, which might
                                        // be pulled in when using fetch as this might result in
                                        // duplicated reads in the output.
                                        break;
                                    }
                                } else {
                                    // Add current BAM record to correct per cluster vector.
                                    if let Some(cluster_bam_records) =
                                        cluster_to_bam_records.get_mut(cluster)
                                    {
                                        cluster_bam_records.push(record_buf);
                                    }
                                }
                            }
                        }
                    }

                    Ok::<_, Box<dyn std::error::Error>>(())
                })?;

            // Sort reads by position (for the current chunk) for each cluster.
            cluster_to_bam_records
                .par_iter_mut()
                .for_each(|(_cluster, cluster_bam_records)| {
                    cluster_bam_records.sort_by(|a, b| {
                        a.reference_sequence_id()
                            .unwrap()
                            .cmp(&b.reference_sequence_id().unwrap())
                            .then(
                                a.alignment_start()
                                    .unwrap()
                                    .get()
                                    .cmp(&b.alignment_start().unwrap().get()),
                            )
                    });
                });

            // Write sorted reads to per cluster BAM files.
            // cluster_to_bam_records
            //     .iter_mut()
            //     .try_for_each(|(cluster, cluster_bam_records)| {
            //         println!(
            //             "Writing chunk {}:{}-{} for cluster {}",
            //             reference_sequence_id, start, end, cluster
            //         );
            //         let (cluster_bam_writer, cluster_header) =
            //             cluster_to_bam_writer_and_header_mapping
            //                 .get_mut(cluster)
            //                 .unwrap();

            //         for record in cluster_bam_records.iter() {
            //             cluster_bam_writer.write_alignment_record(&cluster_header, record)?;
            //         }

            //         cluster_bam_records.clear();

            //         Ok::<_, Box<dyn std::error::Error>>(())
            //     })?;
            cluster_to_bam_records
                .par_iter_mut()
                .try_for_each(|(cluster, cluster_bam_records)| {
                    println!(
                        "Writing chunk {}:{}-{} for cluster {}",
                        reference_sequence_id, start, end, cluster
                    );
                    let (cluster_bam_writer, cluster_header) =
                        cluster_to_bam_writer_and_header_mapping
                            .get_mut(cluster)
                            .unwrap();

                    for record in cluster_bam_records.iter() {
                        cluster_bam_writer
                            .write_alignment_record(&cluster_header, record)
                            .unwrap();
                    }

                    cluster_bam_records.clear();

                    Ok::<_, Box<dyn std::error::Error + Send>>(())
                })
                .unwrap();
        }
    }

    cluster_to_bam_records
        .iter_mut()
        .try_for_each(|(cluster, cluster_bam_records)| {
            let (cluster_bam_writer, cluster_header) = cluster_to_bam_writer_and_header_mapping
                .get_mut(cluster)
                .unwrap();
            cluster_bam_writer.finish(&cluster_header)?;

            Ok::<_, Box<dyn std::error::Error>>(())
        })?;

    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let cmd_line_str = format!(
        "split_bams_per_cluster -s {} -c {} -o {}",
        &cli.sample_to_bam_tsv_path.to_string_lossy(),
        &cli.cluster_to_cb_and_sample_tsv_path.to_string_lossy(),
        &cli.output_prefix.to_string_lossy()
    );

    let bam_to_sample_mapping = match read_sample_to_bam_tsv_file(&cli.sample_to_bam_tsv_path) {
        Ok(sample_to_bam_mapping) => sample_to_bam_mapping,
        Err(e) => {
            println!("Error: {}", e);
            process::exit(1);
        }
    };

    let (cluster_to_samples_mapping, sample_to_cb_input_to_cb_output_and_cluster_mapping) =
        match read_cluster_to_cb_and_sample_tsv_file(&cli.cluster_to_cb_and_sample_tsv_path) {
            Ok(cluster_to_cb_sample_mapping) => cluster_to_cb_sample_mapping,
            Err(e) => {
                println!("Error: {}", e);
                process::exit(1);
            }
        };

    match split_bams_per_cluster(
        &bam_to_sample_mapping,
        &cluster_to_samples_mapping,
        &sample_to_cb_input_to_cb_output_and_cluster_mapping,
        &cli.output_prefix,
        cli.fragment_reads_only,
        &cmd_line_str,
    ) {
        Ok(()) => (),
        Err(e) => {
            println!("Error: {}", e);
            process::exit(1);
        }
    };
}
