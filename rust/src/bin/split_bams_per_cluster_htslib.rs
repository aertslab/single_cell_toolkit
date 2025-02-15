use std::collections::BTreeMap;
use std::error::Error;
use std::path::{Path, PathBuf};
use std::process;

use clap::Parser;
use hashbrown::{HashMap, HashSet};

use rayon::prelude::*;

use itertools::Itertools;
use num_iter::range_step_inclusive;

use rust_htslib::bam::{
    record::Aux, CompressionLevel, Format, Header, HeaderView, IndexedReader, Read, Record, Writer,
};
use rust_htslib::tpool::ThreadPool;

// This lets us write `#[derive(Deserialize)]`.
use serde::Deserialize;

#[derive(Parser, Debug)]
#[command(author, version, long_about = None)]
#[command(name = "split_bams_per_cluster_htslib")]
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
        long = "chroms",
        num_args(0..),
        required = false,
        help = "List of chromosome names to keep reads for in the output BAM files.",
        long_help = "List of chromosome names to keep reads for in the output BAM files.\n\
        If not specified, keep reads for all chromosomes in the output BAM files."
    )]
    chromosomes: Option<Vec<String>>,
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
    #[arg(
        long = "ignore_mate_mapping_quality",
        required = false,
        help = "Ignore mate mapping quality when filtering reads with `--fragment_reads_only`.",
        long_help = "Ignore mate mapping quality when filtering reads with `--fragment_reads_only`.\n\
        Use this option if reads do not contain a `MQ` tag.\n\
        `MQ` tags can be added to reads with `samtools fixmate -m`."
    )]
    ignore_mate_mapping_quality: bool,
    #[arg(
        short = 'C',
        long = "chunk_size",
        required = false,
        default_value_t = 1_000_000,
        help = "Fetch reads from each BAM file in chunks of X bp. Default: 1_000_000.",
        long_help = "Fetch reads from each BAM file for each cluster in chunks of X bp\n\
        (default: 1_000_000) and sort them by position.\n\
        Reduce this value if split_bams_per_cluster_htslib uses too much memory."
    )]
    chunk_size: u64,
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

struct CellBarcodeOutputAndCluster {
    cell_barcode_output: String,
    cluster: String,
}

type CellBarcodeInputToCellBarcodeOutputAndClusterMapping =
    HashMap<String, CellBarcodeOutputAndCluster>;
type SampleToCellBarcodeInputToCellBarcodeOutputAndClusterMapping =
    HashMap<String, CellBarcodeInputToCellBarcodeOutputAndClusterMapping>;

type BamFileToBamIndexedReaderMapping = HashMap<String, IndexedReader>;
type ClusterToBamWriterMapping = HashMap<String, Writer>;

trait ContainsSlice<T> : PartialEq<[T]> {
    fn contains_slice (self: &'_ Self, slice: &'_ [T]) -> bool;
}

impl<T, Item : PartialEq<T>> ContainsSlice<T> for [Item] {
    fn contains_slice (self: &'_ [Item], slice: &'_ [T]) -> bool
    {
        let len = slice.len();
        if len == 0 {
            return true;
        }
        self.windows(len)
            .any(move |sub_slice| sub_slice == slice)
    }
}

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

fn get_hd_and_sq_bam_header_lines(header: &Header) -> Vec<u8> {
    // Only keep "@HD" and "@SQ" lines from BAM header.
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && (x.starts_with(b"@HD\t") || x.starts_with(b"@SQ\t")))
        .collect::<Vec<_>>()
        .join(&b"\n"[..])
}

fn get_hd_coordinate_sorted_bam_header_lines(header: &Header) -> bool {
    // Only keep "@HD" and "@SQ" lines from BAM header.
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && (x.starts_with(b"@HD\t") && x.contains_slice(b"\tSO:coordinate")))
        .collect::<Vec<_>>()
        .is_empty()
}

fn get_sq_bam_header_lines(header: &Header) -> Vec<u8> {
    // Only keep "@HD" and "@SQ" lines from BAM header.
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && (x.starts_with(b"@SQ\t")))
        .collect::<Vec<_>>()
        .join(&b"\n"[..])
}

fn get_non_hd_sq_and_fix_pg_bam_header_lines(header: &Header, sample: &str) -> Vec<u8> {
    // Only keep non-"@HD" and non-"@SQ" lines (@PG and @CO lines) from
    // BAM header and make "ID" and "PP" fields in "@PG" lines unique,
    // by adding sample name.
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && !x.starts_with(b"@HD\t") && !x.starts_with(b"@SQ\t"))
        .map(Vec::<u8>::from)
        .map(|x| {
            match x.starts_with(b"@PG\t") {
                true => {
                    // Split "@PG" line by tab and add sample name to "ID" and "PP" fields
                    // to make them unique over all sample BAM files.
                    x.split(|x| x == &b'\t')
                        .map(|x| {
                            match x.starts_with(b"ID:") {
                                true => {
                                    // Add sample name to "ID" field in "@PG" line.
                                    let mut id_field = Vec::<u8>::from(x);
                                    id_field.extend(format!(".{}", &sample).into_bytes());
                                    id_field
                                }
                                false => x.to_vec(),
                            }
                        })
                        .map(|x| {
                            match x.starts_with(b"PP:") {
                                true => {
                                    // Add sample name to "PP" field in "@PG" line.
                                    let mut pp_field = x;
                                    pp_field.extend(format!(".{}", &sample).into_bytes());
                                    pp_field
                                }
                                false => x.to_vec(),
                            }
                        })
                        .collect::<Vec<Vec<u8>>>()
                        .join(&b'\t')
                }
                false => x,
            }
        })
        .collect::<Vec<Vec<u8>>>()
        .join(&b"\n"[..])
}

fn split_bams_per_cluster(
    bam_to_sample_mapping: &BamToSampleBTreeMapping,
    cluster_to_samples_mapping: &ClusterToSamplesMapping,
    sample_to_cb_input_to_cb_output_and_cluster_mapping: &SampleToCellBarcodeInputToCellBarcodeOutputAndClusterMapping,
    output_prefix: &Path,
    chromosomes: &Option<Vec<String>>,
    fragment_reads_only: bool,
    ignore_mate_mapping_quality: bool,
    chunk_size: u64,
    cmd_line_str: &str,
) -> Result<(), Box<dyn Error>> {
    let bam_thread_pool = ThreadPool::new(16)?;
    let mut bam_file_to_bam_indexed_reader_mapping = BamFileToBamIndexedReaderMapping::new();
    let mut cluster_to_bam_writer_mapping = ClusterToBamWriterMapping::new();

    // Get all cluster names and sort them.
    let mut clusters: Vec<&String> = cluster_to_samples_mapping.keys().collect();
    clusters.sort_unstable();

    let mut cluster_to_bam_records: HashMap<String, Vec<Record>> = HashMap::new();

    let mut merged_header_view: Option<HeaderView> = None;

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
        let mut merged_header = Vec::new();
        let mut hd_and_sq_bam_header_lines: Option<Vec<u8>> = None;
        let mut sq_bam_header_lines: Option<Vec<u8>> = None;

        bam_filenames.iter().enumerate().try_for_each(|(i, bam_filename)| {
            let bam = bam_file_to_bam_indexed_reader_mapping
                .entry(bam_filename.to_string())
                .or_insert({
                    let mut indexed_input_bam = IndexedReader::from_path(Path::new(bam_filename))?;
                    indexed_input_bam.set_thread_pool(&bam_thread_pool)?;
                    indexed_input_bam
                });

            // Read BAM header from current sample BAM file.
            let original_header = Header::from_template(bam.header());

            // Get "@HD" and "@SQ" lines and check if they match exactly with those
            // lines in the first sample BAM file.
            match i {
                0 => {
                    if get_hd_coordinate_sorted_bam_header_lines(&original_header) {
                        Err(format!("BAM file \"{}\" for cluster \"{}\" is not coordinate sorted.", bam_filename, cluster))?
                    }

                    hd_and_sq_bam_header_lines =
                        Some(get_hd_and_sq_bam_header_lines(&original_header));
                    merged_header.extend(&hd_and_sq_bam_header_lines.clone().unwrap());

                    sq_bam_header_lines =
                        Some(get_sq_bam_header_lines(&original_header));
                }
                _ => {
                    if get_hd_coordinate_sorted_bam_header_lines(&original_header) {
                        Err(format!("BAM file \"{}\" for cluster \"{}\" is not coordinate sorted.", bam_filename, cluster))?
                    }

                    if get_sq_bam_header_lines(&original_header)
                        != sq_bam_header_lines.clone().unwrap()
                    {
                        Err(format!("BAM file \"{}\" for cluster \"{}\" has different chromosome order.", bam_filename, cluster))?
                    }
                }
            }

            // Get all "@PG", "@CO" and "@RG" lines.
            let non_hd_sq_and_fix_pg_bam_header_lines = get_non_hd_sq_and_fix_pg_bam_header_lines(
                &original_header,
                bam_to_sample_mapping.get(bam_filename.as_str()).unwrap(),
            );

            merged_header.extend(&b"\n"[..]);
            merged_header.extend(&non_hd_sq_and_fix_pg_bam_header_lines);

            Ok::<(), Box<dyn Error>>(())
        })?;

        // Add "@PG" header line for "split_bams_per_cluster_htslib".
        let pg_header_line = format!(
            "\n@PG\tID:split_bams_per_cluster_htslib\tPN:split_bams_per_cluster_htslib VN:{} CL:{}",
            env!("CARGO_PKG_VERSION"),
            &cmd_line_str
        )
        .into_bytes();
        merged_header.extend(&pg_header_line);

        let merged_header = Header::from_template(&HeaderView::from_bytes(&merged_header));

        // Create per cluster BAM file writer.
        let mut cluster_bam_path = PathBuf::from(output_prefix);
        cluster_bam_path
            .as_mut_os_string()
            .push(format!("{}.bam", cluster));

        let mut cluster_bam_writer =
            Writer::from_path(&cluster_bam_path, &merged_header, Format::Bam)?;

        cluster_bam_writer.set_thread_pool(&bam_thread_pool)?;
        cluster_bam_writer.set_compression_level(CompressionLevel::Fastest)?;
        cluster_to_bam_writer_mapping
            .entry(cluster.to_string())
            .or_insert(cluster_bam_writer);

        // Create empty vector to store BAM records for current cluster.
        cluster_to_bam_records
            .entry(cluster.to_string())
            .or_default();

        // Store merged header view for first cluster for later use to get chromosome names and lengths.
        if merged_header_view.is_none() {
            merged_header_view = Some(HeaderView::from_header(&merged_header));
        }
    }

    let merged_header_view = merged_header_view.unwrap();

    // Filter out BAM files that are not in any requested cluster.
    let bam_to_sample_mapping: BamToSampleBTreeMapping = bam_to_sample_mapping
        .iter()
        .filter(|(k, _v)| bam_file_to_bam_indexed_reader_mapping.get(*k).is_some())
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();

    let cb_tag = b"CB";

    // Loop over each chromosome and fetch reads in chunks from each BAM file and sort
    // them by position before writing them to the per cluster BAM file.
    for tid in 0..merged_header_view.target_count() {
        let chrom_name = std::str::from_utf8(merged_header_view.tid2name(tid)).expect("Chromosome name is not valid UTF-8.").to_string();
        if chromosomes.is_some() {
            if ! chromosomes.as_ref().unwrap().contains(&chrom_name) {
                continue;
            }
        }

        let chrom_end = merged_header_view.target_len(tid).unwrap();

        // Fetch reads from each BAM file for each cluster in chunks of 10_000_000 bp and sort them by position.
        for (start, end) in
            range_step_inclusive(0, chrom_end + chunk_size - 1, chunk_size).tuple_windows()
        {
            // Make sure that the end of the chunk is not larger than the end of the chromosome.
            let end = if end > chrom_end { chrom_end } else { end };

            let start = start as i64;
            let end = end as i64;

            bam_to_sample_mapping
                .iter()
                .try_for_each(|(bam_filename, sample)| {
                    let indexed_input_bam = bam_file_to_bam_indexed_reader_mapping
                        .get_mut(bam_filename)
                        .unwrap();
                    // Fetch chunk from current BAM file (coordinates are 0-based, and end is exclusive).
                    indexed_input_bam.fetch((tid, start, end))?;

                    let sample_to_cb_input_to_cb_output_and_cluster_mapping =
                        sample_to_cb_input_to_cb_output_and_cluster_mapping
                            .get(sample)
                            .unwrap();

                    // Filter reads of current chunk and write them to a per cluster vector.
                    for r in indexed_input_bam.records() {
                        let mut record = r?;

                        if record.pos() < start || record.pos() > end || record.tid() != tid as i32
                        {
                            // Skip reads that are not in the current region, which might
                            // be pulled in when using fetch as this might result in
                            // duplicated reads in the output.
                            continue;
                        }

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
                                if record.is_proper_pair()
                                    && record.tid() == record.mtid()
                                    && record.mapq() >= 30
                                    && record.insert_size().abs() >= 10
                                    && !record.is_secondary()
                                    && !record.is_supplementary()
                                {
                                    if !ignore_mate_mapping_quality {
                                        if let Ok(mate_mapq_aux) = record.aux(b"MQ") {
                                            // So far MQ tags in BAM files have been of I8 or U8 type.
                                            let mate_mapq = match mate_mapq_aux {
                                                Aux::I8(mate_mapq) => mate_mapq as i32,
                                                Aux::I16(mate_mapq) => mate_mapq as i32,
                                                Aux::I32(mate_mapq) => mate_mapq,
                                                Aux::U8(mate_mapq) => mate_mapq as i32,
                                                Aux::U16(mate_mapq) => mate_mapq as i32,
                                                Aux::U32(mate_mapq) => mate_mapq as i32,
                                                _ => -1, // bail!("Value for MQ tag is not an integer."),
                                            };

                                            if mate_mapq >= 30 {
                                                // Keep current BAM record.
                                                keep_read = true;
                                            }
                                        }
                                    } else {
                                        // Keep current BAM record regardless of mate mapping quality.
                                        keep_read = true;
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

                        if let Ok(Aux::String(cb)) = record.aux(cb_tag) {
                            // Add BAM record with updated full barcode name to per
                            // cluster vector, if the barcode was in the list of
                            // filtered barcodes.
                            if let Some(cb_output_and_cluster) =
                                sample_to_cb_input_to_cb_output_and_cluster_mapping.get(cb)
                            {
                                let CellBarcodeOutputAndCluster {
                                    cell_barcode_output: cb_output,
                                    cluster,
                                } = cb_output_and_cluster;

                                // Update CB tag value with full barcode name.
                                record.remove_aux(cb_tag)?;
                                record.push_aux(cb_tag, Aux::String(cb_output))?;

                                // Add current BAM record to correct per cluster vector.
                                if let Some(cluster_bam_records) =
                                    cluster_to_bam_records.get_mut(cluster.as_str())
                                {
                                    cluster_bam_records.push(record);
                                }
                            }
                        }
                    }

                    Ok::<(), rust_htslib::errors::Error>(())
                })?;

            // Sort reads by position (for the current chunk) for each cluster.
            cluster_to_bam_records
                .par_iter_mut()
                .for_each(|(_cluster, cluster_bam_records)| {
                    cluster_bam_records
                        .sort_by(|a, b| a.tid().cmp(&b.tid()).then(a.pos().cmp(&b.pos())));
                });

            // Write sorted reads to per cluster BAM files.
            cluster_to_bam_records
                .iter_mut()
                .try_for_each(|(cluster, cluster_bam_records)| {
                    println!(
                        "Writing chunk {}:{}-{} for cluster {}",
                        chrom_name, start, end, cluster
                    );
                    let cluster_bam_writer =
                        cluster_to_bam_writer_mapping.get_mut(cluster).unwrap();

                    for record in cluster_bam_records.iter() {
                        cluster_bam_writer.write(record)?;
                    }

                    cluster_bam_records.clear();

                    Ok::<(), rust_htslib::errors::Error>(())
                })?;
        }
    }

    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let cmd_line_str = format!(
        "split_bams_per_cluster_htslib -s {} -c {} -o {}",
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
        &cli.chromosomes,
        cli.fragment_reads_only,
        cli.ignore_mate_mapping_quality,
        cli.chunk_size,
        &cmd_line_str,
    ) {
        Ok(()) => (),
        Err(e) => {
            println!("Error: {}", e);
            process::exit(1);
        }
    };
}
