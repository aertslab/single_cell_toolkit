use std::collections::BTreeMap;
use std::fmt;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{bail, Context, Result};
use clap::Parser;
use hashbrown::{HashMap, HashSet};

use rayon::prelude::*;

use itertools::Itertools;
use num_iter::range_step_inclusive;

use rust_htslib::bam::{
    record::Aux, CompressionLevel, Format, Header, HeaderView, IndexedReader, Read, Record, Writer,
};

// This lets us write `#[derive(Deserialize)]`.
use serde::Deserialize;

#[derive(Parser, Debug)]
#[command(author, version, long_about = None)]
#[command(name = "split_bams_per_cluster_htslib")]
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
        \u{20} 2) bam_filename: BAM filename\n\
        with \"sample\\tbam_filename\\n\" as header."
    )]
    sample_to_bam_tsv_path: PathBuf,
    #[arg(
        short = 'c',
        long = "cluster_cb_sample",
        required = true,
        help = "Cluster to original cell barcode, new cell barcode and sample name mapping TSV file.",
        long_help = "Cluster to original cell barcode, new cell barcode and sample name mapping TSV file consisting of 4 columns:\n\
        \u{20} 1) cluster: cluster\n\
        \u{20} 2) cell_barcode_input: input cell barcode (as written in input BAM files)\n\
        \u{20} 3) cell_barcode_output: output cell barcode (as to be written to output cluster BAM file)\n\
        \u{20} 4) sample: sample name (same name as used in sample to BAM filename mapping TSV file)\n\
        with \"cluster\\tcell_barcode_input\\tcell_barcode_output\\tsample\\n\" as header."
    )]
    cluster_to_cb_and_sample_tsv_path: PathBuf,
    #[arg(
        short = 'o',
        long = "output_prefix",
        required = true,
        help = "Output prefix.",
        long_help = "Output prefix used to create output cluster BAM files for each cluster.\n\
        \u{20} e.g. \"./\", \"/full/path/\", \"./my_sample.\", \"/full/path/my_sample__\", ..."
    )]
    output_prefix: PathBuf,
    #[arg(
        long = "chroms",
        num_args(1..),
        required = false,
        value_delimiter = ',',
        help = "List of chromosome names to keep reads for in the output BAM files.",
        long_help = "List of chromosome names to keep reads for in the output BAM files.\n\
        If not specified, keep reads for all chromosomes in the output BAM files.\n\
        \u{20} e.g. --chroms 'chr1,chr3,chr8' or --chroms chr1 chr5 chr8"
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
        short = 'd',
        long = "tag",
        required = false,
        default_value = "CB",
        help = "SAM tag to use to look for the cell barcodes in the BAM files.",
        long_help = "SAM tag to use to look for the cell barcodes in the BAM files."
    )]
    cb_tag: CellBarcodeTag,
    #[arg(
        short = 't',
        long = "threads",
        required = false,
        default_value_t = 16,
        help = "Number of threads to use.",
        long_help = "Number of threads to use for reading input BAM files and writing cluster BAM files.\n\
        Set to 0 to use the number of logical CPUs available on the machine."
    )]
    threads: usize,
    #[arg(
        short = 'C',
        long = "chunk_size",
        required = false,
        default_value_t = 1_000_000,
        help = "Fetch reads from each BAM file in chunks of X bp.",
        long_help = "Fetch reads from each BAM file in chunks of X bp.\n\
        Reduce this value if split_bams_per_cluster_htslib uses too much memory or\n\
        increase this value if split_bams_per_cluster_htslib reads very sparse BAM files."
    )]
    chunk_size: u64,
}

/// A row in the sample-to-BAM-filename TSV file.
///
/// Associates a sample identifier with its input BAM filename.
#[derive(Debug, Deserialize)]
struct SampleBamFilenameRecord {
    sample: Sample,
    bam_filename: BamFilename,
}

/// A row in the cluster-to-cell-barcode-and-sample TSV file.
///
/// Associates a cluster and sample with an input cell barcode from an input
/// sample BAM file and an output cell barcode written to an output cluster BAM file.
#[derive(Debug, Deserialize)]
struct ClusterCbSampleRecord {
    cluster: Cluster,
    cell_barcode_input: CellBarcodeInput,
    cell_barcode_output: CellBarcodeOutput,
    sample: Sample,
}

#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct BamFilename(String);

impl BamFilename {
    fn as_str(&self) -> &str {
        &self.0
    }

    fn as_path(&self) -> &Path {
        Path::new(self.as_str())
    }
}

#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct Sample(String);

impl Sample {
    fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct Cluster(String);

impl Cluster {
    fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct CellBarcodeInput(String);

impl CellBarcodeInput {
    fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct CellBarcodeOutput(String);

impl CellBarcodeOutput {
    fn as_str(&self) -> &str {
        &self.0
    }
}

macro_rules! impl_display_for_string_newtype {
    ($($type:ty),+ $(,)?) => {
        $(
            impl fmt::Display for $type {
                fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                    self.as_str().fmt(formatter)
                }
            }
        )+
    };
}

impl_display_for_string_newtype!(
    BamFilename,
    Sample,
    Cluster,
    CellBarcodeInput,
    CellBarcodeOutput,
);

#[derive(Clone, Debug, Default)]
// A BTreeMap ensures deterministic per-cluster BAM output by reading BAM files
// in key order for every chunk. Records with the same start position therefore
// retain a consistent order across input BAM files.
struct BamToSampleBTreeMapping {
    bam_to_sample: BTreeMap<BamFilename, Sample>,
}

#[derive(Clone, Debug, Default)]
struct SampleSet {
    samples: HashSet<Sample>,
}

#[derive(Clone, Debug, Default)]
struct ClusterToSamplesMapping {
    cluster_to_samples: BTreeMap<Cluster, SampleSet>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct ClusterId {
    id: usize,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct CellBarcodeOutputAndClusterId {
    cell_barcode_output: CellBarcodeOutput,
    cluster_id: ClusterId,
}

#[derive(Clone, Debug, Default)]
struct CellBarcodeInputToCellBarcodeOutputAndClusterIdMapping {
    cb_input_to_cb_output_and_cluster_id: HashMap<CellBarcodeInput, CellBarcodeOutputAndClusterId>,
}

#[derive(Clone, Debug, Default)]
struct SampleToCellBarcodeInputToCellBarcodeOutputAndClusterIdMapping {
    sample_to_cb_and_cluster_id_mapping:
        HashMap<Sample, CellBarcodeInputToCellBarcodeOutputAndClusterIdMapping>,
}

#[derive(Default)]
struct BamFileToBamIndexedReaderMapping {
    bam_to_reader: HashMap<BamFilename, IndexedReader>,
}

struct BamInput {
    sample: Sample,
    reader: IndexedReader,
    per_cluster_records: Vec<Vec<Record>>,
}

struct ClusterOutput {
    cluster: Cluster,
    writer: Writer,
    records: Vec<Record>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct CellBarcodeTag([u8; 2]);

impl CellBarcodeTag {
    fn as_bytes(&self) -> &[u8; 2] {
        &self.0
    }
}

impl Default for CellBarcodeTag {
    fn default() -> Self {
        Self(*b"CB")
    }
}

impl FromStr for CellBarcodeTag {
    type Err = &'static str;

    fn from_str(tag: &str) -> std::result::Result<Self, Self::Err> {
        tag.as_bytes()
            .try_into()
            .map(Self)
            .map_err(|_| "SAM tag for barcode should be exactly 2 characters long.")
    }
}

trait ContainsSlice<T>: PartialEq<[T]> {
    fn contains_slice(&'_ self, slice: &'_ [T]) -> bool;
}

impl<T, Item: PartialEq<T>> ContainsSlice<T> for [Item] {
    fn contains_slice(self: &'_ [Item], slice: &'_ [T]) -> bool {
        let len = slice.len();
        if len == 0 {
            return true;
        }
        self.windows(len).any(move |sub_slice| sub_slice == slice)
    }
}

fn read_sample_to_bam_tsv_file(sample_to_bam_tsv_path: &Path) -> Result<BamToSampleBTreeMapping> {
    let mut bam_to_sample_mapping = BamToSampleBTreeMapping::default();

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
        let sample_bam_filename_record: SampleBamFilenameRecord = result?;

        // Add sample to BAM filename mapping.
        bam_to_sample_mapping
            .bam_to_sample
            .entry(sample_bam_filename_record.bam_filename)
            .or_insert(sample_bam_filename_record.sample);
    }

    Ok(bam_to_sample_mapping)
}

fn read_cluster_to_cb_and_sample_tsv_file(
    cluster_to_cb_and_sample_tsv_path: &Path,
) -> Result<(
    ClusterToSamplesMapping,
    SampleToCellBarcodeInputToCellBarcodeOutputAndClusterIdMapping,
)> {
    let mut cluster_to_samples_mapping = ClusterToSamplesMapping::default();

    let mut sample_to_cb_input_to_cb_output_and_cluster_id_mapping =
        SampleToCellBarcodeInputToCellBarcodeOutputAndClusterIdMapping::default();

    // Build a CSV reader for a plain TSV file.
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .escape(None)
        .double_quote(false)
        .quoting(false)
        .comment(Some(b'#'))
        .from_path(cluster_to_cb_and_sample_tsv_path)?;

    let cluster_cb_sample_records = rdr
        .deserialize()
        .collect::<std::result::Result<Vec<ClusterCbSampleRecord>, csv::Error>>()?;

    let mut clusters = cluster_cb_sample_records
        .iter()
        .map(|record| record.cluster.clone())
        .collect::<Vec<_>>();
    clusters.sort_unstable();
    clusters.dedup();

    let cluster_to_id_mapping = clusters
        .into_iter()
        .enumerate()
        .map(|(id, cluster)| (cluster, ClusterId { id }))
        .collect::<BTreeMap<_, _>>();

    for cluster_cb_sample_record in cluster_cb_sample_records {
        let cluster_id = cluster_to_id_mapping[&cluster_cb_sample_record.cluster];

        // Create cluster to samples mapping.
        cluster_to_samples_mapping
            .cluster_to_samples
            .entry(cluster_cb_sample_record.cluster.clone())
            .or_default()
            .samples
            .insert(cluster_cb_sample_record.sample.clone());

        // Create nested hashmap with multiple levels:
        //   - level 1: Sample name to cell barcode input mapping
        //   - level 2: Cell barcode input to cell barcode output mapping
        //   - level 3: Cell barcode output to cluster mapping
        sample_to_cb_input_to_cb_output_and_cluster_id_mapping
            .sample_to_cb_and_cluster_id_mapping
            .entry(cluster_cb_sample_record.sample)
            .or_default()
            .cb_input_to_cb_output_and_cluster_id
            .entry(cluster_cb_sample_record.cell_barcode_input)
            .or_insert(CellBarcodeOutputAndClusterId {
                cell_barcode_output: cluster_cb_sample_record.cell_barcode_output,
                cluster_id,
            });
    }

    Ok((
        cluster_to_samples_mapping,
        sample_to_cb_input_to_cb_output_and_cluster_id_mapping,
    ))
}

/// Check if "@HD" line in BAM header contains "SO:coordinate".
fn has_coordinate_sorted_bam_header(header: &Header) -> bool {
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .any(|x| !x.is_empty() && (x.starts_with(b"@HD\t") && x.contains_slice(b"\tSO:coordinate")))
}

/// Returns the `@HD` and `@SQ` records from a BAM header.
///
/// The retained records are joined with newline separators in their original order.
fn get_hd_and_sq_bam_header_lines(header: &Header) -> Vec<u8> {
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && (x.starts_with(b"@HD\t") || x.starts_with(b"@SQ\t")))
        .fold(Vec::new(), |mut header_lines, line| {
            if !header_lines.is_empty() {
                header_lines.push(b'\n');
            }
            header_lines.extend_from_slice(line);
            header_lines
        })
}

/// Returns the `@SQ` records from a BAM header.
///
/// The retained records are joined with newline separators in their original order.
fn get_sq_bam_header_lines(header: &Header) -> Vec<u8> {
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && (x.starts_with(b"@SQ\t")))
        .fold(Vec::new(), |mut header_lines, line| {
            if !header_lines.is_empty() {
                header_lines.push(b'\n');
            }
            header_lines.extend_from_slice(line);
            header_lines
        })
}

/// Appends `.{sample}` to `ID:` and `PP:` fields in an `@PG` header line.
///
/// This makes program identifiers unique when headers from multiple sample BAM
/// files are combined while preserving all other fields.
fn suffix_pg_id_and_pp_fields(pg_line: &[u8], sample: &str) -> Vec<u8> {
    pg_line
        .split(|byte| *byte == b'\t')
        .fold(Vec::new(), |mut suffixed_pg_line, field| {
            if !suffixed_pg_line.is_empty() {
                suffixed_pg_line.push(b'\t');
            }

            if field.starts_with(b"ID:") || field.starts_with(b"PP:") {
                suffixed_pg_line.extend_from_slice(field);
                suffixed_pg_line.push(b'.');
                suffixed_pg_line.extend(sample.as_bytes());
            } else {
                suffixed_pg_line.extend_from_slice(field);
            }

            suffixed_pg_line
        })
}

/// Returns non-`@HD` and non-`@SQ` BAM header lines for a sample.
///
/// `@PG` lines have their `ID:` and `PP:` fields suffixed with the sample name
/// so program identifiers remain unique when multiple sample headers are combined.
fn get_non_hd_sq_and_fix_pg_bam_header_lines(header: &Header, sample: &str) -> Vec<u8> {
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && !x.starts_with(b"@HD\t") && !x.starts_with(b"@SQ\t"))
        .fold(Vec::new(), |mut header_lines, line| {
            if !header_lines.is_empty() {
                header_lines.push(b'\n');
            }

            match line.starts_with(b"@PG\t") {
                true => header_lines.extend(suffix_pg_id_and_pp_fields(line, sample)),
                false => header_lines.extend_from_slice(line),
            }

            header_lines
        })
}

/// Returns whether a read is part of a proper scATAC fragment pair.
///
/// A fragment read must be properly paired, have its mate on the same chromosome,
/// have a mapping quality of at least 30, have an absolute insert size of at least
/// 10, and be a primary alignment. Unless `ignore_mate_mapping_quality` is set,
/// the mate mapping quality (`MQ`) tag must also be at least 30.
fn is_fragment_read(record: &Record, ignore_mate_mapping_quality: bool) -> bool {
    // Only keep reads that will be used to create scATAC-seq fragments.
    //   - read is properly paired.
    //   - read and its pair are located on the same chromosome.
    //   - read and its pair have a mapping quality of 30 or higher.
    //   - insert size is at least 10 in absolute value.
    //   - read is primary alignment.
    if !record.is_proper_pair()
        || record.tid() != record.mtid()
        || record.mapq() < 30
        || record.insert_size().abs() < 10
        || record.is_secondary()
        || record.is_supplementary()
    {
        return false;
    }

    if ignore_mate_mapping_quality {
        return true;
    }

    let Ok(mate_mapq_aux) = record.aux(b"MQ") else {
        return false;
    };

    // So far MQ tags in BAM files have been of I8 or U8 type.
    let mate_mapq = match mate_mapq_aux {
        Aux::I8(mate_mapq) => mate_mapq as i32,
        Aux::I16(mate_mapq) => mate_mapq as i32,
        Aux::I32(mate_mapq) => mate_mapq,
        Aux::U8(mate_mapq) => mate_mapq as i32,
        Aux::U16(mate_mapq) => mate_mapq as i32,
        Aux::U32(mate_mapq) => mate_mapq as i32,
        _ => -1,
    };

    mate_mapq >= 30
}

fn split_bams_per_cluster(
    bam_to_sample_mapping: &BamToSampleBTreeMapping,
    cluster_to_samples_mapping: &ClusterToSamplesMapping,
    sample_to_cb_input_to_cb_output_and_cluster_id_mapping: &SampleToCellBarcodeInputToCellBarcodeOutputAndClusterIdMapping,
    output_prefix: &Path,
    chromosomes: &Option<Vec<String>>,
    fragment_reads_only: bool,
    ignore_mate_mapping_quality: bool,
    cb_tag: &CellBarcodeTag,
    threads: usize,
    chunk_size: u64,
    cmd_line_str: &str,
) -> Result<()> {
    let rayon_thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("failed to create Rayon thread pool")?;

    let mut bam_file_to_bam_indexed_reader_mapping = BamFileToBamIndexedReaderMapping::default();
    let mut cluster_outputs =
        Vec::with_capacity(cluster_to_samples_mapping.cluster_to_samples.len());

    let mut merged_header_view: Option<HeaderView> = None;

    // Construct BAM header for each per cluster BAM file by combining headers
    // of each per sample BAM file.
    for (cluster, cluster_samples) in &cluster_to_samples_mapping.cluster_to_samples {
        // Get all BAM filenames for the cluster in BTreeMap key order.
        let bam_filenames = bam_to_sample_mapping
            .bam_to_sample
            .iter()
            .filter_map(|(bam_filename, sample)| {
                cluster_samples
                    .samples
                    .contains(sample)
                    .then_some(bam_filename)
            })
            .collect::<Vec<_>>();

        // Create merged BAM header per cluster BAM file.
        let mut merged_header = Vec::new();
        let mut hd_and_sq_bam_header_lines: Option<Vec<u8>> = None;
        let mut sq_bam_header_lines: Option<Vec<u8>> = None;

        bam_filenames
            .iter()
            .enumerate()
            .try_for_each(|(i, bam_filename)| {
                let bam = bam_file_to_bam_indexed_reader_mapping
                    .bam_to_reader
                    .entry((*bam_filename).clone())
                    .or_insert(IndexedReader::from_path(bam_filename.as_path())?);

                // Read BAM header from current sample BAM file.
                let original_header = Header::from_template(bam.header());

                // Get "@HD" and "@SQ" lines and check if they match exactly with those
                // lines in the first sample BAM file.
                match i {
                    0 => {
                        if !has_coordinate_sorted_bam_header(&original_header) {
                            bail!(
                                "BAM file \"{}\" for cluster \"{}\" is not coordinate sorted.",
                                bam_filename.as_str(),
                                cluster.as_str()
                            );
                        }

                        hd_and_sq_bam_header_lines =
                            Some(get_hd_and_sq_bam_header_lines(&original_header));
                        merged_header.extend(
                            hd_and_sq_bam_header_lines
                                .as_ref()
                                .context("missing @HD/@SQ header lines")?,
                        );

                        sq_bam_header_lines = Some(get_sq_bam_header_lines(&original_header));
                    }
                    _ => {
                        if !has_coordinate_sorted_bam_header(&original_header) {
                            bail!(
                                "BAM file \"{}\" for cluster \"{}\" is not coordinate sorted.",
                                bam_filename.as_str(),
                                cluster.as_str()
                            );
                        }

                        if get_sq_bam_header_lines(&original_header)
                            != *sq_bam_header_lines
                                .as_ref()
                                .context("missing @SQ header lines")?
                        {
                            bail!(
                            "BAM file \"{}\" for cluster \"{}\" has different chromosome order.",
                            bam_filename.as_str(),
                            cluster.as_str()
                        );
                        }
                    }
                }

                // Get all "@PG", "@CO" and "@RG" lines.
                let non_hd_sq_and_fix_pg_bam_header_lines =
                    get_non_hd_sq_and_fix_pg_bam_header_lines(
                        &original_header,
                        bam_to_sample_mapping
                            .bam_to_sample
                            .get(*bam_filename)
                            .unwrap()
                            .as_str(),
                    );

                merged_header.extend(&b"\n"[..]);
                merged_header.extend(&non_hd_sq_and_fix_pg_bam_header_lines);

                Ok::<(), anyhow::Error>(())
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
            .push(format!("{}.bam", cluster.as_str()));

        let mut cluster_bam_writer =
            Writer::from_path(&cluster_bam_path, &merged_header, Format::Bam)?;

        cluster_bam_writer.set_compression_level(CompressionLevel::Fastest)?;
        cluster_outputs.push(ClusterOutput {
            cluster: cluster.clone(),
            writer: cluster_bam_writer,
            records: Vec::new(),
        });

        // Store merged header view for first cluster for later use to get chromosome names and lengths.
        if merged_header_view.is_none() {
            merged_header_view = Some(HeaderView::from_header(&merged_header));
        }
    }

    let merged_header_view = merged_header_view.unwrap();

    // Move initialized readers into BAM inputs with reusable per-cluster record buffers.
    let mut bam_inputs = bam_to_sample_mapping
        .bam_to_sample
        .iter()
        .filter_map(|(bam_filename, sample)| {
            bam_file_to_bam_indexed_reader_mapping
                .bam_to_reader
                .remove(bam_filename)
                .map(|reader| {
                    BamInput {
                        sample: sample.clone(),
                        reader,
                        per_cluster_records: (0..cluster_outputs.len())
                            .map(|_| Vec::new())
                            .collect(),
                    }
                })
        })
        .collect::<Vec<_>>();

    // Loop over each chromosome and fetch reads in chunks from each BAM file and sort
    // them by position before writing them to the per cluster BAM file.
    rayon_thread_pool.install(|| -> Result<()> {
        for tid in 0..merged_header_view.target_count() {
            let chrom_name = std::str::from_utf8(merged_header_view.tid2name(tid))
                .context("Chromosome name is not valid UTF-8.")?
                .to_string();
            if chromosomes.is_some() && !chromosomes.as_ref().unwrap().contains(&chrom_name) {
                continue;
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

                println!(
                    "Read reads for {}:{}-{} from input BAM files...",
                    chrom_name, start, end
                );

                bam_inputs.par_iter_mut().try_for_each(|bam_input| {
                    let indexed_input_bam = &mut bam_input.reader;
                    // Fetch chunk from current BAM file (coordinates are 0-based, and end is exclusive).
                    indexed_input_bam.fetch((tid, start, end))?;

                    let cb_input_to_cb_output_and_cluster_id_mapping =
                        sample_to_cb_input_to_cb_output_and_cluster_id_mapping
                            .sample_to_cb_and_cluster_id_mapping
                            .get(&bam_input.sample)
                            .unwrap();

                    for per_cluster_records in &mut bam_input.per_cluster_records {
                        per_cluster_records.clear();
                    }

                    // Filter reads of current chunk and collect them by output cluster.
                    for r in indexed_input_bam.records() {
                        let mut record = r?;

                        if record.pos() < start || record.pos() > end || record.tid() != tid as i32
                        {
                            // Skip reads that are not in the current region, which might
                            // be pulled in when using fetch as this might result in
                            // duplicated reads in the output.
                            continue;
                        }

                        if fragment_reads_only
                            && !is_fragment_read(&record, ignore_mate_mapping_quality)
                        {
                            // Skip reads that are not from a proper scATAC fragment pair.
                            continue;
                        }

                        if let Ok(Aux::String(cb)) = record.aux(cb_tag.as_bytes()) {
                            // Add BAM record with updated full barcode name to per
                            // cluster records, if the barcode was in the list of
                            // filtered barcodes.
                            if let Some(cb_output_and_cluster_id) =
                                cb_input_to_cb_output_and_cluster_id_mapping
                                    .cb_input_to_cb_output_and_cluster_id
                                    .get(&CellBarcodeInput(cb.to_owned()))
                            {
                                let CellBarcodeOutputAndClusterId {
                                    cell_barcode_output: cb_output,
                                    cluster_id,
                                } = cb_output_and_cluster_id;

                                // Update CB tag value with full barcode name.
                                record.update_aux(
                                    cb_tag.as_bytes(),
                                    Aux::String(cb_output.as_str()),
                                )?;

                                bam_input.per_cluster_records[cluster_id.id].push(record);
                            }
                        }
                    }

                    Ok::<(), rust_htslib::errors::Error>(())
                })?;

                // Merge each BAM input's per-cluster buffers into the corresponding
                // aggregate output buffers. `append` leaves the per-BAM buffers empty
                // while retaining their capacity for the next chunk.
                for bam_input in &mut bam_inputs {
                    for (cluster_output, bam_cluster_records) in cluster_outputs
                        .iter_mut()
                        .zip(&mut bam_input.per_cluster_records)
                    {
                        cluster_output.records.append(bam_cluster_records);
                    }
                }

                println!(
                    "Write reads for {}:{}-{} to per cluster BAM files...",
                    chrom_name, start, end
                );

                // Sort reads by position for each cluster and write them to per cluster BAM files.
                cluster_outputs
                    .par_iter_mut()
                    .try_for_each(|cluster_output| {
                        cluster_output
                            .records
                            .sort_by_key(|record| (record.tid(), record.pos()));

                        for record in &cluster_output.records {
                            cluster_output.writer.write(record)?;
                        }

                        cluster_output.records.clear();

                        Ok::<(), rust_htslib::errors::Error>(())
                    })?;
            }
        }

        Ok(())
    })
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Get full command line so it can be added later to the PG line in the BAM header.
    let cmd_line_str = std::env::args().collect::<Vec<_>>().join(" ");

    let bam_to_sample_mapping = read_sample_to_bam_tsv_file(&cli.sample_to_bam_tsv_path)
        .with_context(|| {
            format!(
                "failed to read sample to BAM filename mapping TSV file \"{}\"",
                cli.sample_to_bam_tsv_path.display()
            )
        })?;

    let (cluster_to_samples_mapping, sample_to_cb_input_to_cb_output_and_cluster_id_mapping) =
        read_cluster_to_cb_and_sample_tsv_file(&cli.cluster_to_cb_and_sample_tsv_path)
            .with_context(|| {
                format!(
                    "failed to read cluster to cell barcode and sample mapping TSV file \"{}\"",
                    cli.cluster_to_cb_and_sample_tsv_path.display()
                )
            })?;

    split_bams_per_cluster(
        &bam_to_sample_mapping,
        &cluster_to_samples_mapping,
        &sample_to_cb_input_to_cb_output_and_cluster_id_mapping,
        &cli.output_prefix,
        &cli.chromosomes,
        cli.fragment_reads_only,
        cli.ignore_mate_mapping_quality,
        &cli.cb_tag,
        cli.threads,
        cli.chunk_size,
        &cmd_line_str,
    )
    .context("failed to split BAM files per cluster")?;

    Ok(())
}
