use std::error::Error;
use std::path::{Path, PathBuf};
use std::process;

use clap::Parser;

use hashbrown::{HashMap, HashSet};

use rust_htslib::bam::{record::Aux, Format, Header, HeaderView, Read, Reader, Writer};
use rust_htslib::tpool::ThreadPool;

// This lets us write `#[derive(Deserialize)]`.
use serde::Deserialize;

#[derive(Parser, Debug)]
#[command(author, version, long_about = None)]
#[command(name = "create_parsebio_per_sample_bams")]
#[command(about = "Create per sample BAM files from ParseBio sublibrary BAM files.")]
struct Cli {
    #[arg(
        short = 'm',
        long = "metadata",
        required = true,
        help = "ParseBio cell metadata CSV file.",
        long_help = "ParseBio cell metadata CSV file."
    )]
    parsebio_cell_metadata_csv_path: PathBuf,
    #[arg(
        short = 's',
        long = "sublibrary_to_bam",
        required = true,
        help = "ParseBio sublibrary to BAM CSV file.",
        long_help = "ParseBio sublibrary to BAM CSV file."
    )]
    parsebio_sublibrary_to_bam_csv_path: PathBuf,
    #[arg(
        short = 'o',
        long = "output_prefix",
        required = true,
        help = "Output prefix.",
        long_help = "Output prefix used to create output BAM files for each ParseBio sublibrary."
    )]
    output_prefix: PathBuf,
}

// ParseBio cell metadata CSV record.
#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParseBioCellMetadataRecord {
    bc_wells: String,
    sample: String,
    species: String,
    gene_count: u32,
    tscp_count: u32,
    mread_count: u32,
    bc1_well: String,
    bc2_well: String,
    bc3_well: String,
    bc1_wind: u8,
    bc2_wind: u8,
    bc3_wind: u8,
}

// ParseBio sublibrary to BAM file CSV file record.
#[derive(Debug, Deserialize)]
struct SublibraryToBamRecord {
    sublibrary: String,
    bam_filename: String,
}

type BarcodeToSampleMapping = HashMap<String, String>;
type SampleToBamWriterMapping = HashMap<String, Writer>;
type SublibraryToBamMapping = HashMap<String, String>;

fn read_parsebio_cell_metadata_csv_file(
    parsebio_cell_metadata_csv_path: &Path,
) -> Result<BarcodeToSampleMapping, Box<dyn Error>> {
    let mut barcode_to_sample_mapping: BarcodeToSampleMapping = BarcodeToSampleMapping::new();

    // Build a CSV reader for a plain CSV file.
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b',')
        .escape(None)
        .double_quote(false)
        .quoting(false)
        .comment(Some(b'#'))
        .from_path(parsebio_cell_metadata_csv_path)?;

    for result in rdr.deserialize() {
        let parsebio_cell_metadata_record: ParseBioCellMetadataRecord = result?;

        barcode_to_sample_mapping
            .entry(parsebio_cell_metadata_record.bc_wells)
            .or_insert(parsebio_cell_metadata_record.sample.clone());
    }

    Ok(barcode_to_sample_mapping)
}

fn read_parsebio_sublibrary_to_bam_csv_file(
    parsebio_sublibrary_to_bam_csv_path: &Path,
) -> Result<SublibraryToBamMapping, Box<dyn Error>> {
    let mut sublibrary_to_bam_mapping: SublibraryToBamMapping = SublibraryToBamMapping::new();

    // Build a CSV reader for a plain CSV file.
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b',')
        .escape(None)
        .double_quote(false)
        .quoting(false)
        .comment(Some(b'#'))
        .from_path(parsebio_sublibrary_to_bam_csv_path)?;

    for result in rdr.deserialize() {
        let sublibrary_to_bam_record: SublibraryToBamRecord = result?;

        sublibrary_to_bam_mapping
            .entry(sublibrary_to_bam_record.sublibrary.clone())
            .or_insert(sublibrary_to_bam_record.bam_filename.clone());
    }

    Ok(sublibrary_to_bam_mapping)
}

fn fix_chromosome_names_in_bam_header(header: &Header) -> Vec<u8> {
    // Only keep "@HD" and "@SQ" lines from BAM header and fix chromosome names
    // in "@SQ" names by removing "hg38_".
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && (x.starts_with(b"@HD\t") || x.starts_with(b"@SQ\t")))
        .map(|x| match x.starts_with(b"@SQ\tSN:") {
            true => match x[7..].strip_prefix(b"hg38_") {
                Some(y) => {
                    let mut sq_line = b"@SQ\tSN:".to_vec();
                    sq_line.extend(Vec::<u8>::from(y));
                    sq_line
                }
                None => Vec::<u8>::from(x),
            },
            false => Vec::<u8>::from(x),
        })
        .collect::<Vec<Vec<u8>>>()
        .join(&b"\n"[..])
}

fn filter_bam_header(header: &Header, sublibrary: &str) -> Vec<u8> {
    // Only keep non-"@HD" and non-"@SQ" lines (@PG and @CO lines) from BAM header
    // and make "ID:STAR" field in "@PG" lines unique, by adding sublibrary name.
    header
        .to_bytes()
        .split(|x| x == &b'\n')
        .filter(|x| !x.is_empty() && !x.starts_with(b"@HD\t") && !x.starts_with(b"@SQ\t"))
        .map(Vec::<u8>::from)
        .map(|x| match x.starts_with(b"@PG\tID:STAR") {
            true => {
                let mut pg_line = Vec::<u8>::from(&x[..11]);
                pg_line.extend(format!(".{}", &sublibrary).into_bytes());
                pg_line.extend(Vec::<u8>::from(&x[11..]));
                pg_line
            }
            false => x,
        })
        .collect::<Vec<Vec<u8>>>()
        .join(&b"\n"[..])
}

// Read all ParseBio sublibrary BAM files and write per sample BAM files.
fn create_parsebio_per_sample_bams(
    barcode_to_sample_mapping: &BarcodeToSampleMapping,
    sublibrary_to_bam_mapping: &SublibraryToBamMapping,
    output_prefix: &Path,
    cmd_line_str: &str,
) -> Result<(), Box<dyn Error>> {
    let bam_thread_pool = ThreadPool::new(16)?;

    // Get all sublibrary names and sort them.
    let mut sublibraries: Vec<&String> = sublibrary_to_bam_mapping.keys().collect();
    sublibraries.sort_unstable();

    let mut merged_header = Vec::new();

    // Construct BAM header for each per sample BAM file by combining headers
    // of each ParseBio sublibrary.
    for (i, sublibrary) in sublibraries.iter().enumerate() {
        let bam = Reader::from_path(Path::new(
            sublibrary_to_bam_mapping.get(sublibrary.as_str()).unwrap(),
        ))?;

        // Read BAM header from current ParseBio sublibrary BAM file.
        let original_header = Header::from_template(bam.header());

        // Get "@HD" and "@SQ" lines if this is the first ParseBio sublibrary.
        if i == 0 {
            let fix_chromosome_names_header = fix_chromosome_names_in_bam_header(&original_header);
            merged_header.extend(fix_chromosome_names_header);
        }

        // Get all "@PG" and "@CO" lines.
        let filtered_header = filter_bam_header(&original_header, sublibrary);

        merged_header.extend(&b"\n"[..]);
        merged_header.extend(&filtered_header);
    }

    // Add "@PG" header line for "create_parsebio_samples_bam".
    let pg_header_line = format!(
        "\n@PG\tID:create_parsebio_samples_bam\tPN:create_parsebio_samples VN:{} CL:{}",
        env!("CARGO_PKG_VERSION"),
        &cmd_line_str
    )
    .into_bytes();
    merged_header.extend(&pg_header_line);

    let merged_header = Header::from_template(&HeaderView::from_bytes(&merged_header));

    let mut samples: HashSet<&String> = HashSet::new();
    for sample in barcode_to_sample_mapping.values() {
        samples.insert(sample);
    }

    let mut sample_to_bam_writer_mapping: SampleToBamWriterMapping =
        SampleToBamWriterMapping::new();

    // Initialize per sample BAM files with BAM header.
    for sample in samples {
        let mut output_bam_path = PathBuf::from(output_prefix);
        output_bam_path
            .as_mut_os_string()
            .push(format!("{}.bam", sample));

        let mut output_bam_writer =
            Writer::from_path(&output_bam_path, &merged_header, Format::Bam)?;

        output_bam_writer.set_thread_pool(&bam_thread_pool)?;
        sample_to_bam_writer_mapping
            .entry(sample.to_owned())
            .or_insert(output_bam_writer);
    }

    let cb_tag = b"CB";

    for sublibrary in sublibraries {
        let mut input_bam = Reader::from_path(Path::new(
            sublibrary_to_bam_mapping.get(sublibrary.as_str()).unwrap(),
        ))?;

        input_bam.set_thread_pool(&bam_thread_pool)?;

        for r in input_bam.records() {
            let mut record = r?;

            if let Ok(Aux::String(cb)) = record.aux(cb_tag) {
                let new_cb = format!("{}__{}", &cb, &sublibrary);

                // Write BAM record with updated full barcode name to per
                // sample BAM file, if the barcode was in the list of
                // filtered barcodes.
                if let Some(sample) = barcode_to_sample_mapping.get(&new_cb) {
                    // Update CB tag value with full barcode name.
                    record.remove_aux(cb_tag)?;
                    record.push_aux(cb_tag, Aux::String(&new_cb))?;

                    if let Some(bam_writer) = sample_to_bam_writer_mapping.get_mut(sample) {
                        // Write current BAM record to correct per sample BAM file.
                        bam_writer.write(&record)?;
                    }
                }
            }
        }
    }

    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let cmd_line_str = format!(
        "create_parsebio_samples_bam {} {} {}",
        &cli.parsebio_cell_metadata_csv_path.to_string_lossy(),
        &cli.parsebio_sublibrary_to_bam_csv_path.to_string_lossy(),
        &cli.output_prefix.to_string_lossy()
    );

    let barcode_to_sample_mapping =
        match read_parsebio_cell_metadata_csv_file(&cli.parsebio_cell_metadata_csv_path) {
            Ok(barcode_to_sample_mapping) => barcode_to_sample_mapping,
            Err(e) => {
                println!("Error: {}", e);
                process::exit(1);
            }
        };

    let sublibrary_to_bam_mapping =
        match read_parsebio_sublibrary_to_bam_csv_file(&cli.parsebio_sublibrary_to_bam_csv_path) {
            Ok(sublibrary_to_bam_mapping) => sublibrary_to_bam_mapping,
            Err(e) => {
                println!("Error: {}", e);
                process::exit(1);
            }
        };

    match create_parsebio_per_sample_bams(
        &barcode_to_sample_mapping,
        &sublibrary_to_bam_mapping,
        &cli.output_prefix,
        &cmd_line_str,
    ) {
        Ok(()) => (),
        Err(e) => {
            println!("Error: {}", e);
            process::exit(1);
        }
    };
}
