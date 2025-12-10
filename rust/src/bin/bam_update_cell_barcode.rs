use std::error::Error;
use std::path::{Path, PathBuf};
use std::process;

use clap::Parser;

use hashbrown::HashMap;

use rust_htslib::bam::{record::Aux, Format, Header, Read, Reader, Writer};

// This lets us write `#[derive(Deserialize)]`.
use serde::Deserialize;

#[derive(Parser, Debug)]
#[command(author, version, long_about = None)]
#[command(name = "bam_update_cell_barcode")]
#[command(about = "Update cell barcodes in a BAM file.")]
struct Cli {
    #[arg(
        short = 'i',
        required = true,
        help = "Input BAM file.",
        long_help = "Input BAM file."
    )]
    input_bam_path: PathBuf,
    #[arg(
        short = 'o',
        required = true,
        help = "Output BAM file.",
        long_help = "Output BAM file in which old cell barcodes are replaced with new cell barcodes."
    )]
    output_bam_path: PathBuf,
    #[arg(
        short = 'b',
        required = true,
        help = "Cell barcode mapping TSV file.",
        long_help = "Cell barcode mapping TSV file consisting of 2 columns:\n\
        \u{20} 1) old cell barcodes\n\
        \u{20} 2) new cell barcodes"
    )]
    barcode_mapping_tsv_path: PathBuf,
    #[arg(
        short='c',
        default_value_t = String::from("CB"),
        help="SAM tag used for old cell barcodes.",
        long_help="SAM tag used for old cell barcodes in input BAM file.",
    )]
    old_bc_tag: String,
    #[arg(
        short='d',
        default_value_t = String::from("CB"),
        help="SAM tag to use for new cell barcodes.",
        long_help="SAM tag to use for writing new cell barcodes to output BAM file.\n\
        If the new SAM tag is the same as the old tag, the old SAM tag value will be overwritten.",
    )]
    new_bc_tag: String,
}

// We don't need to derive `Debug` (which doesn't require Serde), but it's a
// good habit to do it for all your types.
#[derive(Debug, Deserialize)]
struct BarcodeRecord {
    old_barcode: String,
    new_barcode: String,
}

type BarcodeMapping = HashMap<String, String>;

fn read_barcode_mapping_tsv_file(
    barcode_mapping_tsv_path: &Path,
) -> Result<BarcodeMapping, Box<dyn Error>> {
    let mut barcode_mapping: BarcodeMapping = BarcodeMapping::new();

    // Build a CSV reader for a plain TSV file.
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .escape(None)
        .double_quote(false)
        .quoting(false)
        .comment(Some(b'#'))
        .from_path(barcode_mapping_tsv_path)?;

    for result in rdr.deserialize() {
        let barcode_record: BarcodeRecord = result?;

        barcode_mapping.insert(barcode_record.old_barcode, barcode_record.new_barcode);
    }

    Ok(barcode_mapping)
}

fn bam_update_cell_barcode(
    input_bam_path: &Path,
    output_bam_path: &Path,
    old_bc_tag: &str,
    new_bc_tag: &str,
    barcode_mapping: &BarcodeMapping,
) -> Result<(), Box<dyn Error>> {
    let mut input_bam = Reader::from_path(input_bam_path)?;
    let header = Header::from_template(input_bam.header());
    let mut output_bam = Writer::from_path(output_bam_path, &header, Format::Bam)?;

    input_bam.set_threads(2)?;
    output_bam.set_threads(5)?;

    let old_bc_tag = old_bc_tag.as_bytes();
    let new_bc_tag = new_bc_tag.as_bytes();
    let is_same_bc_tag = old_bc_tag == new_bc_tag;

    for r in input_bam.records() {
        let mut record = r?;

        // Check if we have an old barcode tag for the current BAM record.
        if let Ok(old_bc_aux) = record.aux(old_bc_tag) {
            // Get old barcode and check if it is in the barcode mapping file and get the new barcode.
            let new_bc = if let Aux::String(old_bc) = &old_bc_aux {
                barcode_mapping.get(old_bc as &str)
            } else {
                None
            };

            if let Some(new_bc) = new_bc {
                if is_same_bc_tag {
                    // Update existing CB tag and value if the new one has the same tag name.
                    record.update_aux(new_bc_tag, Aux::String(&new_bc))?;
                } else if let Ok(_new_bc_aux) = record.aux(new_bc_tag) {
                    // Remove existing new CB tag if it was there already.
                    record.remove_aux(new_bc_tag)?;

                    // Add new barcode tag and new barcode name to the current BAM record.
                    record.push_aux(new_bc_tag, Aux::String(&new_bc))?;
                } else {
                    // Add new barcode tag and new barcode name to the current BAM record.
                    record.push_aux(new_bc_tag, Aux::String(&new_bc))?;
                }
            }

        }

        // Write current BAM record to a new BAM file.
        output_bam.write(&record)?;
    }

    Ok(())
}

fn main() {
    let cli = Cli::parse();

    // Read barcode mapping TSV file.
    let barcode_mapping = match read_barcode_mapping_tsv_file(&cli.barcode_mapping_tsv_path) {
        Ok(barcode_mapping) => barcode_mapping,
        Err(e) => {
            eprintln!(
                "Error: \"{}\": {}",
                &cli.barcode_mapping_tsv_path.to_string_lossy(),
                e
            );
            process::exit(1);
        }
    };

    // Update cell barcodes in a BAM file.
    match bam_update_cell_barcode(
        &cli.input_bam_path,
        &cli.output_bam_path,
        &cli.old_bc_tag,
        &cli.new_bc_tag,
        &barcode_mapping,
    ) {
        Ok(()) => (),
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    };
}
