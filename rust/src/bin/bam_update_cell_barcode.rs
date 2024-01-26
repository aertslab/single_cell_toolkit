use std::env;
use std::error::Error;
use std::path::Path;
use std::process;

use std::collections::HashMap;

use rust_htslib::bam::{record::Aux, Format, Header, Read, Reader, Writer};

// This lets us write `#[derive(Deserialize)]`.
use serde::Deserialize;

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
    input_bam_filename: &str,
    output_bam_filename: &str,
    old_bc_tag: &str,
    new_bc_tag: &str,
    barcode_mapping: &BarcodeMapping,
) {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
    let header = Header::from_template(input_bam.header());
    let mut output_bam = Writer::from_path(output_bam_filename, &header, Format::Bam).unwrap();

    input_bam
        .set_threads(2)
        .expect("Failed to set number of BAM reading threads to 2.");
    output_bam
        .set_threads(5)
        .expect("Failed to set number of BAM writing threads to 5.");
    let old_bc_tag = old_bc_tag.as_bytes();
    let new_bc_tag = new_bc_tag.as_bytes();
    let is_same_bc_tag = old_bc_tag == new_bc_tag;

    for r in input_bam.records() {
        let mut record = r.unwrap();

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
                    // Remove existing CB tag and value if the new one has the same tag name.
                    record.remove_aux(old_bc_tag).unwrap();
                } else if let Ok(_new_bc_aux) = record.aux(new_bc_tag) {
                    // Remove existing new CB tag if it was there already.
                    record.remove_aux(new_bc_tag).unwrap();
                }

                // Add new barcode tag and new barcode name to the current BAM record.
                record.push_aux(new_bc_tag, Aux::String(new_bc)).unwrap();
            }
        }

        // Write current BAM record to a new BAM file.
        output_bam.write(&record).unwrap();
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 6 {
        eprintln!(
            "Usage: bam_update_cell_barcode input.bam output.bam barcode_mapping.tsv old_bc_tag new_bc_tag"
        );
        process::exit(1);
    }

    // Input BAM file.
    let input_bam_filename = &args[1];
    // Output BAM file.
    let output_bam_filename = &args[2];
    // TSV file with 2 columns: old cell barcode and new cell barcode.
    let barcode_mapping_tsv_filename = &args[3];
    // SAM tag in BAM file that contains the old cell barcode: e.g. "CB".
    let old_bc_tag = &args[4];
    // SAM tag in BAM file to which the new cell barcode will be written:
    // e.g. "CB" (to overwrite the old one) or e.g. "de" (to write to a new tag and keep the old tag too).
    let new_bc_tag = &args[5];

    let barcode_mapping_tsv_path = Path::new(barcode_mapping_tsv_filename);

    let barcode_mapping = match read_barcode_mapping_tsv_file(barcode_mapping_tsv_path) {
        Ok(barcode_mapping) => barcode_mapping,
        Err(e) => {
            println!("{}", e);
            process::exit(1);
        }
    };

    bam_update_cell_barcode(
        input_bam_filename,
        output_bam_filename,
        old_bc_tag,
        new_bc_tag,
        &barcode_mapping,
    );
}
