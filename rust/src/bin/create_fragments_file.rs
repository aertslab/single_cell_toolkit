use bstr::BString;
use clap::Parser;
use itertools::Itertools;

use noodles::bam;
use noodles::bgzf;

use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::Record;

use std::fs::File;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(author, version)]
#[command(name = "create_fragments_file")]
#[command(about = "Create fragments file from BAM file.")]
#[command(
    long_about = "Create BGZF-compressed fragments file from coordinate sorted BAM file:\n\
\u{20} - Filter reads:\n\
\u{20}   - read is properly paired.\n\
\u{20}   - read and its pair are located on the same chromosome.\n\
\u{20}   - read and its pair have a mapping quality of 30 or higher.\n\
\u{20}   - insert size is at least 10.\n\
\u{20}     This also ensures that only a fragment line is generated for\n\
\u{20}     the read of the read pair that is on the positive strand.\n\
\u{20}   - read is primary alignment.\n\
\u{20}   - read has an associated CB tag.\n\
\u{20} - Sort fragments with same chromosome and start by end and cell barcode.\n\
\u{20} - Collapse duplicate fragments (same coordinate and cell barcode).\n\
\u{20} - Write results to a BGZF-compressed fragments file."
)]
struct Cli {
    #[arg(
        short = 'i',
        long = "bam",
        required = true,
        help = "Input BAM file.",
        long_help = "Coordinate sorted input BAM file."
    )]
    input_bam_path: PathBuf,
    #[arg(
        short = 'o',
        long = "fragments",
        required = true,
        help = "Output fragments file.",
        long_help = "Output BGZF-compresssed fragments."
    )]
    output_fragments_path: PathBuf,
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Fragment {
    chromosome: BString,
    start: i64,
    end: i64,
    cb: BString,
}

// Write fragment with fragment counts fast to a BGZF compressed file.
fn write_fragment(
    fragments_file_writer: &mut bgzf::MultithreadedWriter,
    fragment: &Fragment,
    fragment_count: usize,
) {
    fragments_file_writer
        .write_all(&fragment.chromosome)
        .unwrap();
    fragments_file_writer.write_all(b"\t").unwrap();

    // Convert start, end of fragments and fragment count fast from integers
    // to strings with itoa.
    let mut itoa_buffer = itoa::Buffer::new();
    let fragment_start_str = itoa_buffer.format(fragment.start);
    fragments_file_writer
        .write_all(fragment_start_str.as_bytes())
        .unwrap();
    fragments_file_writer.write_all(b"\t").unwrap();

    let fragment_end_str = itoa_buffer.format(fragment.end);
    fragments_file_writer
        .write_all(fragment_end_str.as_bytes())
        .unwrap();
    fragments_file_writer.write_all(b"\t").unwrap();

    fragments_file_writer.write_all(&fragment.cb).unwrap();
    fragments_file_writer.write_all(b"\t").unwrap();

    let fragment_count_str = itoa_buffer.format(fragment_count);
    fragments_file_writer
        .write_all(fragment_count_str.as_bytes())
        .unwrap();
    fragments_file_writer.write_all(b"\n").unwrap();
}

fn create_fragments_file(
    bam_path: &Path,
    fragments_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    // Read input BAM file with 4 decompressing worker threads.
    let bam_decompressing_worker_threads = NonZeroUsize::new(4).unwrap();
    let bam_file = File::open(bam_path)?;

    let bam_bgzf_reader =
        bgzf::MultithreadedReader::with_worker_count(bam_decompressing_worker_threads, bam_file);

    // After decompressing BGZF blocks of BAM file in parallel,
    // decode the actual BAM records.
    let mut bam_reader = bam::io::Reader::from(bam_bgzf_reader);
    let header = bam_reader.read_header()?;

    // Write fragments to BGZF-compressed fragments file with 3 worker threads.
    // Set the compression level to 7, to match the default compression level (6)
    // of bgzip of HTSlib (which has libdeflate level 7 mapped to level 6:
    // https://github.com/samtools/htslib/pull/1488).
    let fragments_file = File::create(fragments_path)?;
    let compression_level = bgzf::writer::CompressionLevel::try_from(7)?;
    let fragments_compressing_worker_threads = NonZeroUsize::new(3).unwrap();

    let mut fragments_file_writer = bgzf::multithreaded_writer::Builder::default()
        .set_compression_level(compression_level)
        .set_worker_count(fragments_compressing_worker_threads)
        .build_from_writer(fragments_file);

    let mut fragments: Vec<Fragment> = Vec::new();

    let mut last_contig: Option<BString> = None;
    let mut last_start: Option<i64> = None;

    for result in bam_reader.records() {
        let record = result?;

        // Write a fragment line if and only if:
        //   - read is properly paired.
        //   - read and its pair are located on the same chromosome.
        //   - read and its pair have a mapping quality of 30 or higher.
        //   - insert size is at least 10.
        //     This also ensures that only a fragment line is generated for
        //     the read of the read pair that is on the positive strand.
        //   - read is primary alignment.
        //   - read has an associated CB tag.
        if record.flags().is_properly_segmented()
            && record.reference_sequence_id().unwrap().unwrap()
                == record.mate_reference_sequence_id().unwrap().unwrap()
            && record.mapping_quality() >= MappingQuality::new(30)
            && record.template_length() >= 10
            && !record.flags().is_secondary()
            && !record.flags().is_supplementary()
        {
            if let Some(mate_mapq) = record.data().get(&Tag::MATE_MAPPING_QUALITY).transpose()? {
                if mate_mapq.as_int() >= Some(30) {
                    if let Some(Value::String(cb)) =
                        record.data().get(&Tag::CELL_BARCODE_ID).transpose()?
                    {
                        // Get zero-based start position of fragment.
                        let start = record.alignment_start().transpose()?.unwrap().get() as i64 - 1;

                        // Construct fragment from current read pair.
                        let fragment = Fragment {
                            chromosome: record
                                .reference_sequence(&header)
                                .transpose()?
                                .map(|(name, _)| name)
                                .unwrap()
                                .to_owned(),
                            start: start + 4,
                            end: start + record.template_length() as i64 - 5,
                            cb: cb.to_owned(),
                        };

                        // If a new chromosome or new fragment start position is seen,
                        // sort the previous fragments with the same chromosome and
                        // start position and count the number of fragments per cell
                        // barcode with the same coordinates and write the result to
                        // BGZF-compressed fragments file.
                        if let Some(last_contig_unwrapped) = last_contig.as_deref() {
                            if last_contig_unwrapped != &fragment.chromosome {
                                // Current fragment is located on a different chromosome
                                // than the fragments collected in `fragments`.

                                fragments.sort();

                                // Write last fragments (all same start position) from
                                // previous chromosome.
                                for (fragment_count, fragment) in
                                    fragments.iter().dedup_with_count()
                                {
                                    write_fragment(
                                        &mut fragments_file_writer,
                                        fragment,
                                        fragment_count,
                                    );
                                }

                                // Start a new BGZF block for fragments from a new chromosome.
                                fragments_file_writer.flush().unwrap();

                                // Set new last chromosome and start position for the next iteration.
                                last_contig = Some(fragment.chromosome.clone());
                                last_start = Some(fragment.start);
                                fragments.clear();

                                fragments_file_writer.flush().unwrap();
                            } else if let Some(last_start_unwrapped) = last_start {
                                if last_start_unwrapped != fragment.start {
                                    // Current fragment is located on the same
                                    // chromosome than the fragments collected
                                    // in `fragments`, but the start coordinate
                                    // is different.

                                    fragments.sort();

                                    // Write previous collected fragments with same
                                    // chromosome and start position.
                                    for (fragment_count, fragment) in
                                        fragments.iter().dedup_with_count()
                                    {
                                        write_fragment(
                                            &mut fragments_file_writer,
                                            fragment,
                                            fragment_count,
                                        );
                                    }

                                    // Set new start position for the next iteration.
                                    last_start = Some(fragment.start);
                                    fragments.clear();
                                }
                            }
                        } else {
                            // Save chromosome and start position from first fragment
                            // for next iteration.
                            last_contig = Some(fragment.chromosome.clone());
                            last_start = Some(fragment.start);
                        }

                        // Push fragment to vector of fragments (which have all the
                        // same chromosome and start position).
                        fragments.push(fragment);
                    }
                }
            }
        }
    }

    // Write last collected fragments with same chromosome and start position.
    fragments.sort();

    for (fragment_count, fragment) in fragments.iter().dedup_with_count() {
        write_fragment(&mut fragments_file_writer, fragment, fragment_count);
    }

    fragments_file_writer.flush().unwrap();

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    create_fragments_file(&cli.input_bam_path, &cli.output_fragments_path)?;

    Ok(())
}
