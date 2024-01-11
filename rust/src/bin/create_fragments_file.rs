use bio_types::genome::AbstractInterval;
use grep_cli;
use itertools::Itertools;
use rust_htslib::bam::{record::Aux, Read, Reader};
use termcolor::ColorChoice;

use itoa;
use std::env;
use std::io::Write;
use std::process;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Fragment {
    chromosome: String,
    start: i64,
    end: i64,
    cb: String,
}

// Write fragment with fragment counts fast to stdout.
fn write_fragment(
    stdout: &mut grep_cli::StandardStream,
    fragment: &Fragment,
    fragment_count: usize,
) {
    stdout.write(fragment.chromosome.as_bytes()).unwrap();
    stdout.write(b"\t").unwrap();

    // Convert start, end of fragments and fragment count fast from integers
    // to strings with itoa.
    let mut itoa_buffer = itoa::Buffer::new();
    let fragment_start_str = itoa_buffer.format(fragment.start);
    stdout.write(fragment_start_str.as_bytes()).unwrap();
    stdout.write(b"\t").unwrap();

    let fragment_end_str = itoa_buffer.format(fragment.end);
    stdout.write(fragment_end_str.as_bytes()).unwrap();
    stdout.write(b"\t").unwrap();

    stdout.write(fragment.cb.as_bytes()).unwrap();
    stdout.write(b"\t").unwrap();

    let fragment_count_str = itoa_buffer.format(fragment_count);
    stdout.write(fragment_count_str.as_bytes()).unwrap();
    stdout.write(b"\n").unwrap();
}

fn create_fragments_file(input_bam_filename: &str) {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();

    input_bam
        .set_threads(3)
        .expect("Failed to set number of BAM reading threads to 3.");

    let mut stdout = grep_cli::stdout(ColorChoice::Never);

    let mut fragments: Vec<Fragment> = Vec::new();

    let mut last_contig: Option<String> = None;
    let mut last_start: Option<i64> = None;
    let mut last_line_flush = 0u32;

    for r in input_bam.records() {
        let record = r.unwrap();

        // Write a fragment line if and only if:
        //   - read is properly paired.
        //   - read and its pair are located on the same chromosome.
        //   - read and its pair have a mapping quality of 30 or higher.
        //   - insert size is at least 10.
        //     This also ensures that only a fragment line is generated for
        //     the read of the read pair that is on the positive strand.
        //   - read is primary alignment.
        //   - read has an associated CB tag.
        if record.is_proper_pair()
            && record.tid() == record.mtid()
            && record.mapq() >= 30
            && record.insert_size() >= 10
            && !record.is_secondary()
            && !record.is_supplementary()
        {
            if let Ok(mate_mapq_aux) = record.aux(b"MQ") {
                // So far MQ tags in BAM files have been of I8 or U8 type.
                let mate_mapq = match mate_mapq_aux {
                    Aux::I8(mate_mapq) => mate_mapq as i32,
                    Aux::I16(mate_mapq) => mate_mapq as i32,
                    Aux::I32(mate_mapq) => mate_mapq as i32,
                    Aux::U8(mate_mapq) => mate_mapq as i32,
                    Aux::U16(mate_mapq) => mate_mapq as i32,
                    Aux::U32(mate_mapq) => mate_mapq as i32,
                    _ => -1, //bail!("Value for MQ tag is not an integer."),
                };

                if mate_mapq >= 30 {
                    if let Ok(Aux::String(cb)) = record.aux(b"CB") {
                        // Construct fragment from current read pair.
                        let fragment = Fragment {
                            chromosome: record.contig().to_owned(),
                            start: record.pos() + 4,
                            end: record.pos() + record.insert_size() - 5,
                            cb: cb.to_owned(),
                        };

                        // If a new chromosome or new fragment start position is seen,
                        // sort the previous fragments with the same chromosome and
                        // start position and count the number of fragments per cell
                        // barcode with the same coordinates and write the result to
                        // stdout.
                        if let Some(last_contig_unwrapped) = last_contig.as_deref() {
                            if last_contig_unwrapped != fragment.chromosome {
                                // Current fragment is located on a different chromosome
                                // than the fragments collected in `fragments`.

                                fragments.sort();

                                // Write last fragments (all same start position) from
                                // previous chromosome.
                                for (fragment_count, fragment) in
                                    fragments.iter().dedup_with_count()
                                {
                                    write_fragment(&mut stdout, fragment, fragment_count);
                                }

                                // Flush stdout and reset last line flush counter.
                                stdout.flush().unwrap();
                                last_line_flush = 0;


                                // Set new last chromosome and start position for the next iteration.
                                last_contig = Some(fragment.chromosome.clone());
                                last_start = Some(fragment.start);
                                fragments.clear();
                            } else {
                                if let Some(last_start_unwrapped) = last_start {
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
                                            write_fragment(&mut stdout, fragment, fragment_count);

                                            last_line_flush += 1;
                                        }

                                        // Flush stdout only if more than 50 lines were
                                        // printed as calling it after every new line
                                        // has a high overhead.
                                        if last_line_flush > 50 {
                                            stdout.flush().unwrap();
                                            last_line_flush = 0;
                                        }

                                        // Set new start position for the next iteration.
                                        last_start = Some(fragment.start);
                                        fragments.clear();
                                    }
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
        write_fragment(&mut stdout, fragment, fragment_count);
    }

    stdout.flush().unwrap();
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!(
            r#"
Purpose: Create fragments file from BAM file.

Usage:
    create_fragments_file sample.bam \
      | bgzip -@ 4 -c /dev/stdin \
      > sample.fragments.raw.tsv.gz

    - Create fragments file from coordinate sorted BAM file:
        - read is properly paired.
        - read and its pair are located on the same chromosome.
        - read and its pair have a mapping quality of 30 or higher.
        - insert size is at least 10.
          This also ensures that only a fragment line is generated for
          the read of the read pair that is on the positive strand.
        - read is primary alignment.
        - read has an associated CB tag.
    - Sort fragments with same chromosome and start by end and cell barcode.
    - Collapse duplicate fragments (same coordinate and cell barcode).
    - Write as bgzipped fragments file.
"#
        );
        process::exit(1);
    }

    let input_bam_filename = &args[1];

    create_fragments_file(input_bam_filename);
}
