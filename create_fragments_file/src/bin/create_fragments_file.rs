use bio_types::genome::AbstractInterval;
use grep_cli;
use rust_htslib::bam::{record::Aux, Read, Reader};
use termcolor::ColorChoice;

use itoa;
use std::env;
use std::io::Write;
use std::process;

fn create_fragments_file(input_bam_filename: &str) {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();

    input_bam
        .set_threads(3)
        .expect("Failed to set number of BAM reading threads to 3.");

    let mut stdout = grep_cli::stdout(ColorChoice::Never);

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
                if let Aux::I32(mate_mapq) = mate_mapq_aux {
                    if mate_mapq >= 30 {
                        if let Ok(Aux::String(cb)) = record.aux(b"CB") {
                            // Write output directly as bytes and do not invoke the formatting machinery of rust.
                            // The code in this block will produce the same output as:
                            //
                            // println!(
                            //     "{}\t{}\t{}\t{}",
                            //     record.contig(),
                            //     record.pos() + 4,
                            //     record.pos() + record.insert_size() - 5,
                            //     cb
                            // );

                            stdout.write(record.contig().as_bytes()).unwrap();
                            stdout.write(b"\t").unwrap();

                            // Convert start and end of fragments fast from integers to strings.
                            let mut pos_buffer = itoa::Buffer::new();
                            let pos_start = pos_buffer.format(record.pos() + 4);
                            stdout.write(pos_start.as_bytes()).unwrap();
                            stdout.write(b"\t").unwrap();

                            let pos_end =
                                pos_buffer.format(record.pos() + record.insert_size() - 5);
                            stdout.write(pos_end.as_bytes()).unwrap();
                            stdout.write(b"\t").unwrap();

                            stdout.write(cb.as_bytes()).unwrap();
                            stdout.write(b"\n").unwrap();
                        }
                    }
                }
            }
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: create_fragments_file input.bam");
        process::exit(1);
    }

    let input_bam_filename = &args[1];

    create_fragments_file(input_bam_filename);
}
