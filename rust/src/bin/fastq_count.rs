use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::env;
use std::path::Path;
use std::path::PathBuf;

fn count_nbr_reads(fastq_filename: &Path) -> u64 {
    let mut nbr_reads: u64 = 0;

    let mut fastq_reader = parse_fastx_file(&fastq_filename)
        .unwrap_or_else(|_| panic!("Invalid FASTQ file \"{}\"", fastq_filename.display()));

    while let Some(_fastq_record) = fastq_reader.next() {
        nbr_reads += 1;
    }

    nbr_reads
}

fn main() {
    if env::args().len() > 1 && !env::args().skip(1).any(|help_arg| help_arg == "-h") {
        rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build_global()
            .unwrap();

        // Only keep FASTQ files ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz") from command line arguments.
        let fastq_file_names: Vec<_> = env::args()
            .skip(1)
            .map(|fastq_file_name| PathBuf::from(&fastq_file_name))
            .filter_map(|fastq_file_path| {
                // Filter out files without extension.
                let fastq_ext = fastq_file_path.extension()?.to_str()?;

                let fastq_file_stem = PathBuf::from(fastq_file_path.file_stem()?);

                match fastq_ext {
                    "fastq" => Some((fastq_file_path, "fastq", fastq_file_stem)),
                    "fq" => Some((fastq_file_path, "fq", fastq_file_stem)),
                    "gz" => {
                        let fastq_file_path_minus_gz = Path::new(&fastq_file_stem);

                        // Filter out files without second extension after ".gz".
                        let fastq_ext = fastq_file_path_minus_gz.extension()?.to_str()?;

                        let fastq_file_stem = PathBuf::from(fastq_file_path_minus_gz.file_stem()?);

                        match fastq_ext {
                            "fastq" => Some((fastq_file_path, "fastq.gz", fastq_file_stem)),
                            "fq" => Some((fastq_file_path, "fq.gz", fastq_file_stem)),
                            _ => None,
                        }
                    }
                    _ => None,
                }
            })
            .collect();

        //        println!("{:?}", &fastq_file_names);

        // Count number of reads in each FASTQ file.
        fastq_file_names
            .into_par_iter()
            .for_each(|fastq_file_name| {
                println!(
                    "{}\t{}",
                    fastq_file_name.0.display(),
                    count_nbr_reads(&fastq_file_name.0)
                )
            });
    } else {
        eprintln!("Usage: {} <input.fq> [<input2.fq> ...]", file!());
        std::process::exit(1)
    }
}
