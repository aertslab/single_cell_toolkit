#!/bin/bash
#
# Copyright (C) 2021 - Gert Hulselmans
#
# Purpose:
#   Extract barcodes and ATAC part of read from BioRad v2.1 FASTQ files,
#   correct barcodes and write new FASTQ files with ATAC part as sequence
#   and barcode info in FASTQ comment field.


decompress_fastq_cat_cmd='cat';
decompress_fastq_zcat_cmd='zcat';


# Number of threads to use to compress each FASTQ output file.
compress_fastq_threads="${compress_fastq_threads:-4}";

# Gzip compression level.
compress_fastq_level="${compress_fastq_level:-6}";

compress_fastq_bgzip_cmd="bgzip -@ ${compress_fastq_threads} -l ${compress_fastq_level} -c";
compress_fastq_pigz_cmd="pigz -p ${compress_fastq_threads} -${compress_fastq_level} -c";
compress_fastq_gzip_cmd="gzip -${compress_fastq_level} -c";



extract_and_correct_biorad_barcode_in_fastq () {
    local fastq_R1_filename="${1}";
    local fastq_R2_filename="${2}";
    local fastq_output_prefix="${3}";
    local interleaved="${4:-true}";
    local compress_fastq_cmd="${5:-bgzip}";

    if [ ${#@} -lt 3 ] ; then
        printf '\nUsage:\n';
        printf '    extract_and_correct_biorad_barcode_in_fastq \\\n';
        printf '        fastq_R1 \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_output_prefix \\\n';
        printf '        <interleaved [true|false]> \\\n';
        printf '        <compress_fastq_cmd <[bgzip|pigz|gzip|stdout|-|uncompressed]>\n\n';
        printf 'Purpose: Extract barcodes and ATAC part of read from BioRad v2.1 FASTQ files,\n';
        printf '         correct barcodes and write new FASTQ files with ATAC part as sequence\n';
        printf '         and barcode info in FASTQ comment field.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R1:   BioRad v2.1 FASTQ R1 filename (uncompressed or gzipped).\n';
        printf '  - fastq_R2:   BioRad v2.1 FASTQ R2 filename (uncompressed or gzipped).\n';
        printf '  - fastq_output_prefix: Output prefix for FASTQ output file(s).\n';
        printf '  - interleaved:\n';
        printf '      - true:   Write one output FASTQ file with reads from R1 and R3 interleaved (default).\n';
        printf '      - false:  Write R1 and R2 output FASTQ file with reads from R1 and R3 respectively.\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'  (default)\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'\n" "${compress_fastq_pigz_cmd}";
        printf "          - \"gzip\":   '%s'\n" "${compress_fastq_gzip_cmd}";
        printf '          - "stdout":  Write uncompressed output to stdout.\n';
        printf '          - "-":       Write uncompressed output to stdout.\n';
        printf '          - "uncompressed":  Write uncompressed FASTQ files.\n';
        printf '          - full custom command\n\n';
        printf '        To change number of compression threads:\n';
        printf '          - export compress_fastq_threads="%s"\n\n' "${compress_fastq_threads}";
        printf '        To change compression level:\n';
        printf '          - export compress_fastq_level="%s"\n\n' "${compress_fastq_level}";
        return $1;
    fi


    if [ ! -e "${fastq_R1_filename}" ] ; then
        printf 'Error: BioRad FASTQ R1 file "%s" could not be found.\n' "${fastq_R1_filename}" >&2;
        return 1;
    fi

    if [ ! -e "${fastq_R2_filename}" ] ; then
        printf 'Error: BioRad FASTQ R2 file "%s" could not be found.\n' "${fastq_R2_filename}" >&2;
        return 1;
    fi


    if [ "${interleaved}" = "true" ] ; then
        # Write to the same file.
        local fastq_R1_output_filename="${fastq_output_prefix%.fastq.gz}_interleaved.fastq.gz";
        local fastq_R2_output_filename="${fastq_output_prefix%.fastq.gz}_interleaved.fastq.gz";
    elif [ "${interleaved}" = "false" ] ; then
        local fastq_R1_output_filename="${fastq_output_prefix%.fastq.gz}_R1.fastq.gz";
        local fastq_R2_output_filename="${fastq_output_prefix%.fastq.gz}_R2.fastq.gz";
    else
        printf 'Error: "interleaved" has an unsupported value: "%s".\n' "${interleaved}" >&2;
        return 1;
    fi

    local corrected_bc_stats_tsv_filename="${fastq_output_prefix%.fastq.gz}.corrected_bc_stats.tsv"


    # Detect if input FASTQ files are gzip compressed or not.
    if [ "${fastq_R2_filename}" != "${fastq_R2_filename%.gz}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_zcat_cmd}";
    else
        local decompress_R2_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi


    case "${compress_fastq_cmd}" in
        bgzip)
            local compress_fastq_cmd="${compress_fastq_bgzip_cmd}";;
        pigz)
            local compress_fastq_cmd="${compress_fastq_pigz_cmd}";;
        gzip)
            local compress_fastq_cmd="${compress_fastq_gzip_cmd}";;
        stdout|-)
            # Write interleaved FASTQ files when writing to stdout.
            local fastq_R1_output_filename='/dev/stdout';
            local fastq_R2_output_filename='/dev/stdout';
            local compress_fastq_cmd="cat";;
        uncompressed)
            # Write uncompressed FASTQ files.
            local fastq_R1_output_filename="${fastq_R1_output_filename%.gz}";
            local fastq_R2_output_filename="${fastq_R2_output_filename%.gz}";
            local compress_fastq_cmd="cat";;
    esac

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";


    # First, extract barcodes and ATAC part of read from BioRad v2.1 R1 FASTQ file,
    # then correct barcodes and write new FASTQ files with ATAC part as sequence
    # and barcode info in FASTQ comment field.
    # Pipe the new R1 FASTQ into mawk, fix read name comments in R2 reads and
    # write both fixed FASTQ files.
    "${script_dir}/run_seq_program.sh" \
        "${script_dir}/extract_and_correct_biorad_barcode_in_fastq.seq" \
            "${fastq_R1_filename}" \
            '/dev/stdout' \
            "${corrected_bc_stats_tsv_filename}" \
      | mawk -F ' ' \
            -v fastq_R2_filename="${fastq_R2_filename}" \
            -v fastq_R1_output_filename="${fastq_R1_output_filename}" \
            -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
            -v decompress_R2_fastq_cmd="${decompress_R2_fastq_cmd}" \
            -v compress_fastq_cmd="${compress_fastq_cmd}" \
            '
            BEGIN {
                read_fastq_R2_cmd = decompress_R2_fastq_cmd " " fastq_R2_filename;

                write_fastq_R1_cmd = compress_fastq_cmd " > " fastq_R1_output_filename;
                write_fastq_R2_cmd = compress_fastq_cmd " > " fastq_R2_output_filename;

                fastq_line_number = 0;

                # Read FASTQ R1 file (which contains the barcode info in the comments
                # and the first ATAC read in the sequence).
                while ( getline fastq_R1_line > 0 ) {
                    fastq_line_number += 1;
                    fastq_part = fastq_line_number % 4;

                    # Read FASTQ R2 file (which contains the second ATAC read).
                    if ( (read_fastq_R2_cmd | getline fastq_R2_line) > 0 ) {
                        if ( fastq_part == 1 ) {
                            # Read name lines.

                            # Find first space position (0 if no comment found) in read name from both input FASTQ files.
                            read_name_R1_space_pos = index(fastq_R1_line, " ");
                            read_name_R2_space_pos = index(fastq_R2_line, " ");

                            # Extract read name from both input FASTQ files.
                            if (read_name_R1_space_pos > 0) {
                                read_name_R1 = substr(fastq_R1_line, 2, read_name_R1_space_pos - 2);
                            } else {
                                read_name_R1 = substr(fastq_R1_line, 2);
                            }

                            if (read_name_R2_space_pos > 0) {
                                read_name_R2 = substr(fastq_R2_line, 2, read_name_R2_space_pos - 2);
                            } else {
                                read_name_R2 = substr(fastq_R2_line, 2);
                            }

                            # Check if read names match between R1 and R2 FASTQ files.
                            if ( read_name_R1 != read_name_R2 ) {
                                print "Error: Read name R1 (\"" read_name_R1 "\") and read name R2 (\"" read_name_R2 "\") are not paired properly (line number: " fastq_line_number ").";
                                exit(1);
                            }

                            # Store R1 readname and comment, so we can write it to R1 and R2 later.
                            read_name_and_comment = fastq_R1_line;
                        } else if ( fastq_part == 2 ) {
                            # Sequence lines.

                            # Store sequence info from R1 and R2 for later use.
                            sequence_R1 = fastq_R1_line;
                            sequence_R2 = fastq_R2_line;
                        } else if ( fastq_part == 0 ) {
                            # Quality lines.

                            # Write the full FASTQ record to the R1 and R2 output FASTQ file with barcode info in the read name comments.
                            # When write_fastq_R1_cmd and write_fastq_R2_cmd are the same, an interleaved FASTQ fille will be written.
                            print read_name_and_comment "\n" sequence_R1 "\n+\n" fastq_R1_line | write_fastq_R1_cmd;
                            print read_name_and_comment "\n" sequence_R2 "\n+\n" fastq_R2_line | write_fastq_R2_cmd;
                        }
                    }
                }

                # Close open file handles.
                close(read_fastq_R2_cmd);
                close(write_fastq_R1_cmd);
                close(write_fastq_R2_cmd);
            }';

    return $?;
}



extract_and_correct_biorad_barcode_in_fastq "${@}";
