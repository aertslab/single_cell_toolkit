#!/bin/bash
#
# Copyright (C) 2020-2021 - Gert Hulselmans
#
# Purpose:
#   Filter FASTQ by list of read names.


decompress_fastq_cat_cmd='cat';
decompress_fastq_zcat_cmd='zcat';
decompress_fastq_igzip_cmd='igzip -c -d';


# Number of threads to use to compress each FASTQ output file.
compress_fastq_threads="${compress_fastq_threads:-4}";

# Gzip compression level for bgzip, pigz and gzip.
compress_fastq_level="${compress_fastq_level:-6}";
# Gzip compression level for igzip (3 is maximum).
compress_fastq_igzip_level="3";

compress_fastq_bgzip_cmd="bgzip -@ ${compress_fastq_threads} -l ${compress_fastq_level} -c";
compress_fastq_pigz_cmd="pigz -p ${compress_fastq_threads} -${compress_fastq_level} -c";
compress_fastq_igzip_cmd="igzip -${compress_fastq_igzip_level} -c";
compress_fastq_gzip_cmd="gzip -${compress_fastq_level} -c";



filter_fastq_by_read_names () {
    local fastq_read_names_filename="${1}";
    local fastq_R1_filename="${2}";
    local fastq_R2_filename="${3}";
    local fastq_R3_filename="${4}";
    local fastq_output_prefix="${5}";
    local compress_fastq_cmd="${6:-pigz}";

    if [ ${#@} -lt 5 ] ; then
        printf '\nUsage:\n';
        printf '    filter_fastq_by_read_names \\\n';
        printf '        fastq_read_names \\\n';
        printf '        fastq_R1 \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_R3 \\\n';
        printf '        fastq_output_prefix \\\n';
        printf '        <compress_fastq_cmd> \\\n\n';
        printf 'Purpose: Filter FASTQ by list of read names.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_read_names:  File with FASTQ read names to keep.\n';
        printf '  - fastq_R1:   FASTQ R1 filename (uncompressed or gzipped).\n';
        printf '  - fastq_R2:   FASTQ R2 filename with barcodes (uncompressed or gzipped).\n';
        printf '  - fastq_R3:   FASTQ R3 filename (uncompressed or gzipped).\n';
        printf '  - fastq_output_prefix: Output prefix for FASTQ output file(s).\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'  (default)\n" "${compress_fastq_pigz_cmd}";
        printf "          - \"igzip\":  '%s'  (very fast, low compression)\n" "${compress_fastq_igzip_cmd}";
        printf "          - \"gzip\":   '%s'\n" "${compress_fastq_gzip_cmd}";
        printf '          - "stdout":  Write uncompressed output to stdout.\n';
        printf '          - "-":       Write uncompressed output to stdout.\n';
        printf '          - "uncompressed":  Write uncompressed FASTQ files.\n';
        printf '          - full custom command\n\n';
        printf '        To change number of compression threads:\n';
        printf '          - export compress_fastq_threads="%s"\n\n' "${compress_fastq_threads}";
        printf '        To change compression level:\n';
        printf '          - export compress_fastq_level="%s"\n\n' "${compress_fastq_level}";
        return 1;
    fi


    local fastq_R1_output_filename="${fastq_output_prefix%.fastq.gz}_R1.fastq.gz";
    local fastq_R2_output_filename="${fastq_output_prefix%.fastq.gz}_R2.fastq.gz";
    local fastq_R3_output_filename="${fastq_output_prefix%.fastq.gz}_R3.fastq.gz";

    if type igzip > /dev/null 2>&1 ; then
        # Decompress gzipped FASTQ files with igzip if installed (6x faster than gzip).
        local decompress_fastq_gzipped_cmd="${decompress_fastq_igzip_cmd}";
    else
        # Decompress gzipped FASTQ files with gzip.
        local decompress_fastq_gzipped_cmd="${decompress_fastq_zcat_cmd}";
    fi


    # Detect if file with FASTQ read names are gzip compressed or not.
    if [ "${fastq_read_names_filename}" != "${fastq_read_names_filename=%.gz}" ] ; then
        local decompress_fastq_read_names_cmd="${decompress_fastq_gzipped_cmd}";
    else
        local decompress_fastq_read_names_cmd="${decompress_fastq_cat_cmd}";
    fi

    # Detect if input FASTQ files are gzip compressed or not.
    if [ "${fastq_R1_filename}" != "${fastq_R1_filename%.gz}" ] ; then
        local decompress_R1_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    else
        local decompress_R1_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi

    if [ "${fastq_R2_filename}" != "${fastq_R2_filename%.gz}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    else
        local decompress_R2_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi

    if [ "${fastq_R3_filename}" != "${fastq_R3_filename%.gz}" ] ; then
        local decompress_R3_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    else
        local decompress_R3_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi


    case "${compress_fastq_cmd}" in
        bgzip)
            local compress_fastq_cmd="${compress_fastq_bgzip_cmd}";;
        pigz)
            local compress_fastq_cmd="${compress_fastq_pigz_cmd}";;
        igzip)
            local compress_fastq_cmd="${compress_fastq_igzip_cmd}";;
        gzip)
            local compress_fastq_cmd="${compress_fastq_gzip_cmd}";;
        stdout|-)
            # Write interleaved FASTQ files when writing to stdout.
            fastq_R1_output_filename='/dev/stdout';
            fastq_R2_output_filename='/dev/stdout';
            fastq_R3_output_filename='/dev/stdout';
            local compress_fastq_cmd="cat";;
        uncompressed)
            # Write uncompressed FASTQ files.
            fastq_R1_output_filename="${fastq_R1_output_filename%.gz}";
            fastq_R2_output_filename="${fastq_R2_output_filename%.gz}";
            fastq_R3_output_filename="${fastq_R3_output_filename%.gz}";
            local compress_fastq_cmd="cat";;
    esac


    mawk \
        -v fastq_read_names_filename="${fastq_read_names_filename}" \
        -v fastq_R1_filename="${fastq_R1_filename}" \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v fastq_R3_filename="${fastq_R3_filename}" \
        -v fastq_R1_output_filename="${fastq_R1_output_filename}" \
        -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
        -v fastq_R3_output_filename="${fastq_R3_output_filename}" \
        -v decompress_fastq_read_names_cmd="${decompress_fastq_read_names_cmd}" \
        -v decompress_R1_fastq_cmd="${decompress_R1_fastq_cmd}" \
        -v decompress_R2_fastq_cmd="${decompress_R2_fastq_cmd}" \
        -v decompress_R3_fastq_cmd="${decompress_R3_fastq_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
    '
    BEGIN {
        read_fastq_read_names_cmd = decompress_fastq_read_names_cmd " " fastq_read_names_filename;

        read_fastq_R1_cmd = decompress_R1_fastq_cmd " " fastq_R1_filename;
        read_fastq_R2_cmd = decompress_R2_fastq_cmd " " fastq_R2_filename;
        read_fastq_R3_cmd = decompress_R3_fastq_cmd " " fastq_R3_filename;

        write_fastq_R1_cmd = compress_fastq_cmd " > " fastq_R1_output_filename;
        write_fastq_R2_cmd = compress_fastq_cmd " > " fastq_R2_output_filename;
        write_fastq_R3_cmd = compress_fastq_cmd " > " fastq_R3_output_filename;

        fastq_line_number = 0;

        while ( (read_fastq_read_names_cmd | getline fastq_read_names_line) > 0 ) {
            if ( fastq_read_names_line != "" && fastq_read_names_line !~ /^#/ ) {
                fastq_read_names_to_keep[fastq_read_names_line] = 1;
            }
        }

        # Close open file handle.
        close(read_fastq_read_names_cmd);

        # Read FASTQ R1 file.
        while ( (read_fastq_R1_cmd | getline fastq_R1_line) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            # Read FASTQ R2 file.
            if ( (read_fastq_R2_cmd | getline fastq_R2_line) > 0 ) {
                # Read FASTQ R3 file.
                if ( (read_fastq_R3_cmd | getline fastq_R3_line) > 0 ) {
                    if ( fastq_part == 1 ) {
                        # Read name lines.

                        # Find first space position (0 if no comment found) in read name from all input FASTQ files.
                        read_name_R1_space_pos = index(fastq_R1_line, " ");
                        read_name_R2_space_pos = index(fastq_R2_line, " ");
                        read_name_R3_space_pos = index(fastq_R3_line, " ");

                        # Extract read name from all input FASTQ files.
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

                        if (read_name_R3_space_pos > 0) {
                            read_name_R3 = substr(fastq_R3_line, 2, read_name_R3_space_pos - 2);
                        } else {
                            read_name_R3 = substr(fastq_R3_line, 2);
                        }

                        # Check if read names match between all 3 FASTQ files.
                        if ( read_name_R1 != read_name_R2 || read_name_R1 != read_name_R3 ) {
                            print "Error: Read name R1 (\"" read_name_R1 "\"), read name R2 (\"" read_name_R2 "\") and R3 (\"" read_name_R3 "\") are not paired properly (line number: " fastq_line_number ").";
                            exit(1);
                        }

                        # Check if this read needs to be kept or not.
                        if ( read_name_R1 in fastq_read_names_to_keep ) {
                            full_read_name_R1_line = fastq_R1_line;
                            full_read_name_R2_line = fastq_R2_line;
                            full_read_name_R3_line = fastq_R3_line;

                            keep_read = 1;
                        } else {
                            keep_read = 0;
                        }
                    } else if ( keep_read == 0) {
                        # Skip unwanted read.
                        continue;
                    } else if ( fastq_part == 2 ) {
                        # Sequence lines.

                        # Store sequence info from R1, R2 and R3 for later use.
                        sequence_R1 = fastq_R1_line;
                        sequence_R2 = fastq_R2_line;
                        sequence_R3 = fastq_R3_line;
                    } else if ( fastq_part == 0 ) {
                        # Quality lines.

                        # Write the full FASTQ record to the R1, R2 and R3 output FASTQ file.
                        print full_read_name_R1_line "\n" sequence_R1 "\n+\n" fastq_R1_line | write_fastq_R1_cmd;
                        print full_read_name_R2_line "\n" sequence_R2 "\n+\n" fastq_R2_line | write_fastq_R2_cmd;
                        print full_read_name_R3_line "\n" sequence_R3 "\n+\n" fastq_R3_line | write_fastq_R3_cmd;
                    }
                }
            }
        }

        # Close open file handles.
        close(read_fastq_R1_cmd);
        close(read_fastq_R2_cmd);
        close(read_fastq_R3_cmd);
        close(write_fastq_R1_cmd);
        close(write_fastq_R2_cmd);
        close(write_fastq_R3_cmd);
    }'

    return $?
}



filter_fastq_by_read_names "${@}";

