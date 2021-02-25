#!/bin/bash
#
# Copyright (C) 2020-2021 - Gert Hulselmans
#
# Purpose:
#   Barcode 10x scATAC FASTQ files by adding the cell barcode from R2
#   in front of the original read name for each read in R1 and R3
#   (separated by barcode_read_name_separator (default: "-")).



decompress_fastq_cat_cmd='cat';
decompress_fastq_zcat_cmd='zcat';


# Number of threads to use to compress each FASTQ output file.
compress_fastq_threads="${compress_fastq_threads:-4}";

# Gzip compression level.
compress_fastq_level="${compress_fastq_level:-6}";

compress_fastq_bgzip_cmd="bgzip -@ ${compress_fastq_threads} -l ${compress_fastq_level} -c";
compress_fastq_pigz_cmd="pigz -p ${compress_fastq_threads} -${compress_fastq_level} -c";
compress_fastq_gzip_cmd="gzip -${compress_fastq_level} -c";



barcode_10x_scatac_fastqs () {
    local fastq_R1_filename="${1}";
    local fastq_R2_filename="${2}";
    local fastq_R3_filename="${3}";
    local fastq_output_prefix="${4}";
    local barcode_read_name_separator="${5:-:}";
    local compress_fastq_cmd="${6:-bgzip}";

    local fastq_R1_output_filename="${fastq_output_prefix%.fastq.gz}_R1.fastq.gz";
    local fastq_R2_output_filename="${fastq_output_prefix%.fastq.gz}_R2.fastq.gz";

    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:   barcode_10x_scatac_fastqs fastq_R1 fastq_R2 fastq_R3 fastq_output_prefix [barcode_read_name_separator] [bgzip|pigz|gzip]\n\n';
        printf 'Purpose: Barcode 10x scATAC FASTQ files by adding the cell barcode from R2\n';
        printf '         in front of the original read name for each read in R1 and R3\n';
        printf '         (separated by barcode_read_name_separator (default: "-")).\n\n';
        printf '         Compression program to use for output FASTQ files:\n';
        printf "           - bgzip:  '%s'  (default)\n" "${compress_fastq_bgzip_cmd}";
        printf "           - pigz:   '%s'\n" "${compress_fastq_pigz_cmd}";
        printf "           - gzip:   '%s'\n" "${compress_fastq_gzip_cmd}";
        printf '           - full custom command\n\n';
        printf '         To change number of compression threads:\n';
        printf '           - export compress_fastq_threads="%s"\n\n' "${compress_fastq_threads}";
        printf '         To change comprssion level:\n';
        printf '           - export compress_fastq_level="%s"\n\n' "${compress_fastq_level}";
        return 1;
    fi

    # Detect if input FASTQ files are gzip compressed or not.
    if [ "${fastq_R1_filename}" != "${fastq_R1_filename%.gz}" ] ; then
        local decompress_R1_fastq_cmd="${decompress_fastq_zcat_cmd}";
    else
        local decompress_R1_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi

    if [ "${fastq_R2_filename}" != "${fastq_R2_filename%.gz}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_zcat_cmd}";
    else
        local decompress_R2_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi

    if [ "${fastq_R3_filename}" != "${fastq_R3_filename%.gz}" ] ; then
        local decompress_R3_fastq_cmd="${decompress_fastq_zcat_cmd}";
    else
        local decompress_R3_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi


    case "${compress_fastq_cmd}" in
        bgzip)  local compress_fastq_cmd="${compress_fastq_bgzip_cmd}";;
        pigz)   local compress_fastq_cmd="${compress_fastq_pigz_cmd}";;
        gzip)   local compress_fastq_cmd="${compress_fastq_gzip_cmd}";;
    esac


    mawk \
        -v fastq_R1_filename="${fastq_R1_filename}" \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v fastq_R3_filename="${fastq_R3_filename}" \
        -v fastq_R1_output_filename="${fastq_R1_output_filename}" \
        -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
        -v barcode_read_name_separator="${barcode_read_name_separator}" \
        -v decompress_R1_fastq_cmd="${decompress_R1_fastq_cmd}" \
        -v decompress_R2_fastq_cmd="${decompress_R2_fastq_cmd}" \
        -v decompress_R3_fastq_cmd="${decompress_R3_fastq_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
    '
    BEGIN {
        read_fastq_R1_cmd = decompress_R1_fastq_cmd " " fastq_R1_filename;
        read_fastq_R2_cmd = decompress_R2_fastq_cmd " " fastq_R2_filename;
        read_fastq_R3_cmd = decompress_R3_fastq_cmd " " fastq_R3_filename;

        write_fastq_R1_cmd = compress_fastq_cmd " > " fastq_R1_output_filename;
        write_fastq_R2_cmd = compress_fastq_cmd " > " fastq_R2_output_filename;

        fastq_line_number = 0;

        # Read FASTQ R2 file (which contains the cell barcodes).
        while ( (read_fastq_R2_cmd | getline fastq_R2_line) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            # Read FASTQ R1 file (which contains read 1).
            if ( (read_fastq_R1_cmd | getline fastq_R1_line) > 0 ) {
                # Read FASTQ R3 file (which contains read 2).
                if ( (read_fastq_R3_cmd | getline fastq_R3_line) > 0 ) {
                    if ( fastq_part == 1 ) {
                        # Extract read name from all input FASTQ files.
                        read_name_R1 = substr(fastq_R1_line, 2, index(fastq_R1_line, " ") - 2);
                        read_name_R2 = substr(fastq_R2_line, 2, index(fastq_R2_line, " ") - 2);
                        read_name_R3 = substr(fastq_R3_line, 2, index(fastq_R3_line, " ") - 2);

                        # Check if read names match between all 3 FASTQ files.
                        if ( read_name_R1 == read_name_R2 == read_name_R3 ) {
                            print "Error: Read name R1 (\"" read_name_R1 "\"), read name R2 (\"" read_name_R2 "\") and R3 (\"" read_name_R3 "\") are not paired properly (line number: " fastq_line_number ")." 
                            exit(1);
                        }

                        # Store full read names for R1 and R3 input FASTQ files.
                        read_name_R1_full = substr(fastq_R1_line, 2);
                        read_name_R3_full = substr(fastq_R3_line, 2);
                    } else if ( fastq_part == 2 ) {
                        # Create output read names, by adding cell barcode from R2 input read before the original read name.
                        read_name_fastq_R1_output = fastq_R2_line barcode_read_name_separator read_name_R1_full;
                        read_name_fastq_R2_output = fastq_R2_line barcode_read_name_separator read_name_R3_full;

                        # Store sequence info from R1 and R3 for later use.
                        seq_fastq_R1_output = fastq_R1_line;
                        seq_fastq_R2_output = fastq_R3_line;
                    } else if ( fastq_part == 0 ) {
                        # Write the full FASTQ record to the R1 and R2 output FASTQ file.
                        print "@" read_name_fastq_R1_output "\n" seq_fastq_R1_output "\n+\n" fastq_R1_line | write_fastq_R1_cmd;
                        print "@" read_name_fastq_R2_output "\n" seq_fastq_R2_output "\n+\n" fastq_R3_line | write_fastq_R2_cmd;
                    }
                }
            }
        }
    }'

    return $?
}



barcode_10x_scatac_fastqs "${@}";

