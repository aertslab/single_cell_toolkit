#!/bin/bash
#
# Copyright (C) 2020 - Gert Hulselmans
#
# Purpose:
#   Debarcode 10x scATAC FASTQ files by adding the cell barcode from R2
#   in from of the original read name for each read in R1 and R3
#   (separated by barcode_read_name_separator (default: "-")).



debarcode_10x_scatac_fastqs () {
    local fastq_R1_filename="${1}";
    local fastq_R2_filename="${2}";
    local fastq_R3_filename="${3}";
    local fastq_output_prefix="${4}";
    local barcode_read_name_separator="${5:-:}";

    local fastq_R1_output_filename="${fastq_output_prefix%.fastq.gz}_R1.fastq.gz";
    local fastq_R2_output_filename="${fastq_output_prefix%.fastq.gz}_R2.fastq.gz";

    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:   debarcode_10x_scatac_fastqs fastq_R1 fastq_R2 fastq_R3 fastq_output_prefix [barcode_read_name_separator]\n\n';
        printf 'Purpose: Debarcode 10x scATAC FASTQ files by adding the cell barcode from R2\n';
        printf '         in from of the original read name for each read in R1 and R3\n';
        printf '         (separated by barcode_read_name_separator (default: "-")).\n\n';
        return 1;
    fi

    mawk \
        -v fastq_R1_filename="${fastq_R1_filename}" \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v fastq_R3_filename="${fastq_R3_filename}" \
        -v fastq_R1_output_filename="${fastq_R1_output_filename}" \
        -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
        -v barcode_read_name_separator="${barcode_read_name_separator}" \
    '
    BEGIN {
        read_fastq_R1_cmd = "zcat " fastq_R1_filename;
        read_fastq_R2_cmd = "zcat " fastq_R2_filename;
        read_fastq_R3_cmd = "zcat " fastq_R3_filename;

        write_fastq_R1_cmd = "bgzip -@ 4 -l 6 -c > " fastq_R1_output_filename;
        write_fastq_R2_cmd = "bgzip -@ 4 -l 6 -c > " fastq_R2_output_filename;

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
                        # Get corrected barcode.
                        CBn = split(fastq_R2_line,CB," ");
                        corrected_barcode_tag = ""
                        if(match(CB[CBn],"^CB")>0) {
                            corrected_barcode_tag = CB[CBn]
                            sub("^CB:Z:","",corrected_barcode_tag)
                        }
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
                        # Create output read names, by adding cell barcode from R2 input read before the original readname.
                        read_name_fastq_R1_output = fastq_R2_line "-" corrected_barcode_tag barcode_read_name_separator read_name_R1_full;
                        read_name_fastq_R2_output = fastq_R2_line "-" corrected_barcode_tag barcode_read_name_separator read_name_R3_full;

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



debarcode_10x_scatac_fastqs "${@}";

