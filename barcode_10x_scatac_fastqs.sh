#!/bin/bash
#
# Copyright (C) 2020-2021 - Gert Hulselmans
#
# Purpose:
#   Barcode 10x scATAC FASTQ files by adding the cell barcode from R2 to each
#   read in R1 and R3, as a comment or in front of the original read name.


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
    local interleaved="${5:-true}";
    local add_barcode_in_comment="${6:-true}";
    local compress_fastq_cmd="${8:-bgzip}";

    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:\n';
        printf '    barcode_10x_scatac_fastqs \\\n';
        printf '        fastq_R1 \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_R3 \\\n';
        printf '        fastq_output_prefix \\\n';
        printf '        <interleaved [true|false]> \\\n';
        printf '        <add_barcode_in_comment [true|false]> \\\n';
        printf '        <barcode_tags_or_separator> \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|gzip|stdout|-]> \\\n\n';
        printf 'Purpose: Barcode 10x scATAC FASTQ files by adding the cell barcode from R2 to each\n';
        printf '         read in R1 and R3, as a comment or in front of the original read name.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R1:   FASTQ R1 fileanme (uncompressed or gzipped).\n';
        printf '  - fastq_R2:   FASTQ R2 fileanme with barcodes (uncompressed or gzipped).\n';
        printf '  - fastq_R3:   FASTQ R3 fileanme (uncompressed or gzipped).\n';
        printf '  - fastq_output_prefix: Output prefix for FASTQ output file(s).\n';
        printf '  - interleaved:\n';
        printf '      - true:   Write one output FASTQ file with reads from R1 and R3 interleaved (default).\n';
        printf '      - false:  Write R1 and R2 output FASTQ file with reads from R1 and R3 respectively.\n';
        printf '  - add_barcode_in_comment:\n';
        printf '      - true:   Add barcode and barcode quality from R2 as read comment in SAM spec format\n';
        printf '                (for usage wtih: "bwa -C") (default).\n';
        printf '      - false:  Add barcode from R2 at the beginning of the read (for usage with bap/sinto).\n';
        printf '  - barcode_tags_or_separator:\n';
        printf '      - If add_barcode_in_comment = "true":\n';
        printf '          Specify barcode tag (2 characters) and barcode quality tag (2 characters)\n';
        printf '          (default: "CR_CY").\n';
        printf '      - If add_barcode_in_comment = "false":\n';
        printf '          Specify string which will be added between barcode name and the original read name\n';
        printf '          (default: "-").\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'  (default)\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'\n" "${compress_fastq_pigz_cmd}";
        printf "          - \"gzip\":   '%s'\n" "${compress_fastq_gzip_cmd}";
        printf '          - "stdout":  Write uncompressed output to stdout.\n';
        printf '          - "-":       Write uncompressed output to stdout.\n';
        printf '          - full custom command\n\n';
        printf '        To change number of compression threads:\n';
        printf '          - export compress_fastq_threads="%s"\n\n' "${compress_fastq_threads}";
        printf '        To change comprssion level:\n';
        printf '          - export compress_fastq_level="%s"\n\n' "${compress_fastq_level}";
        return 1;
    fi


    if [ "${add_barcode_in_comment}" = "true" ] ; then
        local barcode_tags="${7:-CR_CY}";

        if [ "${#barcode_tags}" -ne 5 ] ; then
            printf 'Error: barcode tags field "%s" is wrongly formatted.\n' "${barcode_tags}" >&2;
            return 1;
        fi

        # Extract barcode tag and barcode quality tag.
        local barcode_tag="${barcode_tags:0:2}";
        local barcode_qual_tag="${barcode_tags:3:2}";

        local barcode_read_name_separator='';

        # Redefine for mawk.
        add_barcode_in_comment=1;
    elif [ "${add_barcode_in_comment}" = "false" ] ; then
        local barcode_read_name_separator="${7:--}";

        local barcode_tag='';
        local barcode_qual_tag='';

        # Redefine for mawk.
        add_barcode_in_comment=0;
    else
        printf 'Error: "add_barcode_in_comment" has an unsupported value: "%s".\n' "${add_barcode_in_comment}" >&2;
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
    esac


    mawk \
        -v fastq_R1_filename="${fastq_R1_filename}" \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v fastq_R3_filename="${fastq_R3_filename}" \
        -v fastq_R1_output_filename="${fastq_R1_output_filename}" \
        -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
        -v barcode_tag="${barcode_tag}" \
        -v barcode_qual_tag="${barcode_qual_tag}" \
        -v barcode_read_name_separator="${barcode_read_name_separator}" \
        -v decompress_R1_fastq_cmd="${decompress_R1_fastq_cmd}" \
        -v decompress_R2_fastq_cmd="${decompress_R2_fastq_cmd}" \
        -v decompress_R3_fastq_cmd="${decompress_R3_fastq_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
        -v add_barcode_in_comment="${add_barcode_in_comment} " \
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
                        if ( read_name_R1 == read_name_R2 == read_name_R3 ) {
                            print "Error: Read name R1 (\"" read_name_R1 "\"), read name R2 (\"" read_name_R2 "\") and R3 (\"" read_name_R3 "\") are not paired properly (line number: " fastq_line_number ")." 
                            exit(1);
                        }

                        if (add_barcode_in_comment == 1) {
                            # Store read name comments for R2 input FASTQ file so we can extract existing comments, like a corrected barcode.
                            read_name_R2_comment = substr(fastq_R2_line, length(read_name_R2) + 2);
                        } else {
                            # Store full read names for R1 and R3 input FASTQ files.
                            read_name_R1_full = substr(fastq_R1_line, 2);
                            read_name_R3_full = substr(fastq_R3_line, 2);
                        }
                    } else if ( fastq_part == 2 ) {
                        # Sequence lines.

                        # Store sequence info from R1, R2 and R3 for later use.
                        sequence_R1 = fastq_R1_line;
                        sequence_R2 = fastq_R2_line;
                        sequence_R3 = fastq_R3_line;
                    } else if ( fastq_part == 0 ) {
                        # Quality lines.

                        if (add_barcode_in_comment == 1) {
                            # Store barcode sequence and barcode quality in SAM spec format.
                            read_name_comment = sprintf("%s:Z:%s\t%s:Z:%s", barcode_tag, sequence_R2, barcode_qual_tag, fastq_R2_line);

                            # Split R2 read name comment on spaces and tabs so we can check if there are existing comments in SAM spec format.
                            nbr_splits = split(read_name_R2_comment, read_name_R2_comment_array, / |\t/);

                            for ( i = 1; i <= nbr_splits; i++ ) {
                                if ( $i ~ /^[A-Z]{2}:[ABfHiZ]:/ ) {
                                    # Add SAM spec comment from original barcode file.
                                    read_name_comment = read_name_comment "\t" $i;
                                }
                            }

                            # Write the full FASTQ record to the R1 and R2 output FASTQ file with barcode info in the read name comments.
                            # When write_fastq_R1_cmd and write_fastq_R2_cmd are the same, an interleaved FASTQ fille will be written.
                            print "@" read_name_R1 " " read_name_comment "\n" sequence_R1 "\n+\n" fastq_R1_line | write_fastq_R1_cmd;
                            print "@" read_name_R3 " " read_name_comment "\n" sequence_R3 "\n+\n" fastq_R3_line | write_fastq_R2_cmd;
                        } else {
                            # Write the full FASTQ record to the R1 and R2 output FASTQ file with barcode info in front of the read name.
                            # When write_fastq_R1_cmd and write_fastq_R2_cmd are the same, an interleaved FASTQ fille will be written.
                            print "@" sequence_R2 barcode_read_name_separator read_name_R1 "\n" sequence_R1 "\n+\n" fastq_R1_line | write_fastq_R1_cmd;
                            print "@" sequence_R2 barcode_read_name_separator read_name_R3 "\n" sequence_R3 "\n+\n" fastq_R3_line | write_fastq_R2_cmd;
                        }
                    }
                }
            }
        }
    }'

    return $?
}



barcode_10x_scatac_fastqs "${@}";

