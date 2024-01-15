#!/bin/bash
#
# Copyright (C) 2020-2024 - Gert Hulselmans
#
# Purpose:
#   Barcode 10x scATAC FASTQ files by adding the cell barcode from R2 to each
#   read in R1 and R3, as a comment or in front of the original read name.


set -e
set -o pipefail


decompress_fastq_cat_cmd='cat';
decompress_fastq_zcat_cmd='zcat';
decompress_fastq_igzip_cmd='igzip -c -d';
decompress_fastq_zstd_cmd='zstd -c -d';


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



add_corrected_barcode_to_read_name () {
    local fastq_R1_filename="${1}";
    local fastq_R2_filename="${2}";
    local corrected_bc_filename="${3}";
    local fastq_output_prefix="${4}";
    local interleaved="${5:-false}";
    local compress_fastq_cmd="${6:-bgzip}";

    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:\n';
        printf '    add_corrected_barcode_to_read_name \\\n';
        printf '        fastq_R1 \\\n';
        printf '        fastq_R2 \\\n';
        printf '        corrected_bc_filename \\\n';
        printf '        fastq_output_prefix \\\n';
        printf '        <interleaved [true|false]> \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|gzip|stdout|-|uncompressed]> \\\n\n';
        printf 'Purpose: Barcode scATAC R1 and R2 FASTQ files by adding raw barcode sequence, raw barcode quality\n';
        printf '         and corrected barcode sequence in SAM tag format to each read name in R1 and R2 as comment.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R1:   FASTQ R1 filename.\n';
        printf '  - fastq_R2:   FASTQ R2 filename.\n';
        printf '  - corrected_bc_filename:   File with corrected barcodes.\n';
        printf '  - fastq_output_prefix: Output prefix for FASTQ output file(s).\n';
        printf '  - interleaved:\n';
        printf '      - true:   Write one output FASTQ file with reads from R1 and R2 interleaved.\n';
        printf '      - false:  Write R1 and R2 output FASTQ file with reads from R1 and R2 respectively (default).\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'  (default)\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'\n" "${compress_fastq_pigz_cmd}";
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


    if type igzip > /dev/null 2>&1 ; then
        # Decompress gzipped FASTQ files with igzip if installed (6x faster than gzip).
        local decompress_fastq_gzipped_cmd="${decompress_fastq_igzip_cmd}";
    else
        # Decompress gzipped FASTQ files with gzip.
        local decompress_fastq_gzipped_cmd="${decompress_fastq_zcat_cmd}";
    fi


    # Detect if input FASTQ files are compressed (gzip/zstd) or not.
    if [ "${fastq_R1_filename}" != "${fastq_R1_filename%.gz}" ] ; then
        local decompress_R1_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    elif [ "${fastq_R1_filename}" != "${fastq_R1_filename%.zst}" ] ; then
        local decompress_R1_fastq_cmd="${decompress_fastq_zstd_cmd}";
    else
        local decompress_R1_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi

    if [ "${fastq_R2_filename}" != "${fastq_R2_filename%.gz}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    elif [ "${fastq_R2_filename}" != "${fastq_R2_filename%.zst}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_zstd_cmd}";
    else
        local decompress_R2_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi

    # Detect if input corrected barcode file is compressed (gzip/zstd) or not.
    if [ "${corrected_bc_filename}" != "${corrected_bc_filename%.gz}" ] ; then
        local decompress_corrected_bc_file_cmd="${decompress_fastq_gzipped_cmd}";
    elif [ "${corrected_bc_filename}" != "${corrected_bc_filename%.zst}" ] ; then
        local decompress_corrected_bc_file_cmd="${decompress_fastq_zstd_cmd}";
    else
        local decompress_corrected_bc_file_cmd="${decompress_fastq_cat_cmd}";
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
            local compress_fastq_cmd="cat";;
        uncompressed)
            # Write uncompressed FASTQ files.
            fastq_R1_output_filename="${fastq_R1_output_filename%.gz}";
            fastq_R2_output_filename="${fastq_R2_output_filename%.gz}";
            local compress_fastq_cmd="cat";;
    esac


    if ! type mawk > /dev/null 2>&1 ; then
        printf 'Error: "mawk" not found or executable.\n';
        return 1;
    fi

    if ! type "${decompress_corrected_bc_file_cmd%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${decompress_corrected_bc_file_cmd%% *}";
        return 1;
    fi

    if ! type "${compress_fastq_cmd%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${compress_fastq_cmd%% *}";
        return 1;
    fi


    mawk \
        -v fastq_R1_filename="${fastq_R1_filename}" \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v corrected_bc_filename="${corrected_bc_filename}" \
        -v fastq_R1_output_filename="${fastq_R1_output_filename}" \
        -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
        -v decompress_R1_fastq_cmd="${decompress_R1_fastq_cmd}" \
        -v decompress_R2_fastq_cmd="${decompress_R2_fastq_cmd}" \
        -v decompress_corrected_bc_file_cmd="${decompress_corrected_bc_file_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
    '
    BEGIN {
        read_fastq_R1_cmd = decompress_R1_fastq_cmd " " fastq_R1_filename;
        read_fastq_R2_cmd = decompress_R2_fastq_cmd " " fastq_R2_filename;
        read_corrected_bc_file_cmd = decompress_corrected_bc_file_cmd " " corrected_bc_filename;

        write_fastq_R1_cmd = compress_fastq_cmd " > " fastq_R1_output_filename;
        write_fastq_R2_cmd = compress_fastq_cmd " > " fastq_R2_output_filename;

        read_name_corrected_bc = "";
        corrected_bc_line = "";

        fastq_line_number = 0;
        corrected_bc_line_number = 0;

        # Read FASTQ R1 file (which contains read 1).
        while ( (read_fastq_R1_cmd | getline fastq_R1_line) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            # Read FASTQ R2 file (which contains read 2).
            if ( (read_fastq_R2_cmd | getline fastq_R2_line) > 0 ) {
                if ( fastq_part == 1 ) {
                    # Find first space position (0 if no comment found) in read name from all input FASTQ files.
                    read_name_R1_space_pos = index(fastq_R1_line, " ");
                    read_name_R2_space_pos = index(fastq_R2_line, " ");

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

                    # Remove CB and UMI from end of read name, if it is there.
                    sub(/:[ACGTN+]+$/, "", read_name_R1);
                    sub(/:[ACGTN+]+$/, "", read_name_R2);

                    # Check if read names match between both FASTQ files.
                    if ( read_name_R1 != read_name_R2 ) {
                        print "Error: Read name R1 (\"" read_name_R1 "\") and read name R2 (\"" read_name_R2 "\") are not paired properly (line number: " fastq_line_number ").";
                        exit(1);
                    }

                    # Keep getting the next corrected barcode line till the read name matches with the
                    # read name in the FASTQ files.
                    # This allows to add corrected barcodes after quality trimming R1 and R2 FASTQ files,
                    # while the barcode correction was done on the original barcode FASTQ file (which
                    # can have more reads as none were removed as can happen during quality trimming).
                    while ( read_name_corrected_bc != read_name_R1 ) {
                        if ( (read_corrected_bc_file_cmd | getline corrected_bc_line) > 0 ) {
                            corrected_bc_line_number += 1;

                            # Split corrected barcode line on spaces:
                            #   - read name: strip "@" from the start.
                            #   - corrected barcode info as SAM tags separated by "\t":
                            #       CR:Z:raw_bc\tCY:Z:raw_bc_qual\tCB:Z:corrected_bc_seq

                            # Find first space position (0 if no comment found) in read name from all input FASTQ files.
                            read_name_corrected_bc_space_pos = index(corrected_bc_line, " ");

                            # Extract read name from all input FASTQ files.
                            if (read_name_corrected_bc_space_pos > 0) {
                                read_name_corrected_bc = substr(corrected_bc_line, 2, read_name_corrected_bc_space_pos - 2);
                                corrected_bc_sam_tags = substr(corrected_bc_line, read_name_corrected_bc_space_pos + 1);
                            } else {
                                read_name_corrected_bc = substr(corrected_bc_line, 2);
                                corrected_bc_sam_tags = "";
                            }

                            # Remove CB and UMI from end of read name, if it is there.
                            sub(/:[ACGTN+]+$/, "", read_name_corrected_bc);
                        }
                    }
                } else if ( fastq_part == 2 ) {
                    # Sequence lines.

                    # Store sequence info from R1 and R2 for later use.
                    sequence_R1 = fastq_R1_line;
                    sequence_R2 = fastq_R2_line;
                } else if ( fastq_part == 0 ) {
                    # Quality lines.

                    # Write the full FASTQ record to the R1 and R2 output FASTQ file with barcode info in the read name comments.
                    # When write_fastq_R1_cmd and write_fastq_R2_cmd are the same, an interleaved FASTQ fille will be written.
                    print "@" read_name_R1 " " corrected_bc_sam_tags "\n" sequence_R1 "\n+\n" fastq_R1_line | write_fastq_R1_cmd;
                    print "@" read_name_R2 " " corrected_bc_sam_tags "\n" sequence_R2 "\n+\n" fastq_R2_line | write_fastq_R2_cmd;
                }
            }
        }

        # Close open file handles.
        close(read_fastq_R1_cmd);
        close(read_fastq_R2_cmd);
        close(read_corrected_bc_file_cmd);
        close(write_fastq_R1_cmd);
        close(write_fastq_R2_cmd);
    }'

    return $?
}



add_corrected_barcode_to_read_name "${@}";

