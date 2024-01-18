#!/bin/bash
#
# Copyright (C) 2020-2024 - Gert Hulselmans
#
# Purpose:
#   Split ScaleBio ATAC FASTQ R1 and R2 files and corrected tagmentation
#   barcode file per tagmentation barcode.


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



split_scalebio_atac_fastqs () {
    local fastq_R1_filename="${1}";
    local fastq_R2_filename="${2}";
    local corrected_bc_filename="${3}";
    local fastq_output_dir="${4}";
    local compress_fastq_cmd="${5:-bgzip}";


    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:\n';
        printf '    split_scalebio_atac_fastqs \\\n';
        printf '        fastq_R1 \\\n';
        printf '        fastq_R2 \\\n';
        printf '        corrected_bc_filename \\\n';
        printf '        fastq_output_dir \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|gzip|uncompressed]> \\\n\n';
        printf 'Purpose: Split ScaleBio ATAC FASTQ R1 and R2 files and corrected tagmentation\n'
        printf '         barcode file per tagmentation barcode.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R1:   FASTQ R1 filename.\n';
        printf '  - fastq_R2:   FASTQ R2 filename.\n';
        printf '  - corrected_bc_filename:   File with corrected barcodes.\n';
        printf '  - fastq_output_dir: Output dir for FASTQ output file(s).\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'  (default)\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'\n" "${compress_fastq_pigz_cmd}";
        printf "          - \"igzip\":  '%s'  (very fast, low compression)\n" "${compress_fastq_igzip_cmd}";
        printf "          - \"gzip\":   '%s'\n" "${compress_fastq_gzip_cmd}";
        printf '          - "uncompressed":  Write uncompressed FASTQ files.\n';
        printf '          - full custom command\n\n';
        printf '        To change number of compression threads:\n';
        printf '          - export compress_fastq_threads="%s"\n\n' "${compress_fastq_threads}";
        printf '        To change compression level:\n';
        printf '          - export compress_fastq_level="%s"\n\n' "${compress_fastq_level}";
        return 1;
    fi

    local fastq_R1_basename="${fastq_R1_filename##*/}";

    # Get current extglob setting.
    local extglob_setting=$(shopt -p extglob);

    # Enable extglob setting.
    shopt -s extglob;

    # Try to extract the suffix after the sample name from the FASTQ file name,
    # so the same suffix can be used for the output FASTQ files.
    local fastq_output_basename_prefix="${fastq_R1_basename%_S+([0-9])_R1?(_00[1-9]).fastq?(.gz)}"
    local fastq_output_basename_prefix="${fastq_output_basename_prefix%_R1?(_00[1-9]).fastq?(.gz)}"

    # Restore extglob setting.
    eval ${extglob_setting};

    if [ "${#fastq_output_basename_prefix}" -eq "${#fastq_R1_basename}" ] ; then
        # FASTQ sample name could not be extracted. Set FASTQ R1 output prefix manually.
        local fastq_output_R1_suffix="_R1.fastq.gz";
    else
        # FASTQ sample name and suffix could be extracted.
        local fastq_output_R1_suffix="${fastq_R1_basename:${#fastq_output_basename_prefix}}";
    fi

    # Create output suffixes for FASTQ R2 and R3, based on FASTQ R1 suffix.
    local fastq_output_R2_suffix="${fastq_output_R1_suffix/_R1/_R2}";
    local fastq_output_R3_suffix="${fastq_output_R1_suffix/_R1/_R3}";


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
        uncompressed)
            # Write uncompressed FASTQ files.
            fastq_output_R1_suffix="${fastq_output_R1_suffix%.gz}";
            fastq_output_R2_suffix="${fastq_output_R2_suffix%.gz}";
            fastq_output_R3_suffix="${fastq_output_R3_suffix%.gz}";
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
        -v fastq_output_prefix="${fastq_output_dir%/}/${fastq_output_basename_prefix}" \
        -v fastq_output_R1_suffix="${fastq_output_R1_suffix}" \
        -v fastq_output_R2_suffix="${fastq_output_R2_suffix}" \
        -v fastq_output_R3_suffix="${fastq_output_R3_suffix}" \
        -v decompress_R1_fastq_cmd="${decompress_R1_fastq_cmd}" \
        -v decompress_R2_fastq_cmd="${decompress_R2_fastq_cmd}" \
        -v decompress_corrected_bc_file_cmd="${decompress_corrected_bc_file_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
    '
    BEGIN {
        read_fastq_R1_cmd = decompress_R1_fastq_cmd " " fastq_R1_filename;
        read_fastq_R2_cmd = decompress_R2_fastq_cmd " " fastq_R2_filename;
        read_corrected_bc_file_cmd = decompress_corrected_bc_file_cmd " " corrected_bc_filename;

        partial_write_fastq_cmd = compress_fastq_cmd " > " fastq_output_prefix;

        read_name_corrected_bc = "";
        corrected_bc_line = "";

        fastq_line_number = 0;
        corrected_bc_line_number = 0;

        read_len_R1_prev = -1;
        read_len_R2_prev = -1;
        read_len_corrected_bc_prev = -1;

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

                    read_len_R1 = length(read_name_R1);

                    if ( read_len_R1 != read_len_R1_prev ) {
                        # Try to remove CB and UMI from end of read name, if current
                        # read name length is different from the previous one.
                        sub(/:[ACGTN+]+$/, "", read_name_R1);
                        read_len_R1 = length(read_name_R1);
                        read_len_R1_prev = read_len_R1;
                    }

                    read_len_R2 = length(read_name_R2);

                    if ( read_len_R2 != read_len_R2_prev ) {
                        # Try to remove CB and UMI from end of read name, if current
                        # read name length is different from the previous one.
                        sub(/:[ACGTN+]+$/, "", read_name_R2);
                        read_len_R2 = length(read_name_R2);
                        read_len_R2_prev = read_len_R2;
                    }

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

                            read_len_corrected_bc = length(read_name_corrected_bc);

                            if ( read_len_corrected_bc != read_len_corrected_bc_prev ) {
                                # Try to remove CB and UMI from end of read name, if current
                                # read name length is different from the previous one.
                                sub(/:[ACGTN+]+$/, "", read_name_corrected_bc);
                                read_len_corrected_bc = length(read_name_corrected_bc);
                                read_len_corrected_bc_prev = read_len_corrected_bc;
                            }

                            # Split corrected barcode info on "\t".
                            n_corrected_bc_sam_tags = split(corrected_bc_sam_tags, corrected_bc_sam_tags_array, "\t");

                            if (n_corrected_bc_sam_tags == 3) {
                                cr_sam_tag = corrected_bc_sam_tags_array[1];
                                cy_sam_tag = corrected_bc_sam_tags_array[2];
                                tb_sam_tag = corrected_bc_sam_tags_array[3];

                                if ( cr_sam_tag !~ /^cr:Z:/ ) {
                                    print "Error: Expected \"cr:Z:\" instead of \"" substr(cr_sam_tag, 1, 5) "\" for uncorrected barcode sequence with read name (\"" read_name_corrected_bc "\") (line number: " corrected_bc_line_number ").";
                                    exit(1);
                                }

                                if ( cy_sam_tag !~ /^cy:Z:/ ) {
                                    print "Error: Expected \"cy:Z:\" instead of \"" substr(cy_sam_tag, 1, 5) "\" for uncorrected barcode quality with read name (\"" read_name_corrected_bc "\") (line number: " corrected_bc_line_number ").";
                                    exit(1);
                                }
                                if ( cr_sam_tag !~ /^cr:Z:/ ) {
                                    print "Error: Expected \"tb:Z:\" instead of \"" substr(cr_sam_tag, 1, 5) "\" for corrected tagmentation barcode with read name (\"" read_name_corrected_bc "\") (line number: " corrected_bc_line_number ").";
                                    exit(1);
                                }

                                # Extract 10x/HyDrop ATAC sequence and quality and ScaleBio ATAC tagmentation barcode.
                                sequence_CB = substr(cr_sam_tag, 6);
                                quality_CB = substr(cy_sam_tag, 6);
                                tagmentation_barcode = substr(tb_sam_tag, 6);
                            } else {
                                # No corrected tagmentation barcode found.
                                tagmentation_barcode = "";
                            }
                        }
                    }
                } else if ( fastq_part == 2 ) {
                    # Sequence lines.

                    # Store sequence info from R1 and R2 for later use.
                    sequence_R1 = fastq_R1_line;
                    sequence_R2 = fastq_R2_line;
                } else if ( fastq_part == 0 ) {
                    # Quality lines.

                    # Only write read if a corrected tagmentation barcode was found.
                    if ( tagmentation_barcode != "" ) {
                        # Keep track of seen tagmentation barcodes, so filehandles can be closed later.
                        tagmentation_barcodes[tagmentation_barcode] += 1;

                        # Write the full FASTQ record to the R1 and R2 output FASTQ file with barcode info in the read name comments.
                        # When write_fastq_R1_cmd and write_fastq_R2_cmd are the same, an interleaved FASTQ fille will be written.
                        print "@" read_name_R1 "\n" sequence_R1 "\n+\n" fastq_R1_line | partial_write_fastq_cmd "_" tagmentation_barcode fastq_output_R1_suffix;
                        print "@" read_name_corrected_bc "\n" sequence_CB "\n+\n" quality_CB | partial_write_fastq_cmd "_" tagmentation_barcode fastq_output_R2_suffix;
                        print "@" read_name_R2 "\n" sequence_R2 "\n+\n" fastq_R2_line | partial_write_fastq_cmd "_" tagmentation_barcode fastq_output_R3_suffix;
                    }
                }
            }
        }

        # Close open file handles.
        close(read_fastq_R1_cmd);
        close(read_fastq_R2_cmd);
        close(read_corrected_bc_file_cmd);

        for ( tagmentation_barcode in tagmentation_barcodes ) {
            close(partial_write_fastq_cmd "_" tagmentation_barcode fastq_output_R1_suffix);
            close(partial_write_fastq_cmd "_" tagmentation_barcode fastq_output_R2_suffix);
            close(partial_write_fastq_cmd "_" tagmentation_barcode fastq_output_R3_suffix);
        }
    }'

    return $?
}



split_scalebio_atac_fastqs "${@}";

