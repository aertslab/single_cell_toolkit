#!/bin/bash
#
# Copyright (C) 2020-2024 - Gert Hulselmans
#
# Purpose:
#   Create corrected scifi-RNA barcode FASTQ file for usage with STAR solo.
#   A corrected scifi-RNA barcode FASTQ record consists of 13 bp RT barcode,
#   16 bp 10x ATAC barcode and 8 bp UMI.


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


SEQC_RUN="seqc run -release";
CODON_RUN="codon run -plugin seq -release";
CODON_OR_SEQ_RUN="${CODON_OR_SEQ_RUN:-${CODON_RUN}}";



create_scifi_rna_barcode_fastq_for_star_solo () {
    local tenx_atac_bc_whitelist_filename="${1}";
    local fastq_R1_filename="${2}";
    local fastq_I2_filename="${3}";
    local fastq_with_corrected_BC_filename="${4}";
    local max_mismatches="${5:-1}";
    local min_frac_bcs_to_find="${6:-0.5}";
    local compress_fastq_cmd="${7:-bgzip}";

    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:\n';
        printf '    create_scifi_rna_barcode_fastq_for_star_solo \\\n';
        printf '        tenx_atac_bc_whitelist_file \\\n';
        printf '        fastq_R1 \\\n';
        printf '        fastq_I2 \\\n';
        printf '        fastq_with_corrected_BC \\\n';
        printf '        max_mismatches \\\n';
        printf '        min_frac_bcs_to_find \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|gzip|stdout|-|uncompressed]> \\\n\n';
        printf 'Purpose: Create corrected scifi-RNA barcode FASTQ file for usage with STAR solo.\n';
        printf '         A corrected scifi-RNA barcode FASTQ record consists of 13 bp RT barcode,\n';
        printf '         16 bp 10x ATAC barcode and 8 bp UMI.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R1:\n';
        printf '        FASTQ R1 filename which contains 8 bp UMI and 13 bp RT barcode.\n';
        printf '  - fastq_I2:\n';
        printf '        FASTQ I2 filename which contains 16 bp 10x ATAC barcodes.\n';
        printf '  - fastq_BC:\n'
        printf '        FASTQ barcode output file with 13 bp RT barcode, 16 bp 10x ATAC barcode and 8 bp UMI.\n';
        printf '    max_mismatches:\n';
        printf '        Maximum amount of mismatches allowed between raw barcode and whitelists.\n';
        printf '        Default: 1\n';
        printf '    min_frac_bcs_to_find:\n';
        printf '        Minimum fraction of reads that need to have a barcode that matches the\n';
        printf '        whitelist.\n';
        printf '        Default: 0.5\n';
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

    if [ "${fastq_I2_filename}" != "${fastq_I2_filename%.gz}" ] ; then
        local decompress_I2_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    elif [ "${fastq_I2_filename}" != "${fastq_I2_filename%.zst}" ] ; then
        local decompress_I2_fastq_cmd="${decompress_fastq_zstd_cmd}";
    else
        local decompress_I2_fastq_cmd="${decompress_fastq_cat_cmd}";
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
            fastq_BC_output_filename='/dev/stdout';
            local compress_fastq_cmd="cat";;
        uncompressed)
            # Write uncompressed FASTQ files.
            fastq_BC_output_filename="${fastq_BC_output_filename%.gz}";
            local compress_fastq_cmd="cat";;
    esac


    if ! type mawk > /dev/null 2>&1 ; then
        printf 'Error: "mawk" not found or executable.\n';
        return 1;
    fi

    if ! type "${CODON_OR_SEQ_RUN%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${CODON_OR_SEQ_RUN%% *}";
        return 1;
    fi

    if ! type "${compress_fastq_cmd%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${compress_fastq_cmd%% *}";
        return 1;
    fi


    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Extract:
    #   - 8 bp UMI and 13 bp RT barcode from R1 FASTQ file.
    #   - 16 bp 10x ATAC barcode from I2 FASTQ file.
    # and write them in 13 bp RT barcode, 16 bp 10x ATAC barcode and 8 bp UMI order
    # to FASTQ barcode output file after barcode correcting RT and 10x ATAC barcode.
    mawk \
        -v fastq_R1_filename="${fastq_R1_filename}" \
        -v fastq_I2_filename="${fastq_I2_filename}" \
        -v decompress_R1_fastq_cmd="${decompress_R1_fastq_cmd}" \
        -v decompress_I2_fastq_cmd="${decompress_I2_fastq_cmd}" \
    '
    BEGIN {
        read_fastq_R1_cmd = decompress_R1_fastq_cmd " " fastq_R1_filename;
        read_fastq_I2_cmd = decompress_I2_fastq_cmd " " fastq_I2_filename;

        read_name_corrected_bc = "";
        corrected_bc_line = "";

        fastq_line_number = 0;
        corrected_bc_line_number = 0;

        # Read FASTQ R1 file (which contains 8 bp UMI and 13 bp RT barcode).
        while ( (read_fastq_R1_cmd | getline fastq_R1_line) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            # Read FASTQ I2 file (which contains 16 bp 10x ATAC cell barcode).
            if ( (read_fastq_I2_cmd | getline fastq_I2_line) > 0 ) {
                if ( fastq_part == 1 ) {
                    # Find first space position (0 if no comment found) in read name from all input FASTQ files.
                    read_name_R1_space_pos = index(fastq_R1_line, " ");
                    read_name_I2_space_pos = index(fastq_I2_line, " ");

                    # Extract read name from all input FASTQ files.
                    if (read_name_R1_space_pos > 0) {
                        read_name_R1 = substr(fastq_R1_line, 2, read_name_R1_space_pos - 2);
                    } else {
                        read_name_R1 = substr(fastq_R1_line, 2);
                    }

                    if (read_name_I2_space_pos > 0) {
                        read_name_I2 = substr(fastq_I2_line, 2, read_name_I2_space_pos - 2);
                    } else {
                        read_name_I2 = substr(fastq_I2_line, 2);
                    }

                    # Check if read names match between both FASTQ files.
                    if ( read_name_R1 != read_name_I2 ) {
                        print "Error: Read name R1 (\"" read_name_R1 "\") and read name I2 (\"" read_name_I2 "\") are not paired properly (line number: " fastq_line_number ")." > "/dev/stderr";
                        exit(1);
                    }
                } else if ( fastq_part == 2 ) {
                    # Sequence lines.

                    # Store sequence info from R1 and I2 for later use.
                    sequence_R1 = fastq_R1_line;
                    sequence_I2 = fastq_I2_line;
                } else if ( fastq_part == 0 ) {
                    # Quality lines.

                    if ( ( length(sequence_R1) < 21 ) || ( length(fastq_R1_line) < 21 ) ) {
                        print "Error: Read name R1 (\"" read_name_R1 "\") sequence and/or quality is less than 21 bp (line number: " fastq_line_number ")." > "/dev/stderr";
                        exit(1);
                    }

                    if ( ( length(sequence_I2) < 16 ) || ( length(fastq_I2_line) < 16 ) ) {
                        print "Error: Read name I2 (\"" read_name_I2 "\") sequence and/or quality is less than 16 bp (line number: " fastq_line_number ")." > "/dev/stderr";
                        exit(1);
                    }

                    UMI_seq = substr(sequence_R1, 1, 8);
                    UMI_qual = substr(fastq_R1_line, 1, 8);
                    RT_BC_seq = substr(sequence_R1, 9, 13);
                    RT_BC_qual = substr(fastq_R1_line, 9, 13);

                    # Write the RT barcode, 10x scATAC barcode and UMI to FASTQ record.
                    print "@" read_name_R1 "\n" RT_BC_seq sequence_I2 UMI_seq "\n+\n" RT_BC_qual fastq_I2_line UMI_qual;
                }
            }
        }

        # Close open file handles.
        close(read_fastq_R1_cmd);
        close(read_fastq_I2_cmd);
    }' \
    | ${CODON_OR_SEQ_RUN} \
        "${script_dir}/extract_and_correct_scifi_rna_barcode_from_fastq.seq" \
            "${tenx_atac_bc_whitelist_filename}" \
            "/dev/stdin" \
            "/dev/stdout" \
            "${fastq_with_corrected_BC_filename%.gz}.corrected_bc_stats.tsv" \
            "${max_mismatches}" \
            "${min_frac_bcs_to_find}" \
      | ${compress_fastq_cmd} \
      > "${fastq_with_corrected_BC_filename}";

    return $?
}



create_scifi_rna_barcode_fastq_for_star_solo "${@}";
