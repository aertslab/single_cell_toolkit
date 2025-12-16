#!/bin/bash
#
# Copyright (C) 2020-2025 - Gert Hulselmans
#
# Purpose:
#   Create corrected ISSAAC-seq barcode FASTQ file for usage with STAR solo.
#   A corrected ISSAAC-seq barcode FASTQ record consists of 30 bp barcode,
#   and 10 bp UMI.


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


CODON_RUN_DEFAULT_CMD="codon run -plugin seq -release";
CODON_RUN_CMD="${CODON_RUN_CMD:-${CODON_RUN_DEFAULT_CMD}}";



create_issaac_seq_barcode_fastq_for_star_solo () {
    local hydrop_rna_bc_whitelist_filename="${1}";
    local fastq_R2_filename="${2}";
    local fastq_I2_filename="${3}";
    local fastq_with_corrected_BC_filename="${4}";
    local max_mismatches="${5:-1}";
    local min_frac_bcs_to_find="${6:-0.5}";
    local compress_fastq_cmd="${7:-bgzip}";

    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:\n';
        printf '    create_issaac_seq_barcode_fastq_for_star_solo \\\n';
        printf '        hydrop_rna_bc_whitelist_file \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_I2 \\\n';
        printf '        fastq_with_corrected_BC \\\n';
        printf '        max_mismatches \\\n';
        printf '        min_frac_bcs_to_find \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|gzip|stdout|-|uncompressed]> \\\n\n';
        printf 'Purpose: Create corrected ISSAAC-seq barcode FASTQ file for usage with STAR solo.\n';
        printf '         A corrected ISSAAC-seq barcode FASTQ record consists of 30 bp barcode,\n';
        printf '         and 10 bp UMI.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R2:\n';
        printf '        FASTQ R2 filename which contains 10 bp UMI.\n';
        printf '  - fastq_I2:\n';
        printf '        FASTQ I2 filename which contains 30 bp HyDrop RNA barcodes.\n';
        printf '  - fastq_BC:\n'
        printf '        FASTQ barcode output file with 30 bp barcode and 10 bp UMI.\n';
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
    if [ "${fastq_R2_filename}" != "${fastq_R2_filename%.gz}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    elif [ "${fastq_R2_filename}" != "${fastq_R2_filename%.zst}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_zstd_cmd}";
    else
        local decompress_R2_fastq_cmd="${decompress_fastq_cat_cmd}";
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

    if ! type "${CODON_RUN_CMD%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${CODON_RUN_CMD%% *}";
        return 1;
    fi

    if ! type "${compress_fastq_cmd%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${compress_fastq_cmd%% *}";
        return 1;
    fi


    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Extract:
    #   - 10 bp UMI from R2 FASTQ file.
    #   - 30 bp HyDrop RNA barcode from I2 FASTQ file (3x96_5L HyDrop design).
    # and write them in 30 bp barcode and 10 bp UMI order
    # to FASTQ barcode output file after barcode correcting the HyDrop RNA barcode.
    mawk \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v fastq_I2_filename="${fastq_I2_filename}" \
        -v decompress_R2_fastq_cmd="${decompress_R2_fastq_cmd}" \
        -v decompress_I2_fastq_cmd="${decompress_I2_fastq_cmd}" \
    '
    BEGIN {
        read_fastq_R2_cmd = decompress_R2_fastq_cmd " " fastq_R2_filename;
        read_fastq_I2_cmd = decompress_I2_fastq_cmd " " fastq_I2_filename;

        read_name_corrected_bc = "";
        corrected_bc_line = "";

        fastq_line_number = 0;
        corrected_bc_line_number = 0;

        # Read FASTQ R2 file (which contains 8 bp UMI and 13 bp RT barcode).
        while ( (read_fastq_R2_cmd | getline fastq_R2_line) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            # Read FASTQ I2 file (which contains 30 bp HyDrop RNA cell barcode).
            if ( (read_fastq_I2_cmd | getline fastq_I2_line) > 0 ) {
                if ( fastq_part == 1 ) {
                    # Find first space position (0 if no comment found) in read name from all input FASTQ files.
                    read_name_R2_space_pos = index(fastq_R2_line, " ");
                    read_name_I2_space_pos = index(fastq_I2_line, " ");

                    # Extract read name from all input FASTQ files.
                    if (read_name_R2_space_pos > 0) {
                        read_name_R2 = substr(fastq_R2_line, 2, read_name_R2_space_pos - 2);
                    } else {
                        read_name_R2 = substr(fastq_R2_line, 2);
                    }

                    if (read_name_I2_space_pos > 0) {
                        read_name_I2 = substr(fastq_I2_line, 2, read_name_I2_space_pos - 2);
                    } else {
                        read_name_I2 = substr(fastq_I2_line, 2);
                    }

                    # Check if read names match between both FASTQ files.
                    if ( read_name_R2 != read_name_I2 ) {
                        print "Error: Read name R2 (\"" read_name_R2 "\") and read name I2 (\"" read_name_I2 "\") are not paired properly (line number: " fastq_line_number ")." > "/dev/stderr";
                        exit(1);
                    }
                } else if ( fastq_part == 2 ) {
                    # Sequence lines.

                    # Store sequence info from R2 and I2 for later use.
                    sequence_R2 = fastq_R2_line;
                    sequence_I2 = fastq_I2_line;
                } else if ( fastq_part == 0 ) {
                    # Quality lines.

                    if ( ( length(sequence_R2) < 10 ) || ( length(fastq_R2_line) < 10 ) ) {
                        print "Error: Read name R2 (\"" read_name_R2 "\") sequence and/or quality is less than 10 bp (line number: " fastq_line_number ")." > "/dev/stderr";
                        exit(1);
                    }

                    if ( ( length(sequence_I2) < 40 ) || ( length(fastq_I2_line) < 40 ) ) {
                        print "Error: Read name I2 (\"" read_name_I2 "\") sequence and/or quality is less than 40 bp (line number: " fastq_line_number ")." > "/dev/stderr";
                        exit(1);
                    }

                    UMI_seq = substr(sequence_R2, 1, 10);
                    UMI_qual = substr(fastq_R2_line, 1, 10);

                    # HyDrop RNA barcode design: "3x96_5L" = "v2": 10_BC1-5_L-10_BC2-5_L-10_BC3
                    BC_seq = substr(sequence_I2, 1, 10) substr(sequence_I2, 16, 10) substr(sequence_I2, 31, 10);
                    BC_qual = substr(fastq_I2_line, 1, 10) substr(fastq_R2_line, 16, 10) substr(fastq_I2_line, 31, 10);

                    # Write the HyDrop RNA barcode and UMI to FASTQ record.
                    print "@" read_name_R2 "\n" BC_seq UMI_seq "\n+\n" RT_BC_qual fastq_I2_line UMI_qual;
                }
            }
        }

        # Close open file handles.
        close(read_fastq_R2_cmd);
        close(read_fastq_I2_cmd);
    }' \
    | ${CODON_RUN_CMD} \
        "${script_dir}/extract_and_correct_issaac_seq_barcode_from_fastq.codon" \
            "${hydrop_rna_bc_whitelist_filename}" \
            "/dev/stdin" \
            "/dev/stdout" \
            "${fastq_with_corrected_BC_filename%.gz}.corrected_bc_stats.tsv" \
            "${max_mismatches}" \
            "${min_frac_bcs_to_find}" \
      | ${compress_fastq_cmd} \
      > "${fastq_with_corrected_BC_filename}";

    return $?
}



create_issaac_seq_barcode_fastq_for_star_solo "${@}";
