#!/bin/bash
#
# Copyright (C) 2020-2023 - Gert Hulselmans
#
# Purpose:
#   Read ScaleBio index 2 FASTQ file with raw barcode reads and correct them
#   according to the provided whitelist.
#
#   For each FASTQ record, return:
#     - read name
#     - raw barcode sequences (CR:Z:tagmentation_bc_seq-tenx_atac_bc_seq)
#     - raw barcode qualities (CY:Z:tagmentation_bc_qual tenx_atac_bc_qual)
#     - corrected barcode sequences (CB:z:corrected_tagmentation_bc-corrected_10x_bc-bc_suffix)
#       (if raw barcode sequences were both correctable.)


set -e
set -o pipefail



extract_and_correct_scalebio_atac_barcode_from_fastq () {
    local tenx_atac_bc_whitelist_filename="${1}";
    local fastq_with_raw_bc_filename="${2}";
    local corrected_bc_filename="${3}";
    local bc_suffix="${4:-1}";
    local max_mismatches="${5:-1}";
    local min_frac_bcs_to_find="${6:-0.5}";

    if [ ${#@} -lt 3 ] ; then
        printf 'Usage:\n';
        printf '    extract_and_correct_scalebio_atac_barcode_from_fastq \\\n';
        printf '        tenx_atac_bc_whitelist_file fastq_with_raw_bc_file corrected_bc_file \\\n';
        printf '        [bc_suffix] [max_mismatches] [min_frac_bcs_to_find]\n\n';
        printf 'Purpose:\n';
        printf '    Read ScaleBio index 2 FASTQ file with raw barcode reads and correct them\n';
        printf '    according to the provided whitelist.\n\n';
        printf '    For each FASTQ record, return:\n';
        printf '      - read name\n';
        printf '      - raw barcode sequences (CR:Z:tagmentation_bc_seq-tenx_atac_bc_seq)\n';
        printf '      - raw barcode qualities (CY:Z:tagmentation_bc_qual tenx_atac_bc_qual)\n';
        printf '      - corrected barcode sequences\n';
        printf '        (CB:z:corrected_tagmentation_bc-corrected_10x_bc-bc_suffix)\n';
        printf '        (if raw barcode sequences were both correctable.)\n\n';
        printf 'Arguments:\n';
        printf '    tenx_atac_bc_whitelist_file:\n';
        printf '        File with 10x ATAC barcode whitelist to use to correct raw barcode\n';
        printf '        sequences.\n';
        printf '    fastq_with_raw_bc_file:\n';
        printf '        FASTQ ScaleBio index 2 file with raw tagmentation and 10x ATAC barcode\n';
        printf '        reads to correct.\n';
        printf '    corrected_bc_file:\n';
        printf '        Output file with read name, raw barcode sequence, raw barcode quality\n';
        printf '        and corrected barcode sequence (if correctable) for each FASTQ record\n';
        printf '        in FASTQ index 2 file.\n';
        printf '    bc_suffix:\n';
        printf '        Barcode suffix to add after corrected barcode sequence.\n';
        printf '        Default: "1"\n';
        printf '    max_mismatches\n';
        printf '        Maximum amount of mismatches allowed between raw barcode and whitelists.\n';
        printf '        Default: 1\n';
        printf '    min_frac_bcs_to_find\n';
        printf '        Minimum fraction of reads that need to have a barcode that matches the\n';
        printf '        whitelist.\n';
        printf '        Default: 0.5\n';
        return 1;
    fi

    if [ ! -e "${tenx_atac_bc_whitelist_filename}" ] ; then
        printf 'Error: 10x ATAC barcode whitelist file "%s" could not be found.\n' "${tenx_atac_bc_whitelist_filename}" >&2;
        return 1;
    fi

    if [ ! -e "${fastq_with_raw_bc_filename}" ] ; then
        printf 'Error: FASTQ file with raw barcodes "%s" could not be found.\n' "${fastq_with_raw_bc_filename}" >&2;
        return 1;
    fi

    local first_barcode='';

    # Read first barcode from barcode whitelist file.
    if [ "${tenx_atac_bc_whitelist_filename%.gz}" == "${tenx_atac_bc_whitelist_filename}" ] ; then
        # Uncompressed file.
        read -r first_barcode < "${tenx_atac_bc_whitelist_filename}";
    else
        # Unset pipefail.
        set +o pipefail

        # Gzip compressed file.
        first_barcode=$(zcat "${tenx_atac_bc_whitelist_filename}" | head -n 1);

        # Set pipefail.
        set -o pipefail
    fi

    # Get length of first barcode.
    local -i bc_length="${#first_barcode}";

    # Check if the first barcode is not empty and only contains ACGT and no other characters.
    if [[ ${bc_length} -eq 0 || "${first_barcode}" != "${first_barcode//[^ACGT]/}" ]] ; then
        printf 'Error: The first line of "%s" does not contain a valid barcode.\n' "${tenx_atac_bc_whitelist_filename}";
        return 1;
    fi

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Correct tagmentation and 10x ATAC barcodes in index 2 FASTQ file:
    seqc run \
        -release \
        "${script_dir}/extract_and_correct_scalebio_atac_barcode_from_fastq.seq" \
            "${tenx_atac_bc_whitelist_filename}" \
            "${fastq_with_raw_bc_filename}" \
            "/dev/stdout" \
            "${corrected_bc_filename}.corrected_bc_stats.tsv" \
            "${bc_suffix}" \
            "${max_mismatches}" \
            "${min_frac_bcs_to_find}" \
      | pigz -p 4 \
      > "${corrected_bc_filename}";
}



extract_and_correct_scalebio_atac_barcode_from_fastq "${@}";
