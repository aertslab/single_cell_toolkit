#!/bin/bash

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
        printf 'Usage: extract_and_correct_scalebio_atac_barcode_from_fastq tenx_atac_bc_whitelist_file fastq_with_raw_bc_file corrected_bc_file [bc_suffix] [max_mismatches] [min_frac_bcs_to_find]\n';
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
