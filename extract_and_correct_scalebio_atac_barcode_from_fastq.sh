#!/bin/bash
#
# Copyright (C) 2020-2024 - Gert Hulselmans
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


SEQC_RUN="seqc run -release";
CODON_RUN="codon run -plugin seq -release";
CODON_OR_SEQ_RUN="${CODON_OR_SEQ_RUN:-${CODON_RUN}}";



extract_and_correct_scalebio_atac_barcode_from_fastq () {
    local tenx_or_hydrop_atac_bc_whitelist_filename="${1}";
    local fastq_with_raw_bc_filename="${2}";
    local corrected_bc_filename="${3}";
    local scalebio_bc_type="${4}";
    local bc_suffix="${5:-1}";
    local max_mismatches="${6:-1}";
    local min_frac_bcs_to_find="${7:-0.5}";

    if [ ${#@} -lt 4 ] ; then
        printf 'Usage:\n';
        printf '    extract_and_correct_scalebio_atac_barcode_from_fastq \\\n';
        printf '        10x_or_hydrop_atac_bc_whitelist_file fastq_with_raw_bc_file corrected_bc_file \\\n';
        printf '        scalebio_bc_type [bc_suffix] [max_mismatches] [min_frac_bcs_to_find]\n\n';
        printf 'Purpose:\n';
        printf '    Read ScaleBio index 2 FASTQ file with raw barcode reads and correct them\n';
        printf '    according to the provided whitelist.\n\n';
        printf '    For each FASTQ record, return:\n';
        printf '      - read name\n';
        printf '      - raw barcode sequences (CR:Z:tagmentation_bc_seq-10x_or_hydrop_atac_bc_seq)\n';
        printf '      - raw barcode qualities (CY:Z:tagmentation_bc_qual 10x_or_hydrop_atac_bc_qual)\n';
        printf '      - corrected barcode sequences\n';
        printf '        (CB:z:corrected_tagmentation_bc-corrected_10x_or_hydrop_atac_bc-bc_suffix)\n';
        printf '        (if raw barcode sequences were both correctable.)\n\n';
        printf 'Arguments:\n';
        printf '    10x_or_hydrop_atac_bc_whitelist_file:\n';
        printf '        File with 10x or HyDrop ATAC barcode whitelist to use to correct raw\n';
        printf '        barcode sequences.\n';
        printf '    fastq_with_raw_bc_file:\n';
        printf '        FASTQ ScaleBio index 2 file with raw tagmentation and 10x ATAC barcode\n';
        printf '        or HyDrop ATAC barcode reads to correct.\n';
        printf '    corrected_bc_file:\n';
        printf '        Output file with read name, raw barcode sequence, raw barcode quality\n';
        printf '        and corrected barcode sequence (if correctable) for each FASTQ record\n';
        printf '        in FASTQ index 2 file.\n';
        printf '    scalebio_bc_type:\n';
        printf '        ScaleBio tagmentation barcode type set to look for:\n';
        printf '        "scalebio" (original ScaleBio), "scalebioih" (ScaleBio inhouse),\n';
        printf '        "scalebioih6" (ScaleBio inhouse trial) and "all".\n';
        printf '    bc_suffix:\n';
        printf '        Barcode suffix to add after corrected barcode sequence.\n';
        printf '        Default: "1"\n';
        printf '    max_mismatches:\n';
        printf '        Maximum amount of mismatches allowed between raw barcode and whitelists.\n';
        printf '        Default: 1\n';
        printf '    min_frac_bcs_to_find:\n';
        printf '        Minimum fraction of reads that need to have a barcode that matches the\n';
        printf '        whitelist.\n';
        printf '        Default: 0.5\n';
        return 1;
    fi

    if [ ! -e "${tenx_or_hydrop_atac_bc_whitelist_filename}" ] ; then
        printf 'Error: 10x or HyDrop ATAC barcode whitelist file "%s" could not be found.\n' "${tenx_or_hydrop_atac_bc_whitelist_filename}" >&2;
        return 1;
    fi

    if [ ! -e "${fastq_with_raw_bc_filename}" ] ; then
        printf 'Error: FASTQ file with raw barcodes "%s" could not be found.\n' "${fastq_with_raw_bc_filename}" >&2;
        return 1;
    fi

    local first_barcode='';

    # Read first barcode line from barcode whitelist file.
    if [ "${tenx_or_hydrop_atac_bc_whitelist_filename%.gz}" == "${tenx_or_hydrop_atac_bc_whitelist_filename}" ] ; then
        # Uncompressed file.
        read -r first_barcode < "${tenx_or_hydrop_atac_bc_whitelist_filename}";
    else
        # Unset pipefail.
        set +o pipefail

        # Gzip compressed file.
        first_barcode=$(zcat "${tenx_or_hydrop_atac_bc_whitelist_filename}" | head -n 1);

        # Set pipefail.
        set -o pipefail
    fi

    # Only get the first barcode of the first line in case there are multiple columns like in HyDrop barcode files.
    first_barcode="${first_barcode%%$'\t'*}";

    # Get length of first barcode.
    local -i bc_length="${#first_barcode}";

    # Check if the first barcode is not empty and only contains ACGT and no other characters.
    if [[ ${bc_length} -eq 0 || "${first_barcode}" != "${first_barcode//[^ACGT]/}" ]] ; then
        printf 'Error: The first line of "%s" does not contain a valid barcode.\n' \
            "${tenx_or_hydrop_atac_bc_whitelist_filename}";
        return 1;
    fi

    if [[ ${bc_length} -eq 16 ]] ; then
        local tenx_or_hydrop="10x";
    elif [[ ${bc_length} -eq 30 ]] ; then
        local tenx_or_hydrop="hydrop";
    else
        printf 'Error: Barcode whitelist file "%s" does not contain 10x or HyDrop ATAC barcodes.' \
            "${tenx_or_hydrop_atac_bc_whitelist_filename}";
        return 1;
    fi

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Correct tagmentation and 10x ATAC barcodes in index 2 FASTQ file:
    ${CODON_OR_SEQ_RUN} \
        "${script_dir}/extract_and_correct_scalebio_atac_barcode_from_fastq.seq" \
            "${tenx_or_hydrop_atac_bc_whitelist_filename}" \
            "${tenx_or_hydrop}" \
            "${fastq_with_raw_bc_filename}" \
            "/dev/stdout" \
            "${corrected_bc_filename%.zst}.corrected_bc_stats.tsv" \
            "${scalebio_bc_type}" \
            "${bc_suffix}" \
            "${max_mismatches}" \
            "${min_frac_bcs_to_find}" \
      | zstd -6 -T4 -q -f -o "${corrected_bc_filename%.zst}.zst";
}



extract_and_correct_scalebio_atac_barcode_from_fastq "${@}";
