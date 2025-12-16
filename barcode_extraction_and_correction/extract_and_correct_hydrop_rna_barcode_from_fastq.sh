#!/bin/bash
#
# Copyright (C) 2020-2025 - Gert Hulselmans
#
# Purpose:
#   Read HyDrop RNA read 1 FASTQ file with raw barcode reads and correct them
#   according to the provided whitelist.
#
#   For each FASTQ record, return:
#     - read name
#     - raw barcode sequences (CR:Z:bc_raw_seq)
#     - raw barcode qualities (CY:Z:bc_raw_qual)
#     - corrected barcode sequences (CB:z:bc_corrected_seq-bc_suffix)
#       (if raw barcode sequence was correctable.)


set -e
set -o pipefail


CODON_RUN_DEFAULT_CMD="codon run -plugin seq -release";
CODON_RUN_CMD="${CODON_RUN_CMD:-${CODON_RUN_DEFAULT_CMD}}";



extract_and_correct_hydrop_rna_barcode_from_fastq () {
    local hydrop_rna_bc_whitelist_filename="${1}";
    local fastq_with_raw_bc_filename="${2}";
    local corrected_bc_filename="${3}";
    local bc_suffix="${4:-1}";
    local max_mismatches="${5:-1}";
    local min_frac_bcs_to_find="${6:-0.5}";

    if [ ${#@} -lt 3 ] ; then
        printf 'Usage:\n';
        printf '    extract_and_correct_hydrop_rna_barcode_from_fastq \\\n';
        printf '        hydrop_rna_bc_whitelist_file fastq_with_raw_bc_file corrected_bc_file \\\n';
        printf '        [bc_suffix] [max_mismatches] [min_frac_bcs_to_find]\n\n';
        printf 'Purpose:\n';
        printf '    Read read 1 FASTQ file with raw HyDrop RNA barcode reads and correct them\n';
        printf '    according to the provided whitelist.\n\n';
        printf '    For each FASTQ record, return:\n';
        printf '      - read name\n';
        printf '      - raw barcode sequences (CR:Z:bc_raw_seq)\n';
        printf '      - raw barcode qualities (CY:Z:bc_raw_qual)\n';
        printf '      - corrected barcode sequences (CB:z:bc_corrected_seq-bc_suffix)\n';
        printf '        (if raw barcode sequence was correctable.)\n\n';
        printf 'Arguments:\n';
        printf '    hydrop_rna_bc_whitelist_file:\n';
        printf '        File with HyDrop RNA barcode whitelist to use to correct raw barcode sequences.\n';
        printf '    fastq_with_raw_bc_file:\n';
        printf '        FASTQ read 1 file with raw barcode reads to correct.\n';
        printf '    corrected_bc_file:\n';
        printf '        Output file with read name, raw barcode sequence, raw barcode quality\n';
        printf '        and corrected barcode sequence (if correctable) for each FASTQ record\n';
        printf '        in FASTQ read 1 file.\n';
        printf '    bc_suffix:\n';
        printf '        Barcode suffix to add after corrected barcode sequence.\n';
        printf '        Default: "1"\n';
        printf '    max_mismatches\n';
        printf '        Maximum amount of mismatches allowed between raw barcode and whitelist.\n';
        printf '        Default: 1\n';
        printf '    min_frac_bcs_to_find\n';
        printf '        Minimum fraction of reads that need to have a barcode that matches the\n';
        printf '        whitelist.\n';
        printf '        Default: 0.5\n';
        return 1;
    fi

    if [ ! -e "${hydrop_rna_bc_whitelist_filename}" ] ; then
        printf 'Error: HyDrop RNA barcode whitelist file "%s" could not be found.\n' "${hydrop_rna_bc_whitelist_filename}" >&2;
        return 1;
    fi

    if [ ! -e "${fastq_with_raw_bc_filename}" ] ; then
        printf 'Error: FASTQ file with raw barcodes "%s" could not be found.\n' "${fastq_with_raw_bc_filename}" >&2;
        return 1;
    fi


    if ! type "${CODON_RUN_CMD%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${CODON_RUN_CMD%% *}";
        return 1;
    fi

    if ! type zstd > /dev/null 2>&1 ; then
        printf 'Error: "zstd" not found or executable.\n';
        return 1;
    fi


    local first_barcode='';

    # Read first barcode line from barcode whitelist file.
    if [ "${hydrop_rna_bc_whitelist_filename%.gz}" == "${hydrop_rna_bc_whitelist_filename}" ] ; then
        # Uncompressed file.
        read -r first_barcode < "${hydrop_rna_bc_whitelist_filename}";
    else
        # Unset pipefail.
        set +o pipefail

        # Gzip compressed file.
        first_barcode=$(zcat "${hydrop_rna_bc_whitelist_filename}" | head -n 1);

        # Set pipefail.
        set -o pipefail
    fi

    # Only get the first barcode of the first line in case there are multiple columns like in HyDrop barcode files.
    first_barcode="${first_barcode%%$'\t'*}";

    # Get length of first barcode.
    local -i bc_length="${#first_barcode}";

    # Check if the first barcode is not empty and only contains ACGT and no other characters.
    if [[ ${bc_length} -eq 0 || "${first_barcode}" != "${first_barcode//[^ACGT]/}" ]] ; then
        printf 'Error: The first line of "%s" does not contain a valid barcode.\n' "${hydrop_rna_bc_whitelist_filename}";
        return 1;
    fi

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Extract and correct HyDrop RNA barcodes in read 1 FASTQ file:
    ${CODON_RUN_CMD} \
        "${script_dir}/extract_and_correct_hydrop_rna_barcode_from_fastq.codon" \
            "${hydrop_rna_bc_whitelist_filename}" \
            "${fastq_with_raw_bc_filename}" \
            "/dev/stdout" \
            "${corrected_bc_filename%.zst}.corrected_bc_stats.tsv" \
            "${bc_suffix}" \
            "${max_mismatches}" \
            "${min_frac_bcs_to_find}" \
      | zstd -6 -T4 -q -f -o "${corrected_bc_filename%.zst}.zst";
}



extract_and_correct_hydrop_rna_barcode_from_fastq "${@}";
