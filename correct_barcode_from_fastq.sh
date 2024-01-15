#!/bin/bash
#
# Copyright (C) 2020-2024 - Gert Hulselmans
#
# Purpose:
#   Read index 2 FASTQ file with raw barcode reads and correct them
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


SEQC_RUN="seqc run -release";
CODON_RUN="codon run -plugin seq -release";
CODON_OR_SEQ_RUN="${CODON_OR_SEQ_RUN:-${CODON_RUN}}";



correct_barcode_from_fastq () {
    local bc_whitelist_filename="${1}";
    local bc_remapping_filename="${2}";
    local fastq_with_raw_bc_filename="${3}";
    local corrected_bc_filename="${4}";
    local bc_suffix="${5:-1}";
    local max_mismatches="${6:-1}";
    local min_frac_bcs_to_find="${7:-0.5}";

    if [ ${#@} -lt 4 ] ; then
        printf 'Usage:\n';
        printf '    correct_barcode_from_fastq \\\n';
        printf '        bc_whitelist_file bc_remapping_file fastq_with_raw_bc_file \\\n';
        printf '        corrected_bc_file [bc_suffix] [max_mismatches] [min_frac_bcs_to_find]\n\n';
        printf 'Purpose:\n';
        printf '    Read index 2 FASTQ file with raw barcode reads and correct them\n';
        printf '    according to the provided whitelist.\n\n';
        printf '    For each FASTQ record, return:\n';
        printf '      - read name\n';
        printf '      - raw barcode sequences (CR:Z:bc_raw_seq)\n';
        printf '      - raw barcode qualities (CY:Z:bc_raw_qual)\n';
        printf '      - corrected barcode sequences (CB:z:bc_corrected_seq-bc_suffix)\n';
        printf '        (if raw barcode sequence was correctable.)\n\n';
        printf 'Arguments:\n';
        printf '    bc_whitelist_file:\n';
        printf '        File with barcode whitelist to use to correct raw barcode sequences.\n';
        printf '    bc_remapping_file:\n';
        printf '        File with barcodes to use for mapping corrected barcodes to other\n';
        printf '        barcodes (e.g. map 10x multiome ATAC barcodes to 10X multiome RNA\n';
        printf '        barcodes). Set to "false" or "none" to disable remapping.\n';
        printf '    fastq_with_raw_bc_file:\n';
        printf '        FASTQ index 2 file with raw barcode reads to correct.\n';
        printf '    corrected_bc_file:\n';
        printf '        Output file with read name, raw barcode sequence, raw barcode quality\n';
        printf '        and corrected barcode sequence (if correctable) for each FASTQ record\n';
        printf '        in FASTQ index 2 file.\n';
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

    if [ ! -e "${bc_whitelist_filename}" ] ; then
        printf 'Error: Barcode whitelist file "%s" could not be found.\n' "${bc_whitelist_filename}" >&2;
        return 1;
    fi

    case "${bc_remapping_filename,,}" in
        ""|false|none)
            ;;
        *)
            if [ ! -e "${bc_remapping_filename}" ] ; then
                printf 'Error: Barcode remapping file "%s" could not be found.\n' "${bc_remapping_filename}" >&2;
                return 1;
            fi;;
    esac

    if [ ! -e "${fastq_with_raw_bc_filename}" ] ; then
        printf 'Error: FASTQ file with raw barcodes "%s" could not be found.\n' "${fastq_with_raw_bc_filename}" >&2;
        return 1;
    fi

    local first_barcode='';

    # Read first barcode line from barcode whitelist file.
    if [ "${bc_whitelist_filename%.gz}" == "${bc_whitelist_filename}" ] ; then
        # Uncompressed file.
        read -r first_barcode < "${bc_whitelist_filename}";
    else
        # Unset pipefail.
        set +o pipefail

        # Gzip compressed file.
        first_barcode=$(zcat "${bc_whitelist_filename}" | head -n 1);

        # Set pipefail.
        set -o pipefail
    fi

    # Only get the first barcode of the first line in case there are multiple columns like in HyDrop barcode files.
    first_barcode="${first_barcode%%$'\t'*}";

    # Get length of first barcode.
    local -i bc_length="${#first_barcode}";

    # Check if the first barcode is not empty and only contains ACGT and no other characters.
    if [[ ${bc_length} -eq 0 || "${first_barcode}" != "${first_barcode//[^ACGT]/}" ]] ; then
        printf 'Error: The first line of "%s" does not contain a valid barcode.\n' "${bc_whitelist_filename}";
        return 1;
    fi

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Correct barcodes in index 2 FASTQ file:
    #   - Replace "type K = Kmer[16]" with the correct barcode length in the correct_barcode_from_fastq.seq script.
    ${CODON_OR_SEQ_RUN} \
        -D bc_length="${bc_length}" \
        "${script_dir}/correct_barcode_from_fastq.seq" \
            "${bc_whitelist_filename}" \
            "${bc_remapping_filename}" \
            "${fastq_with_raw_bc_filename}" \
            "/dev/stdout" \
            "${corrected_bc_filename%.zst}.corrected_bc_stats.tsv" \
            "${bc_suffix}" \
            "${max_mismatches}" \
            "${min_frac_bcs_to_find}" \
      | zstd -6 -T4 -q -f -o "${corrected_bc_filename%.zst}.zst";
}



correct_barcode_from_fastq "${@}";
