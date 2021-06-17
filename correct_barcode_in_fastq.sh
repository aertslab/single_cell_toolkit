#!/bin/bash


correct_barcode_in_fastq () {
    local bc_whitelist_filename="${1}";
    local fastq_with_raw_bc_filename="${2}";
    local fastq_with_corrected_bc_filename="${3}";
    local max_mismatches="${4:-1}";

    if [ ${#@} -lt 3 ] ; then
        printf 'Usage: correct_barcode_in_fastq bc_whitelist_file fastq_with_raw_bc_file fastq_with_corrected_bc_file max_mismatches\n';
        return 1;
    fi

    if [ ! -e "${bc_whitelist_filename}" ] ; then
        printf 'Error: Barcode whitelist file "%s" could not be found.\n' "${bc_whitelist_filename}" >&2;
        return 1;
    fi

    if [ ! -e "${fastq_with_raw_bc_filename}" ] ; then
        printf 'Error: FASTQ file with raw barcodes "%s" could not be found.\n' "${fastq_with_raw_bc_filename}" >&2;
        return 1;
    fi

    local first_barcode='';

    # Read first barcode from barcode whitelist file.
    if [ "${bc_whitelist_filename%.gz}" == "{bc_whitelist_filename}" ] ; then
        # Uncompressed file.
        read -r first_barcode < "${bc_whitelist_filename}";
    else
        # Gzip compressed file.
        first_barcode=$(zcat "${bc_whitelist_filename}" | head -n 1);
    fi

    # Get length of first barcode.
    local -i bc_length="${#first_barcode}";

    # Check if the first barcode is not empty and only contains ACGT and no other characters.
    if [[ ${bc_length} -eq 0 || "${first_barcode}" != "${first_barcode//[^ACGT]/}" ]] ; then
        printf 'Error: The first line of "%s" does not contain a valid barcode.\n' "${bc_whitelist_filename}";
        return 1;
    fi

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Correct barcodes in FASTQ file:
    #   - Replace "type K = Kmer[16]" with the correct barcode length in the correct_barcode_in_fastq.seq script.
    seqc run \
        -D bc_length="${bc_length}" \
        -release \
        "${script_dir}/correct_barcode_in_fastq.seq" \
            "${bc_whitelist_filename}" \
            "${fastq_with_raw_bc_filename}" \
            "/dev/stdout" \
            "${fastq_with_corrected_bc_filename}.corrected.bc_stats.tsv" \
            "${max_mismatches}" \
      | pigz -p 4 \
      > "${fastq_with_corrected_bc_filename}";
}



correct_barcode_in_fastq "${@}";
