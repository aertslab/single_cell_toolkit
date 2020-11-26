#!/bin/bash


correct_barcode_in_fastq () {
    local bc_whitelist_filename="${1}";
    local fastq_with_raw_bc_filename="${2}";
    local fastq_with_corrected_bc_filename="${3}";

    if [ ${#@} -ne 3 ] ; then
        printf 'Usage: correct_barcode_in_fastq bc_whitelist_file fastq_with_raw_bc_file fastq_with_corrected_bc_file\n';
        return 1;
    fi

    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Correct barcodes in FASTQ file.
    "${script_dir}/run_seq_program.sh" \
        "${script_dir}/correct_barcode_in_fastq.seq" \
        "${bc_whitelist_filename}" \
        "${fastq_with_raw_bc_filename}" \
        "/dev/stdout" \
        "${fastq_with_corrected_bc_filename}.corrected.bc_stats.tsv" \
      | pigz -p 4 \
      > "${fastq_with_corrected_bc_filename}";
}



correct_barcode_in_fastq "${@}";
