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

    # Create a FIFO file for the FASTQ output, so correct_barcode_in_fastq_fast.seq can write data
    local fastq_with_corrected_bc_fifo_filename="${fastq_with_corrected_bc_filename%.gz}_fifo";
    rm -f "${fastq_with_corrected_bc_fifo_filename}";
    mkfifo "${fastq_with_corrected_bc_fifo_filename}";


    "${script_dir}/run_seq_program.sh" \
        "${script_dir}/correct_barcode_in_fastq_fast.seq" \
        "${bc_whitelist_filename}" \
        "${fastq_with_raw_bc_filename}" \
        "/dev/stdout" \
        "${fastq_with_corrected_bc_filename}.corrected.bc_stats.tsv" \
      | pigz -p 4 \
      > "${fastq_with_corrected_bc_filename}";
}



correct_barcode_in_fastq "${@}";

