#!/bin/bash
#
# Copyright (C) 2022-2024 - Gert Hulselmans
#
# Purpose:
#   Get first X number of reads from FASTQ file.


set -e
set -o pipefail



fastq_subset () {
    local fastq_filename="${1}";
    local output_dir="${2}";
    local -i no_reads="${3}";

    if [ ${#@} -ne 3 ] ; then
        printf 'Usage:   fastq_subset fastq_filename output_dir no_reads\n\n';
        printf 'Purpose: Get first X number of reads from FASTQ file.\n\n';

        return 1;
    fi

    if [ ! -e "${fastq_filename}" ] ; then
        printf 'Error: FASTQ filename "%s" does not exist.\n' "${fastq_filename}";
        return 1;
    fi

    if [ ! -d "${output_dir}" ] ; then
        printf 'Error: Output dir "%s" does not exist.\n' "${output_dir}";
        return 1;
    fi

    if [ $(dirname $(realpath "${fastq_filename}")) == $(dirname $(realpath "${output_dir}")) ] ; then
        printf 'Error: Specified output dir "%s" would cause overwriting of FASTQ file "%s".\n' "${output_filename}" "${fastq_filename}";
        return 1;
    fi

    if ! type igzip > /dev/null 2>&1 ; then
         printf 'Error: "igzip" (from ISA-L) is not installed.\n' >&2;
         return 1;
    fi

    if ! type bgzip > /dev/null 2>&1 ; then
         printf 'Error: "bgzip" (from HTSlib) is not installed.\n' >&2;
         return 1;
    fi

    igzip -cd "${fastq_filename}" \
      | head -n $((4*${no_reads})) \
      | bgzip -@ 4 -c \
      > "${output_dir}/$(basename "${fastq_filename}")";
}



fastq_subset "${@}";
