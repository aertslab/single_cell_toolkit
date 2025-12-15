#!/bin/bash
#
# Copyright (C) 2021-2024 - Gert Hulselmans
#
# Purpose:
#   Limit AVITI FASTQ quality scores to max Q40 as at the moment
#   CellRanger/CellRangerATAC/CellRangerARC have problems with
#   scores higher than Q40.


decompress_fastq_zcat_cmd='zcat';
decompress_fastq_igzip_cmd='igzip -c -d';


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



limit_aviti_fastq_to_max_q40 () {
    local fastq_input_filename="${1}";
    local output_dir="${2}";
    local compress_fastq_cmd="${3:-bgzip}";

    if [ ${#@} -lt 2 ] ; then
        printf '\nUsage:\n';
        printf '    limit_aviti_fastq_to_max_q40 \\\n';
        printf '        fastq \\\n';
        printf '        output_dir \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|igzip|gzip]> \\\n\n';
        printf 'Purpose: Limit AVITI FASTQ quality scores to max Q40, so they can work with\n'
        printf '         CellRanger/CellRangerATAC/CellRangerARC.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq:      FASTQ input filename.\n';
        printf '  - output_dir: Output directory.\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'  (default)\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'\n" "${compress_fastq_pigz_cmd}";
        printf "          - \"igzip\":  '%s'  (very fast, low compression)\n" "${compress_fastq_igzip_cmd}";
        printf "          - \"gzip\":   '%s'\n" "${compress_fastq_gzip_cmd}";
        printf '          - full custom command\n\n';
        printf '        To change number of compression threads:\n';
        printf '          - export compress_fastq_threads="%s"\n\n' "${compress_fastq_threads}";
        printf '        To change compression level:\n';
        printf '          - export compress_fastq_level="%s"\n\n' "${compress_fastq_level}";
        return 1;
    fi

    if [ "${fastq_input_filename}" = "${fastq_input_filename%.gz}" ] ; then
        printf 'Error: FASTQ file "%s" is not gzipped.\n' "${fastq_input_filename}";
        return 1;
    fi

    if [ ! -e "${fastq_input_filename}" ] ; then
        printf 'Error: FASTQ file "%s" does not exist.\n' "${fastq_input_filename}";
        return 1;
    fi

    if [ ! -d "${output_dir}" ] ; then
        printf 'Error: Output dir "%s" does not exist.\n' "${output_dir}";
        return 1;
    fi

    local fastq_output_filename="${output_dir%/}/${fastq_input_filename##.*/}";

    if [ -e "${fastq_output_filename}" ] ; then
        if type realpath > /dev/null 2>&1 ; then
            if [ "$(realpath "${fastq_input_filename}")" = "$(realpath "${fastq_output_filename}")" ] ; then
                printf 'Error: Output FASTQ file "%s" would overwrite input FASTQ file.\n' "${fastq_output_filename}";
                return 1;
            fi
        else
            if [ "${fastq_input_filename}" = "${fastq_output_filename}" ] ; then
                printf 'Error: Output FASTQ file "%s" would overwrite input FASTQ file.\n' "${fastq_output_filename}";
                return 1;
            fi
        fi
    fi

    if type igzip > /dev/null 2>&1 ; then
        # Decompress gzipped FASTQ files with igzip if installed (6x faster than gzip).
        local decompress_fastq_gzipped_cmd="${decompress_fastq_igzip_cmd}";
    else
        # Decompress gzipped FASTQ files with gzip.
        local decompress_fastq_gzipped_cmd="${decompress_fastq_zcat_cmd}";
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
    esac


    if ! type mawk > /dev/null 2>&1 ; then
         printf 'Error: "mawk" is not installed.\n' >&2;
         return 1;
    fi

    if ! type "${compress_fastq_cmd%% *}" > /dev/null 2>&1 ; then
         printf 'Error: "%s" is not installed.\n' "${compress_fastq_cmd%% *}" >&2;
         return 1;
    fi


    ${decompress_fastq_gzipped_cmd} "${fastq_input_filename}" \
      | mawk '
            {
                if (NR % 4 == 0) {
                    # Change all Phred scores higher than 40 to 40.
                    gsub(/[J-Z]/, "I" ,$0);
                }

                print $0;
            }' \
      | ${compress_fastq_cmd} \
      > "${fastq_output_filename}";

    return $?
}



limit_aviti_fastq_to_max_q40 "${@}";
