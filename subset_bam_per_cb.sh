#!/bin/bash
#
# Copyright (C) 2024 - Gert Hulselmans
#
# Purpose:
#   Subset BAM file per cell barcode.


set -e
set -o pipefail



subset_bam_file_per_cb_chunk () {
    local bam_file="${1}";
    local barcodes_file="${2}";
    local bam_output_prefix="${3%.}";

    if [ ${#@} -ne 3 ] ; then
        printf 'Usage:\n';
        printf '  subset_bam_file_per_cb_chunk bam_file barcodes_file bam_output_prefix\n\n';
        printf 'Arguments:\n';
        printf '  bam_file:          BAM file to subset per provided cell barcode.\n';
        printf '  barcodes_file:     File with cell barcodes to subset input BAM file.\n';
        printf '  bam_output_prefix: Prefix used for output per CB BAM files.\n';
        return 1;
    fi

    if [ "${bam_file##*.}" != "bam" ] ; then
        printf 'Error: BAM file should have ".bam" extension.\n';
        return 1;
    fi

    if [ ! -f "${bam_file}" ] ; then
        printf 'Error: BAM file "%s" does not exist.\n' "${bam_file}";
        return 1;
    fi

    if [ ! -f "${barcodes_file}" ] ; then
        printf 'Error: Barcodes file "%s" does not exist.\n' "${barcodes_file}";
        return 1;
    fi

    if [ ! -d $(dirname "${bam_output_prefix}") ] ; then
        printf 'Error: Directory "%s" does not exist.\n' "$(dirname "${bam_output_prefix}")";
        return 1;
    fi

    # First extract all requested cell barcodes from BAM file
    # before writing individual per CB BAM files.
    samtools view -@ 4 -h -D CB:"${barcodes_file}" "${bam_file}" \
      | mawk \
        -v bam_output_prefix="${bam_output_prefix}" \
        -F '\t' \
        '
        BEGIN {
            header = "";
            CB_tag_found = 0;
        }
        {
            if ( $0 ~ /^@/ ) {
                # Collect SAM header.
                if ( header == "" ) {
                    header = $0;
                } else {
                    header = header "\n" $0;
                }
            } else {
                CB_tag_found = 0;

                # Search for CB tag.
                for ( i = 12; i <= NF; i++ ) {
                    if ( substr($i, 1, 5) == "CB:Z:" ) {
                        CB = substr($i, 6);
                        CB_tag_found = 1;
                        break;
                    }
                }

                if ( CB_tag_found == 1 ) {
                    if (! (CB in CBs_found) ) {
                        # Write BAM header if this is the first time this CB was seen.
                        print header | "samtools view -O bam -o " bam_output_prefix "." CB ".bam -";
                        CBs_found[CB] = 1;
                    }

                    # Write read to CB specific BAM file.
                    print $0 | "samtools view -O bam -o " bam_output_prefix "." CB ".bam -";
                }
            }
        }
        END {
            # Close all file handles.
            for ( CB in CBs_found ) {
                close("samtools view -O bam -o " bam_output_prefix "." CB ".bam -");
            }
        }
        '
}



subset_bam_file_per_cb () {
    local bam_file="${1}";
    local barcodes_file="${2}";
    local bam_output_prefix="${3%.}";

    # Number of barcodes to process in one invocation of subset_bam_file_per_cb_chunk.
    # This is used to limit the number of open file handles and the number of parallel
    # spawned samtools processes.
    local chunk_size="${4:-1000}";

    if [ ${#@} -lt 3 ] ; then
        printf 'Usage:\n';
        printf '  subset_bam_file_per_cb bam_file barcodes_file bam_output_prefix [chunk_size]\n\n';
        printf 'Arguments:\n';
        printf '  bam_file:          BAM file to subset per provided cell barcode.\n';
        printf '  barcodes_file:     File with cell barcodes to subset input BAM file.\n';
        printf '  bam_output_prefix: Prefix used for output per CB BAM files.\n';
        printf '  chunk_size:        Number of cell barcodes to process simultaniously.\n'
        printf '                     If more cell barcodes are provided, the input BAM file\n';
        printf '                     will be read multiple times. It is recommended to not\n';
        printf '                     set this value too high as the same number of parallel\n';
        printf '                     spawned samtools processes will be created for writing\n';
        printf '                     the per CB BAM files.\n';
        printf '                     Default: 1000.\n';
        return 1;
    fi

    if [ "${bam_file##*.}" != "bam" ] ; then
        printf 'Error: BAM file should have ".bam" extension.\n';
        return 1;
    fi

    if [ ! -f "${bam_file}" ] ; then
        printf 'Error: BAM file "%s" does not exist.\n' "${bam_file}";
        return 1;
    fi

    if [ ! -f "${barcodes_file}" ] ; then
        printf 'Error: Barcodes file "%s" does not exist.\n' "${barcodes_file}";
        return 1;
    fi

    if [ ! -d $(dirname "${bam_output_prefix}") ] ; then
        printf 'Error: Directory "%s" does not exist.\n' "$(dirname "${bam_output_prefix}")";
        return 1;
    fi

    if [ $(cat ${barcodes_file} | wc -l) -le ${chunk_size} ] ; then
        # If barcodes file contains less than ${chunk_size} lines, process all barcodes at once.
        subset_bam_file_per_cb_chunk "${bam_file}" "${barcodes_file}" "${bam_output_prefix}";
        return $?;
    else
        # Split barcodes file in chunks of ${chunk_size} lines.
        local subset_bam_file_per_cb_tmp_dir=$(mktemp -t -d subset_bam_file_per_cb_XXXXXX);
        split -l ${chunk_size} -d "${barcodes_file}" ${subset_bam_file_per_cb_tmp_dir}/barcodes_part

        # Subset BAM file in per chunk of barcodes.
        for barcodes_part_file in ${subset_bam_file_per_cb_tmp_dir}/barcodes_part* ; do
            subset_bam_file_per_cb_chunk "${bam_file}" "${barcodes_part_file}" "${bam_output_prefix}";
        done

        # Remove temporary files.
        rm ${subset_bam_file_per_cb_tmp_dir}/barcodes_part*;
        rmdir "${subset_bam_file_per_cb_tmp_dir}";
    fi
}



subset_bam_file_per_cb "${@}";
