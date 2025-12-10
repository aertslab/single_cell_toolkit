#!/bin/bash
#
# Copyright (C) 2025 - Gert Hulselmans
#
# Purpose:
#   Convert fragments file to cut sites file.



fragments_to_cut_sites () {
    local fragments_tsv_file="${1}";
    local cut_sites_tsv_file="${2}";

    if [ ${#@} -ne 2 ] ; then
        printf 'Usage:\n';
        printf '    fragments_to_cut_sites fragments_tsv_file cut_sites_tsv_file\n\n';
        printf 'Purpose:\n';
        printf '    Convert fragments file to cut sites file.\n';
        return 1;
    fi

    if ! type bgzip > /dev/null 2>&1 ; then
        printf 'Error: "bgzip" (HTSlib) not found or executable.\n';
        return 1;
    fi

    if ! type igzip > /dev/null 2>&1 ; then
        printf 'Error: "igzip" (ISA-L) not found or executable.\n';
        return 1;
    fi

    if ! type mawk > /dev/null 2>&1 ; then
        printf 'Error: "mawk" not found or executable.\n';
        return 1;
    fi

    # Convert fragments file to cut sites.
    igzip -cd "${fragments_tsv_file}" \
        | mawk \
            -F '\t' \
            -v 'OFS=\t' \
            '
            BEGIN {
                chrom_prev = "";
                sort_cmd = "LC_ALL=C sort -k 2,2n -k 3,3n";
            }
            {
                if ( $0 !~ /^#/ ) {
                    chrom = $1;

                    if ( chrom != chrom_prev ) {
                        close(sort_cmd);
                    }

                    print $1, $2, $2 + 1, $4, $5 | sort_cmd;
                    print $1, $3 - 1, $3, $4, $5 | sort_cmd;

                    chrom_prev = chrom;
                }
            }
            END {
                close(sort_cmd);
            }
            ' \
        | bgzip -@ 4 -c \
        > "${cut_sites_tsv_file}";
}



fragments_to_cut_sites "${@}";
