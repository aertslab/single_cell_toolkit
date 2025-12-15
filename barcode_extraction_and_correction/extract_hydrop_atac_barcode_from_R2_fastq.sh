#!/bin/bash
#
# Copyright (C) 2021-2024 - Gert Hulselmans
#
# Purpose:
#   Extract HyDrop ATAC barcodes from R2 fastq read and write a new FASTQ
#   file which contains only the barcode.


set -e
set -o pipefail


decompress_fastq_cat_cmd='cat';
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



extract_hydrop_atac_barcode_from_R2_fastq () {
    local fastq_R2_filename="${1}";
    local fastq_R2_BCONLY_filename="${2}";
    local hydrop_atac_barcode_design="${3:-3x96_5L}";
    local compress_fastq_cmd="${4:-bgzip}";

    if [ ${#@} -lt 2 ] ; then
        printf '\nUsage:\n';
        printf '    extract_hydrop_atac_barcode_from_R2_fastq \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_R2_BCONLY \\\n';
        printf '        <hydrop_atac_barcode_design [3x96_5L|v2|2x384_10L|v1|3x96_10L|v0|10X]> \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|igzip|gzip|stdout|-|uncompressed]> \\\n\n';
        printf 'Purpose: Extract HyDrop ATAC barcodes from R2 fastq read and write a new FASTQ\n';
        printf '         file which contains only the barcode.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R2:   FASTQ R2 filename with barcodes (uncompressed or gzipped).\n';
        printf '  - fastq_R2_BCONLY: Output FASTQ R2 filename with extracted HyDrop ATAC barcodes.\n';
        printf '  - hydrop_atac_barcode_design: "3x96_5L" = "v2" (default) or "2x384_10L" = "v1" or "3x96_10L" = "v0" or "10X".\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'  (default)\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'\n" "${compress_fastq_pigz_cmd}";
        printf "          - \"igzip\":  '%s'  (very fast, low compression)\n" "${compress_fastq_igzip_cmd}";
        printf "          - \"gzip\":   '%s'\n" "${compress_fastq_gzip_cmd}";
        printf '          - "stdout":  Write uncompressed output to stdout.\n';
        printf '          - "-":       Write uncompressed output to stdout.\n';
        printf '          - "uncompressed":  Write uncompressed FASTQ files.\n';
        printf '          - full custom command\n\n';
        printf '        To change number of compression threads:\n';
        printf '          - export compress_fastq_threads="%s"\n\n' "${compress_fastq_threads}";
        printf '        To change compression level:\n';
        printf '          - export compress_fastq_level="%s"\n\n' "${compress_fastq_level}";
        return 1;
    fi


    case "${hydrop_atac_barcode_design}" in
        '3x96_10L'|'v0')
            # 10_BC1-10_L-10_BC2-10_L-10_BC3
            hydrop_atac_barcode_design='3x96_10L';
            hydrop_atac_barcode_splits=3;
            local split1_start=1;
            local split1_end=10;
            local split2_start=21;
            local split2_end=30;
            local split3_start=41;
            local split3_end=50;
            ;;
        '2x384_10L'|'v1')
            # 10_BC1-10_L-10_BC2
            hydrop_atac_barcode_design='2x384_10L';
            hydrop_atac_barcode_splits=2;
            local split1_start=1;
            local split1_end=10;
            local split2_start=21;
            local split2_end=30;
            # Not used:
            local split3_start=41;
            local split3_end=50;
            ;;
        '3x96_5L'|'v2')
            # 10_BC1-5_L-10_BC2-5_L-10_BC3
            hydrop_atac_barcode_design='3x96_5L';
            hydrop_atac_barcode_splits=3;
            local split1_start=1;
            local split1_end=10;
            local split2_start=16;
            local split2_end=25;
            local split3_start=31;
            local split3_end=40;
            ;;
        '10X')
            # 10X ATAC BC.
            hydrop_atac_barcode_design='10X';
            hydrop_atac_barcode_splits=1;
            local split1_start=1;
            local split1_end=16;
            ;;
        *)
            printf 'Error: Unsupported HyDrop ATAC barcode design "%s".\n' "${hydrop_atac_barcode_design}";
            printf '       Choose: "3x96_5L" = "v2" (default) or "2x384_10L" = "v1" or "3x96_10L" = "v0".\n';
            return 1;;
    esac


    if type igzip > /dev/null 2>&1 ; then
        # Decompress gzipped FASTQ files with igzip if installed (6x faster than gzip).
        local decompress_fastq_gzipped_cmd="${decompress_fastq_igzip_cmd}";
    else
        # Decompress gzipped FASTQ files with gzip.
        local decompress_fastq_gzipped_cmd="${decompress_fastq_zcat_cmd}";
    fi


    # Detect if input FASTQ files are gzip compressed or not.
    if [ "${fastq_R2_filename}" != "${fastq_R2_filename%.gz}" ] ; then
        local decompress_R2_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    else
        local decompress_R2_fastq_cmd="${decompress_fastq_cat_cmd}";
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
        stdout|-)
            # Write interleaved FASTQ files when writing to stdout.
            fastq_R2_BCONLY_filename='/dev/stdout';
            local compress_fastq_cmd="cat";;
        uncompressed)
            # Write uncompressed FASTQ files.
            fastq_R2_BCONLY_filename="${fastq_R2_BCONLY_filename%.gz}";
            local compress_fastq_cmd="cat";;
    esac


    if ! type mawk > /dev/null 2>&1 ; then
        printf 'Error: "mawk" not found or executable.\n';
        return 1;
    fi

    if ! type "${compress_fastq_cmd%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${compress_fastq_cmd%% *}";
        return 1;
    fi


    ${decompress_R2_fastq_cmd} "${fastq_R2_filename}" \
      | mawk \
            -v "split1_start=${split1_start}" \
            -v "split1_end=${split1_end}" \
            -v "split2_start=${split2_start}" \
            -v "split2_end=${split2_end}" \
            -v "split3_start=${split3_start}" \
            -v "split3_end=${split3_end}" \
            -v "hydrop_atac_barcode_splits=${hydrop_atac_barcode_splits}" \
            -v "hydrop_atac_barcode_design=${hydrop_atac_barcode_design}" \
            '
            BEGIN {
                split1_length = split1_end - split1_start + 1;
                if (hydrop_atac_barcode_splits >= 2) {
                    split2_length = split2_end - split2_start + 1;
                }
                if (hydrop_atac_barcode_splits == 3) {
                    split3_length = split3_end - split3_start + 1;
                }
            }
            {
                if (NR % 2 == 1) {
                    # Read name or "+" line.
                    print $0;
                } else {
                    if (hydrop_atac_barcode_splits == 2) {
                        # Extract HyDrop ATAC barcode info from sequence or quality line for 2x384_10L design.
                        print substr($0, split1_start, split1_length) substr($0, split2_start, split2_length);
                    } else if (hydrop_atac_barcode_splits == 3) {
                        # Extract HyDrop ATAC barcode info from sequence or quality line for 3x96_5L or 3x96_10L design.
                        print substr($0, split1_start, split1_length) substr($0, split2_start, split2_length) substr($0, split3_start, split3_length);
                    } else if (hydrop_atac_barcode_splits == 1) {
                        # Extract 10X ATAC barcode info from sequence or quality line.
                        print substr($0, split1_start, split1_length);
                    }
                }
            }' \
      | ${compress_fastq_cmd} \
      > "${fastq_R2_BCONLY_filename}";

    return $?
}



extract_hydrop_atac_barcode_from_R2_fastq "${@}";

