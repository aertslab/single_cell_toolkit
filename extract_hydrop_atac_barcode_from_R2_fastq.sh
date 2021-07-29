#!/bin/bash
#
# Copyright (C) 2021 - Gert Hulselmans
#
# Purpose:
#   Extract HyDrop ATAC barcodes from R2 fastq read and write a new FASTQ
#   file which contains only the barcode.


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
    local hydrop_atac_barcode_design="${3:-2x384}";
    local compress_fastq_cmd="${4:-pigz}";

    local split1_start=1;
    local split1_end=10;
    local split2_start=21;
    local split2_end=30;
    local split3_start=41;
    local split3_end=50;

    if [ ${#@} -lt 2 ] ; then
        printf '\nUsage:\n';
        printf '    extract_hydrop_atac_barcode_from_R2_fastq \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_R2_BCONLY \\\n';
        printf '        <hydrop_atac_barcode_design [3x96|2x384]> \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|igzip|gzip|stdout|-|uncompressed]> \\\n\n';
        printf 'Purpose: Extract HyDrop ATAC barcodes from R2 fastq read and write a new FASTQ\n';
        printf '         file which contains only the barcode.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R2:   FASTQ R2 filename with barcodes (uncompressed or gzipped).\n';
        printf '  - fastq_R2_BCONLY: Output FASTQ R2 filename with extracted HyDrop ATAC barcodes.\n';
        printf '  - hydrop_atac_barcode_design: "3x96" or "2x384" (default).\n';
        printf '  - compress_fastq_cmd:\n';
        printf '      - Compression program to use for output FASTQ files:\n';
        printf "          - \"bgzip\":  '%s'\n" "${compress_fastq_bgzip_cmd}";
        printf "          - \"pigz\":   '%s'  (default)\n" "${compress_fastq_pigz_cmd}";
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
        '3x96')
            ;;
        '2x384')
            ;;
        *)
            printf 'Error: Unsupported HyDrop ATAC barcode design "%s". Choose: "3x96" or "2x384".' "${hydrop_atac_barcode_design}";
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


    ${decompress_R2_fastq_cmd} "${fastq_R2_filename}" \
      | mawk \
            -v "split1_start=${split1_start}" \
            -v "split1_end=${split1_end}" \
            -v "split2_start=${split2_start}" \
            -v "split2_end=${split2_end}" \
            -v "split3_start=${split3_start}" \
            -v "split3_end=${split3_end}" \
            -v "hydrop_atac_barcode_design=${hydrop_atac_barcode_design}" \
            '
            BEGIN {
                split1_length = split1_end - split1_start + 1;
                split2_length = split2_end - split2_start + 1;
                split3_length = split3_end - split3_start + 1;

                if (hydrop_atac_barcode_design == "3x96") {
                    hydrop_atac_barcode_splits = 3;
                } else if (hydrop_atac_barcode_design == "2x384") {
                    hydrop_atac_barcode_splits = 2;
                }
            }
            {
                if (NR % 2 == 1) {
                    # Read name or "+" line.
                    print $0;
                } else {
                    if (hydrop_atac_barcode_splits == 2) {
                        # Extract HyDrop ATAC barcode info from sequence or quality line for 2x384 design.
                        print substr($0, split1_start, split1_length) substr($0, split2_start, split2_length);
                    } else if (hydrop_atac_barcode_splits == 3) {
                        # Extract HyDrop ATAC barcode info from sequence or quality line for 3x96 design.
                        print substr($0, split1_start, split1_length) substr($0, split2_start, split2_length) substr($0, split3_start, split3_length);
                    }
                }
            }' \
      | ${compress_fastq_cmd} \
      > "${fastq_R2_BCONLY_filename}";

    return $?
}



extract_hydrop_atac_barcode_from_R2_fastq "${@}";

