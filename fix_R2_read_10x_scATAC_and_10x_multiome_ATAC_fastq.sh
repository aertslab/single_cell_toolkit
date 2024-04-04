#!/bin/bash
#
# Copyright (C) 2021-2024 - Gert Hulselmans
#
# Purpose:
#   Fix R2 read FASTQ files for demultiplexing runs with 10x scATAC and
#   10x multiome ATAC samples.
#      - 10x scATAC: Take first 16 bp.
#      - 10x multiome ATAC: Skip 8 bp and take next 16 bp.


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



fix_R2_read_10x_scATAC_and_10x_multiome_ATAC_fastq () {
    local fastq_R2_filename="${1}";
    local check_or_fix="${2}";
    local compress_fastq_cmd="${3:-bgzip}";

    if [ ${#@} -lt 2 ] ; then
        printf '\nUsage:\n';
        printf '    fix_R2_read_10x_scATAC_and_10x_multiome_ATAC_fastq \\\n';
        printf '        fastq_R2 \\\n';
        printf '        check|fix \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|igzip|gzip]> \\\n\n';
        printf 'Purpose: Fix R2 read FASTQ files for demultiplexing runs with 10x scATAC\n';
        printf '         and 10x multiome ATAC samples.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R2:   FASTQ R2 filename with 10x or HyDrop ATAC barcodes (gzipped).\n';
        printf '  - check|fix:  Check if given FASTQ file contains 10x or HyDrop ATAC\n';
        printf '                barcodes and optionally fix FASTQ files.\n'
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


    # Detect if input FASTQ files are gzip compressed or not.
    if [ "${fastq_R2_filename}" = "${fastq_R2_filename%.gz}" ] ; then
        printf 'Error: FASTQ file "%s" is not gzipped.\n' "${fastq_R2_filename}";
        return 1;
    fi

    if [ ! -e "${fastq_R2_filename}" ] ; then
        printf 'Error: FASTQ file "%s" does not exist.\n' "${fastq_R2_filename}";
        return 1;
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


    local is_scATAC_or_multiome_ATAC=$(
        mawk \
            -v fastq_R2_filename="${fastq_R2_filename}" \
            -v decompress_fastq_cmd="${decompress_fastq_gzipped_cmd}" \
            -v compress_fastq_cmd="${compress_fastq_cmd}" \
        '
        BEGIN {
            read_fastq_R2_cmd = decompress_fastq_cmd " " fastq_R2_filename;

            fastq_line_number = 0;

            # Read FASTQ input file.
            while ( (read_fastq_R2_cmd | getline) > 0 ) {
                fastq_line_number += 1;
                fastq_part = fastq_line_number % 4;

                if ( fastq_line_number == 100000 ) {
                    # Go to end block after reading 100000 lines.
                    exit;
                }

                if ( fastq_part == 1 ) {
                    if ( $2 !~ /^2:/ ) {
                        # FASTQ file is not R2 read.
                        exit;
                    }
                } else if ( fastq_part == 2 ) {
                    # Sequence lines.
                    seq_length = length($0);

                    if ( seq_length != 24 && seq_length != 30  && seq_length != 40) {
                        # FASTQ file does not need to be corrected.
                        # i5 sequence lengths indicate what was sequenced in the run:
                        #   - 24 bp = At least scATAC + scMultiomeATAC
                        #   - 30 bp = At least one HyDrop ATAC 384 sample
                        #   - 40 bp = At least one HyDrop ATAC 96 ligation sample
                        exit;
                    }

                    if ( substr($0, 17, 8) == "GTGTAGAT" ) {
                        is_10x_scATAC += 1;
                    }

                    if ( substr($0, 1, 8) == "CAGACGCG" ) {
                        is_10x_multiome_ATAC += 1;
                    } else if ( substr($0, 1, 8) == "CAGACGGG" ) {
                        is_10x_multiome_ATAC += 1;
                    }
                }
            }
        }
        END {
            # Number of FASTQ reads read (25000 or less if small FASTQ file).
            nbr_fastq_reads = fastq_line_number / 4;

            # If 60% of the reads matched the exact contant sequence for 10x scATAC
            # or 10x multiome ATAC, assign that single cell technique.
            if ( ( is_10x_scATAC / nbr_fastq_reads ) >= 0.6 ) {
                print "is_10x_scATAC";
            } else if ( ( is_10x_multiome_ATAC / nbr_fastq_reads ) >= 0.6 ) {
                print "is_10x_multiome_ATAC";
            } else {
                printf("is_unknown: %f%% 10x scATAC, %f%% 10x multiome ATAC\n", is_10x_scATAC / nbr_fastq_reads * 100.0, is_10x_multiome_ATAC / nbr_fastq_reads * 100.0);
            }
        }
        '
    );


    if [ "${is_scATAC_or_multiome_ATAC}" == "is_10x_scATAC" ] ; then
        local is_multiome_atac=0;
        printf 'FASTQ "%s" contains R2 reads of 10x scATAC.\n' "${fastq_R2_filename}";
    elif [ "${is_scATAC_or_multiome_ATAC}" == "is_10x_multiome_ATAC" ] ; then
        local is_multiome_atac=1;
        printf 'FASTQ "%s" contains R2 reads of 10x scmultiome ATAC.\n' "${fastq_R2_filename}";
    else
        # Uncomment for debugging.
        #printf 'FASTQ "%s" contains %s\n' "${fastq_R2_filename}" "${is_scATAC_or_multiome_ATAC}";

        # Unknown: Do not modify the R2 read.
        return 0;
    fi

    if [ "${check_or_fix}" != "fix" ] ; then
        # Quit here if "fix" was not requested and just report if the FASTQ files
        # contained ScaleBio ATAC with 10x or HyDrop ATAC barcodes.
        return 0;
    fi

    local fastq_R2_output_filename="${fastq_R2_filename}.fixed.fq.gz";

    mawk \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
        -v decompress_fastq_cmd="${decompress_fastq_gzipped_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
        -v is_multiome_atac="${is_multiome_atac}" \
        -F ' ' \
    '
    BEGIN {
        read_fastq_R2_cmd = decompress_fastq_cmd " " fastq_R2_filename;

        write_fastq_R2_cmd = compress_fastq_cmd " > " fastq_R2_output_filename;

        fastq_line_number = 0;

        # Read FASTQ R2 file (which contains the cell barcodes).
        while ( (read_fastq_R2_cmd | getline) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            if ( fastq_part == 1 ) {
                read_name = substr($0, 2);
            } else if ( fastq_part == 2 ) {
                # Sequence lines.

                # Extract correct part of sequence line depending on the single-cell technique.
                if ( is_multiome_atac == 1 ) {
                    # 10x multiome ATAC: Skip 8bp spacer and take next 16 bp.
                    sequence = substr($0, 9, 16);
                } else {
                    # 10x scATAC: First 16 bp
                    sequence = substr($0, 1, 16);
                }
            } else if ( fastq_part == 0 ) {
                # Quality lines.

                # Extract correct part of quality line depending on the single-cell technique.
                if ( is_multiome_atac == 1 ) {
                    # 10x multiome ATAC: Skip 8bp spacer and take next 16 bp.
                    qual = substr($0, 9, 16);
                } else {
                    # 10x scATAC: First 16 bp
                    qual = substr($0, 1, 16);
                }

                # Write the full FASTQ record to output FASTQ file.
                print "@" read_name "\n" sequence "\n+\n" qual | write_fastq_R2_cmd;
            }
        }

        # Close open file handles.
        close(read_fastq_R2_cmd);
        close(write_fastq_R2_cmd);
    }'

    if [ $? -eq 0 ] ; then
        # Rename fixed R2 FASTQ file to the original FASTQ R2 filename.
        mv "${fastq_R2_output_filename}" "${fastq_R2_filename}";
    fi

    return $?
}



fix_R2_read_10x_scATAC_and_10x_multiome_ATAC_fastq "${@}";
