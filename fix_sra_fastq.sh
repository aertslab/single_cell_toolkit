#!/bin/bash
#
# Copyright (C) 2021 - Gert Hulselmans
#
# Purpose:
#   Fix read names of FASTQ files generated with SRA Toolkit.
#   If an Illumina read name is specified as first comment in the FASTQ
#   comment field, it will be used as new read name, else the SRR read
#   name will be used (without comment fields).


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



fix_sra_fastq () {
    local fastq_input_filename="${1}";
    local fastq_output_filename="${2}";
    local compress_fastq_cmd="${3:-pigz}";

    if [ ${#@} -lt 2 ] ; then
        printf '\nUsage:\n';
        printf '    fix_sra_fastq \\\n';
        printf '        fastq_input \\\n';
        printf '        fastq_output \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|gzip|stdout|-|uncompressed]> \\\n\n';
        printf 'Purpose: Fix read names of FASTQ files generated with SRA Toolkit.\n';
        printf '         If an Illumina read name is specified as first comment in the FASTQ\n';
        printf '         comment field, it will be used as new read name, else the SRR read\n';
        printf '         name will be used (without comment fields).\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_input:   FASTQ input filename created with SRA toolkit (uncompressed or gzipped).\n';
        printf '  - fastq_output:  FASTQ output filename with fixed read name (uncompressed or gzipped).\n';
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


    if type igzip > /dev/null 2>&1 ; then
        # Decompress gzipped FASTQ files with igzip if installed (6x faster than gzip).
        local decompress_fastq_gzipped_cmd="${decompress_fastq_igzip_cmd}";
    else
        # Decompress gzipped FASTQ files with gzip.
        local decompress_fastq_gzipped_cmd="${decompress_fastq_zcat_cmd}";
    fi


    # Detect if input FASTQ files are gzip compressed or not.
    if [ "${fastq_input_filename}" != "${fastq_input_filename%.gz}" ] ; then
        local decompress_input_fastq_cmd="${decompress_fastq_gzipped_cmd}";
    else
        local decompress_input_fastq_cmd="${decompress_fastq_cat_cmd}";
    fi

    if [ "${fastq_output_filename}" == "${fastq_output_filename%.gz}" ] ; then
        case "${fastq_output_filename}" in
            -|stdout|/dev/stdout)
                # Write uncompressed FASTQ file to stdout.
                fastq_output_filename='/dev/stdout';
                local compress_fastq_cmd="cat";;
        esac
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
            # Write uncompressed FASTQ file to stdout.
            fastq_output_filename='/dev/stdout';
            local compress_fastq_cmd="cat";;
        uncompressed)
            # Write uncompressed FASTQ file.
            fastq_output_filename="${fastq_R1_output_filename%.gz}";
            local compress_fastq_cmd="cat";;
    esac


    mawk \
        -v fastq_input_filename="${fastq_input_filename}" \
        -v fastq_output_filename="${fastq_output_filename}" \
        -v decompress_input_fastq_cmd="${decompress_input_fastq_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
        -F ' ' \
    '
    BEGIN {
        read_fastq_input_cmd = decompress_input_fastq_cmd " " fastq_input_filename;

        write_fastq_output_cmd = compress_fastq_cmd " > " fastq_output_filename;

        fastq_line_number = 0;

        # Read FASTQ input file.
        while ( (read_fastq_input_cmd | getline) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            if ( fastq_part == 1 ) {
                # Read name lines with comments.

                if ( $1 ~ /^@SRR/ ) {
                    # FASTQ file was generated with SRA Toolkit (fastq-dump or fasterq-dump).
                    if (NF >= 3 && $2 ~ /^[A-Za-z0-9]+:[A-Za-z0-9]+:/ ) {
                        # Assume it is the original read name given by the Illumina sequencer.
                        read_name = $2;
                    } else {
                        # Use the read name given by SRA Toolkit.
                        read_name = substr($1, 2);
                    }
                } else {
                    # Just keep the FASTQ name and trow away the comment.
                    read_name = substr($1, 2);
                }
            } else if ( fastq_part == 2 ) {
                # Sequence lines.

                # Store sequence line.
                sequence = $0;
            } else if ( fastq_part == 0 ) {
                # Quality lines.

                # Write the full FASTQ record to output FASTQ file.
                print "@" read_name "\n" sequence "\n+\n" $0 | write_fastq_output_cmd;
            }
        }

        # Close open file handles.
        close(read_fastq_input_cmd);
        close(write_fastq_output_cmd);
    }'

    return $?
}



fix_sra_fastq "${@}";
