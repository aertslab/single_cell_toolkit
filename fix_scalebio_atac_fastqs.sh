#!/bin/bash
#
# Copyright (C) 2020-2024 - Gert Hulselmans
#
# Purpose:
#   Fix ScaleBio ATAC R2 and R3 reads after standard demultiplexing so
#   afterwards R2 contains 16 bp 10X ATAC barcode and 8 bp tagmentation
#   barcode or 40 bp HyDrop ATAC barcode and 8 bp tagmentation barcode
#   and R3 contains only the second ATAC read.


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



fix_scalebio_atac_fastqs () {
    local fastq_R2_filename="${1}";
    local fastq_R3_filename="${2}";
    local check_or_fix="${3}";
    local compress_fastq_cmd="${4:-bgzip}";

    if [ ${#@} -lt 3 ] ; then
        printf '\nUsage:\n';
        printf '    fix_scalebio_atac_fastqs \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_R3 \\\n';
        printf '        check|fix \\\n';
        printf '        <compress_fastq_cmd [bgzip|pigz|gzip|stdout|-|uncompressed]> \\\n\n';
        printf 'Purpose: Fix ScaleBio ATAC R2 and R3 reads after standard demultiplexing so\n';
        printf '         afterwards R2 contains 16 bp 10X ATAC barcode and 8 bp tagmentation\n';
        printf '         barcode or 40 bp HyDrop ATAC barcode and 8 bp tagmentation barcode\n';
        printf '         and R3 contains only the second ATAC read.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R2:   FASTQ R2 filename with 10x barcodes (gzipped).\n';
        printf '  - fastq_R3:   FASTQ R3 filename with tagmentation barcode,\n';
        printf '                spacer and ATAC read (gzipped).\n';
        printf '  - check|fix:  Check if given FASTQ files contain ScaleBio ATAC with 10x or\n';
        printf '                or HyDrop ATAC barcodes and optionally fix FASTQ files.\n'
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


    if type igzip > /dev/null 2>&1 ; then
        # Decompress gzipped FASTQ files with igzip if installed (6x faster than gzip).
        local decompress_fastq_gzipped_cmd="${decompress_fastq_igzip_cmd}";
    else
        # Decompress gzipped FASTQ files with gzip.
        local decompress_fastq_gzipped_cmd="${decompress_fastq_zcat_cmd}";
    fi


    # Detect if input FASTQ files are gzip compressed or not.
    if [ "${fastq_R2_filename}" = "${fastq_R2_filename%.gz}" ] ; then
        printf 'Error: FASTQ file "%s" is not gzipped.\n' "${fastq_R2_filename}";
        return 1;
    fi

    if [ "${fastq_R3_filename}" = "${fastq_R3_filename%.gz}" ] ; then
        printf 'Error: FASTQ file "%s" is not gzipped.\n' "${fastq_R3_filename}";
        return 1;
    fi

    if [ ! -e "${fastq_R2_filename}" ] ; then
        printf 'Error: FASTQ file "%s" does not exist.\n' "${fastq_R2_filename}";
        return 1;
    fi

    if [ ! -e "${fastq_R3_filename}" ] ; then
        printf 'Error: FASTQ file "%s" does not exist.\n' "${fastq_R3_filename}";
        return 1;
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

    # Check if input R3 (real read2) FASTQ file contains ScaleBio ATAC tagmentation barcodes.
    local is_scalebio_ATAC=$(
        mawk \
            -v fastq_R3_filename="${fastq_R3_filename}" \
            -v decompress_fastq_cmd="${decompress_fastq_gzipped_cmd}" \
        '
        BEGIN {
            read_fastq_R3_cmd = decompress_fastq_cmd " " fastq_R3_filename;

            fastq_line_number = 0;
            is_scalebio_ATAC = 0;

            # ScaleATAC tagmentation barcode whitelist (rev comp).
            scalebio_atac_tagmentation_barcodes["GAACCGCG"] = 0;
            scalebio_atac_tagmentation_barcodes["AGGTTATA"] = 0;
            scalebio_atac_tagmentation_barcodes["TCATCCTT"] = 0;
            scalebio_atac_tagmentation_barcodes["TGGCCGGT"] = 0;
            scalebio_atac_tagmentation_barcodes["CAATTAAC"] = 0;
            scalebio_atac_tagmentation_barcodes["ATAATGTG"] = 0;
            scalebio_atac_tagmentation_barcodes["TCTGTTGG"] = 0;
            scalebio_atac_tagmentation_barcodes["CTCACCAA"] = 0;
            scalebio_atac_tagmentation_barcodes["TATTAGCT"] = 0;
            scalebio_atac_tagmentation_barcodes["ATGTAAGT"] = 0;
            scalebio_atac_tagmentation_barcodes["GCACGGAC"] = 0;
            scalebio_atac_tagmentation_barcodes["GGTACCTT"] = 0;
            scalebio_atac_tagmentation_barcodes["ATCCACTG"] = 0;
            scalebio_atac_tagmentation_barcodes["GCTTGTCA"] = 0;
            scalebio_atac_tagmentation_barcodes["CAAGCTAG"] = 0;
            scalebio_atac_tagmentation_barcodes["TAAGTGGT"] = 0;
            scalebio_atac_tagmentation_barcodes["CGGACAAC"] = 0;
            scalebio_atac_tagmentation_barcodes["ATATGGAT"] = 0;
            scalebio_atac_tagmentation_barcodes["GCTCATTG"] = 0;
            scalebio_atac_tagmentation_barcodes["ATCTGCCA"] = 0;
            scalebio_atac_tagmentation_barcodes["CTTGGTAT"] = 0;
            scalebio_atac_tagmentation_barcodes["GATCTATC"] = 0;
            scalebio_atac_tagmentation_barcodes["AGCTCGCT"] = 0;
            scalebio_atac_tagmentation_barcodes["CGGAACTG"] = 0;
            # ScaleATAC tagmentation inhouse barcode whitelist (rev comp).
            scalebio_atac_tagmentation_barcodes["GGTCACGA"] = 0;
            scalebio_atac_tagmentation_barcodes["TCTCTACT"] = 0;
            scalebio_atac_tagmentation_barcodes["GCAGAATT"] = 0;
            scalebio_atac_tagmentation_barcodes["GCGGCACA"] = 0;
            scalebio_atac_tagmentation_barcodes["CTGCTTCC"] = 0;
            scalebio_atac_tagmentation_barcodes["AGTTCAGG"] = 0;
            scalebio_atac_tagmentation_barcodes["AACTGTAG"] = 0;
            scalebio_atac_tagmentation_barcodes["ATGAGGCC"] = 0;
            scalebio_atac_tagmentation_barcodes["GACCTGAA"] = 0;
            scalebio_atac_tagmentation_barcodes["TCGATATC"] = 0;
            scalebio_atac_tagmentation_barcodes["AAGATACT"] = 0;
            scalebio_atac_tagmentation_barcodes["CGCCGATC"] = 0;
            scalebio_atac_tagmentation_barcodes["CCATTCGA"] = 0;
            scalebio_atac_tagmentation_barcodes["CTCTCGTC"] = 0;
            scalebio_atac_tagmentation_barcodes["AACGTTCC"] = 0;
            scalebio_atac_tagmentation_barcodes["TTACAGGA"] = 0;
            scalebio_atac_tagmentation_barcodes["TGGATCGA"] = 0;
            scalebio_atac_tagmentation_barcodes["TCCAACGC"] = 0;
            scalebio_atac_tagmentation_barcodes["GCGCAAGC"] = 0;
            scalebio_atac_tagmentation_barcodes["CCGTGAAG"] = 0;
            scalebio_atac_tagmentation_barcodes["TAAGGTCA"] = 0;
            scalebio_atac_tagmentation_barcodes["TTGCCTAG"] = 0;
            scalebio_atac_tagmentation_barcodes["GGAGCGTC"] = 0;
            scalebio_atac_tagmentation_barcodes["CTAGCGCT"] = 0;

            hydrop_ligation1 = "GGGAC";
            hydrop_ligation2 = "GTCAG";

            # Read FASTQ input file.
            while ( (read_fastq_R3_cmd | getline) > 0 ) {
                fastq_line_number += 1;
                fastq_part = fastq_line_number % 4;

                if ( fastq_line_number == 100000 ) {
                    # Go to end block after reading 100000 lines.
                    exit;
                }

                if ( fastq_part == 1 ) {
                    if ( $2 !~ /^3:/ ) {
                        # FASTQ file is not R3 read.
                        exit;
                    }
                } else if ( fastq_part == 2 ) {
                    # Sequence lines.
                    seq_length = length($0);

                    if ( seq_length <= 8 + 19 + 30) {
                        # FASTQ file is likely not ScaleATAC as it should have:
                        #   - 8 bp tagmentation barcode
                        #   - 19 bp spacer
                        #   - ATAC sequence read (assume 30 bp is the minimum).
                        continue;
                    }

                    if ( substr($0, 1, 8)  in scalebio_atac_tagmentation_barcodes ) {
                        is_scalebio_ATAC += 1;
                    }
                }
            }
        }
        END {
            # Number of FASTQ reads read (25000 or less if small FASTQ file).
            nbr_fastq_reads = fastq_line_number / 4;

            # If 60% of the reads matched an exact contant sequence from ScaleATAC
            # tagmentation barcodes assume it contains ScaleATAC.
            if ( ( is_scalebio_ATAC / nbr_fastq_reads ) >= 0.6 ) {
                print "is_scalebio_ATAC";
            } else {
                printf("is_unknown: %f%% ScaleATAC\n", is_scalebio_ATAC / nbr_fastq_reads * 100.0);
            }
        }
        '
    );


    if [ "${is_scalebio_ATAC}" == "is_scalebio_ATAC" ] ; then
        # Check if input R2 (index 2) FASTQ file contain 10x or HyDrop ATAC barcodes.
        local is_hydrop_scalebio_ATAC=$(
            mawk \
                -v fastq_R2_filename="${fastq_R2_filename}" \
                -v decompress_fastq_cmd="${decompress_fastq_gzipped_cmd}" \
            '
            BEGIN {
                read_fastq_R2_cmd = decompress_fastq_cmd " " fastq_R2_filename;

                fastq_line_number = 0;
                is_hydrop_ATAC = 0;

                # ScaleATAC tagmentation barcode whitelist (rev comp).
                hydrop_ligation1 = "GGGAC";
                hydrop_ligation2 = "GTCAG";

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
                            # FASTQ file is not R3 read.
                            exit;
                        }
                    } else if ( fastq_part == 2 ) {
                        # Sequence lines.
                        seq_length = length($0);

                        if ( seq_length < 40) {
                            # FASTQ file is likely not HyDropATAC as it should have:
                            #   - 3x 10bp barcode separated with 2 5bp spacers.
                            continue;
                        }

                        if ( ( substr($0, 11, 5) == hydrop_ligation1 ) && ( substr($0, 26, 5) == hydrop_ligation2 ) ) {
                            is_hydrop_ATAC += 1;
                        }
                    }
                }
            }
            END {
                # Number of FASTQ reads read (25000 or less if small FASTQ file).
                nbr_fastq_reads = fastq_line_number / 4;

                # If 60% of the reads matched the 2 5bp hydrop ATAC spacers
                # assume it contains HyDrop ScaleATAC.
                if ( ( is_hydrop_ATAC / nbr_fastq_reads ) >= 0.6 ) {
                    print "is_hydrop_scalebio_ATAC";
                } else {
                    print "is_10x_scalebio_ATAC";
                }
            }
            '
        );

        if [ "${is_hydrop_scalebio_ATAC}" == "is_hydrop_scalebio_ATAC" ] ; then
            printf 'FASTQ "%s" contains ScaleBio ATAC tagmentation barcodes.\nFASTQ "%s" contains HyDrop ATAC barcodes.\n' "${fastq_R3_filename}" "${fastq_R2_filename}";
            is_hydrop_scalebio_ATAC=1;
        else
            printf 'FASTQ "%s" contains ScaleBio ATAC tagmentation barcodes.\nFASTQ "%s" contains 10x ATAC barcodes.\n' "${fastq_R3_filename}" "${fastq_R2_filename}";
            is_hydrop_scalebio_ATAC=0;
        fi
    else
        # Uncomment for debugging.
        #printf 'FASTQ "%s" contains %s\n' "${fastq_R3_filename}" "${is_scATAC_or_multiome_ATAC}";

        # Unknown: Do not modify the R2 and R3 read.
        return 0;
    fi

    if [ "${check_or_fix}" != "fix" ] ; then
        # Quit here if "fix" was not requested and just report if the FASTQ files
        # contained ScaleBio ATAC with 10x or HyDrop ATAC barcodes.
        return 0;
    fi

    local fastq_R2_output_filename="${fastq_R2_filename}.fixed.fq.gz";
    local fastq_R3_output_filename="${fastq_R3_filename}.fixed.fq.gz";


    mawk \
        -v fastq_R2_filename="${fastq_R2_filename}" \
        -v fastq_R3_filename="${fastq_R3_filename}" \
        -v fastq_R2_output_filename="${fastq_R2_output_filename}" \
        -v fastq_R3_output_filename="${fastq_R3_output_filename}" \
        -v decompress_fastq_cmd="${decompress_fastq_gzipped_cmd}" \
        -v compress_fastq_cmd="${compress_fastq_cmd}" \
        -v is_hydrop_scalebio_ATAC="${is_hydrop_scalebio_ATAC}" \
    '
    BEGIN {
        read_fastq_R2_cmd = decompress_fastq_cmd " " fastq_R2_filename;
        read_fastq_R3_cmd = decompress_fastq_cmd " " fastq_R3_filename;

        write_fastq_R2_cmd = compress_fastq_cmd " > " fastq_R2_output_filename;
        write_fastq_R3_cmd = compress_fastq_cmd " > " fastq_R3_output_filename;

        fastq_line_number = 0;

        # Read FASTQ R2 file (which contains the cell barcodes).
        while ( (read_fastq_R2_cmd | getline fastq_R2_line) > 0 ) {
            fastq_line_number += 1;
            fastq_part = fastq_line_number % 4;

            # Read FASTQ R3 file (which contains read 2).
            if ( (read_fastq_R3_cmd | getline fastq_R3_line) > 0 ) {
                if ( fastq_part == 1 ) {
                    # Read name lines.

                    # Find first space position (0 if no comment found) in read name from all input FASTQ files.
                    read_name_R2_space_pos = index(fastq_R2_line, " ");
                    read_name_R3_space_pos = index(fastq_R3_line, " ");

                    # Extract read name from all input FASTQ files.
                    if (read_name_R2_space_pos > 0) {
                        read_name_R2 = substr(fastq_R2_line, 2, read_name_R2_space_pos - 2);
                    } else {
                        read_name_R2 = substr(fastq_R2_line, 2);
                    }

                    if (read_name_R3_space_pos > 0) {
                        read_name_R3 = substr(fastq_R3_line, 2, read_name_R3_space_pos - 2);
                    } else {
                        read_name_R3 = substr(fastq_R3_line, 2);
                    }

                    # Check if read names match between all 2 FASTQ files.
                    if ( read_name_R2 != read_name_R3 ) {
                        print "Error: Read name R2 (\"" read_name_R2 "\") and R3 (\"" read_name_R3 "\") are not paired properly (line number: " fastq_line_number ").";
                        exit(1);
                    }

                    # Store full read names for R2 and R3 input FASTQ files.
                    read_name_R2_full = substr(fastq_R2_line, 2);
                    read_name_R3_full = substr(fastq_R3_line, 2);
                } else if ( fastq_part == 2 ) {
                    # Sequence lines.

                    # Store sequence info from R2 and R3 for later use.
                    sequence_R2 = fastq_R2_line;
                    sequence_R3 = fastq_R3_line;

                    # Check if R2 was not already corrected.
                    sequence_R2_length = length(sequence_R2)
                    if ( sequence_R2_length == 24 ) {
                        print "Error: Read name R2 (\"" read_name_R2 "\") has a sequence length of 24 instead of 16 and is probably already fixed (line number: " fastq_line_number ").";
                        exit(1);
                    } else if ( sequence_R2_length == 48 ) {
                        print "Error: Read name R2 (\"" read_name_R2 "\") has a sequence length of 48 instead of 40 and is probably already fixed (line number: " fastq_line_number ").";
                        exit(1);
                    } else if ( ( sequence_R2_length != 16 ) && ( sequence_R2_length != 40 ) ) {
                        print "Error: Read name R2 (\"" read_name_R2 "\") has a sequence length of " sequence_R2_length " instead of 16 (line number: " fastq_line_number ").";
                        exit(1);
                    }
                } else if ( fastq_part == 0 ) {
                    # Quality lines.

                    # Write the full FASTQ record to the R1 and R2 output FASTQ file with barcode info in front of the read name.
                    # When write_fastq_R1_cmd and write_fastq_R2_cmd are the same, an interleaved FASTQ fille will be written.
                    if ( is_hydrop_scalebio_ATAC == 1 ) {
                        # ScaleATAC with HyDrop ATAC barcodes.
                        print "@" read_name_R2_full "\n" substr(sequence_R2, 1, 40) substr(sequence_R3, 1, 8) "\n+\n" substr(fastq_R2_line, 1, 40) substr(fastq_R3_line, 1, 8) | write_fastq_R2_cmd;
                    } else {
                        # ScaleATAC with 10x ATAC barcodes.
                        print "@" read_name_R2_full "\n" substr(sequence_R2, 1, 16) substr(sequence_R3, 1, 8) "\n+\n" substr(fastq_R2_line, 1, 16) substr(fastq_R3_line, 1, 8) | write_fastq_R2_cmd;
                    }

                    print "@" read_name_R3_full "\n" substr(sequence_R3, 28) "\n+\n" substr(fastq_R3_line, 28) | write_fastq_R3_cmd;
                }
            }
        }

        # Close open file handles.
        close(read_fastq_R2_cmd);
        close(read_fastq_R3_cmd);
        close(write_fastq_R2_cmd);
        close(write_fastq_R3_cmd);
    }'

    if [ $? -eq 0 ] ; then
        # Rename fixed R2 and R3 FASTQ file to the original FASTQ R2 and R3 filename.
        mv "${fastq_R2_output_filename}" "${fastq_R2_filename}";
        mv "${fastq_R3_output_filename}" "${fastq_R3_filename}";
    fi


    return $?
}



fix_scalebio_atac_fastqs "${@}";

