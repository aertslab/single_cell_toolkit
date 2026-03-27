#!/bin/bash
#
# Copyright (C) 2025 - Gert Hulselmans
#
# Purpose:
#   Extract CB and UMI (10X sc/snRNA) from "LINKAGE_GROUP" field of SRA file and convert to FASTQ:
#     - Some old 10X sc/snRNA dataset on GEO do not have the read with
#       the CB and UMI, and thus can not be dumped with fasterq-dump.
#       For some of those datasets the CB and UMI is stored in the SRA
#       file in the "LINKAGE_GROUP" field of each read.
#     - As only CB and UB tags are available in the "LINKAGE_GROUP" field
#       original quality scores are not available and quality scores are
#       set to "E" (phred score 36) in the FASTQ file.
#     - If no CB is found for a read, create a fake CB and UMI with all Gs.
#       and use minimal quality comment line.
#     - Add human readable sample name to output FASTQ filenames.
#     - Add suffix to FASTQ filenames so CellRanger can work with those
#       FASTQ files.



extract_cb_umi_as_fastq_from_sra_file () {
    local sra_file="${1}";
    local output_dir="${2}";

    if [ ${#@} -ne 2 ] ; then
        printf '\nUsage:\n';
        printf '    extract_cb_umi_as_fastq_from_sra_file sra_file output_dir\n\n';
        printf 'Purpose:\n';
        printf '    Extract CB and UMI (10X sc/snRNA) from "LINKAGE_GROUP" field of SRA file and convert to FASTQ:\n';
        printf '      - Some old 10X sc/snRNA dataset on GEO do not have the read with\n';
        printf '        the CB and UMI, and thus can not be dumped with fasterq-dump.\n';
        printf '        For some of those datasets the CB and UMI is stored in the SRA\n';
        printf '        file in the "LINKAGE_GROUP" field of each read.\n';
        printf '      - As only CB and UB tags are available in the "LINKAGE_GROUP" field\n';
        printf '        original quality scores are not available and quality scores are\n';
        printf '        set to "E" (phred score 36) in the FASTQ file.\n';
        printf '      - If no CB is found for a read, create a fake CB and UMI with all Gs.\n';
        printf '        and use minimal quality comment line.\n';
        printf '      - Add human readable sample name to output FASTQ filenames.\n';
        printf '      - Add suffix to FASTQ filenames so CellRanger can work with those\n';
        printf '        FASTQ files.\n\n';
        printf 'Parameters:\n';
        printf '    sra_file:    SRA input filename.\n'
        printf '                 e.g.: downloaded with "prefetch -O ./ --max-size 100G ${SRR_number}"\n';
        printf '    output_dir:  Output dir to which the final FASTQ files will be written.\n\n';
        return 1;
    fi

    if ! type bgzip > /dev/null 2>&1 ; then
         printf 'Error: "bgzip" (from HTSlib) is not installed.\n' >&2;
         return 1;
    fi

    if ! type mawk > /dev/null 2>&1 ; then
         printf 'Error: "mawk" is not installed.\n' >&2;
         return 1;
    fi

    if ! type vdb-dump > /dev/null 2>&1 ; then
         printf 'Error: "vdb-dump" (from SRA-Toolkit) is not installed.\n' >&2;
         return 1;
    fi

    if vdb-dump -R 1 -C LINKAGE_GROUP "${sra_file}" 2>&1 \
            | grep -q 'column not found while updating cursor within virtual database module' ; then
        printf 'Error: SRA file "%s" does not contain a "LINKAGE_GROUP" field.\n' "${sra_file}";
        return 1;
    fi

    # Get SRR number from SRA filename:
    #   - Remove ".sra" extension.
    #   - Remove parent dirs.
    local srr_number="${sra_file%.sra}";
    srr_number="${srr_number##*/}";

    echo "Get sample name for \"${srr_number}\"..."

    # Get sample name for SRR number.
    local sample_name=$(srr_number_to_sample_name "${srr_number}");

    if [ $? -ne 0 ] ; then
        # Exit if EDirect tools are not found.
        return 1;
    fi

    echo "Sample name for \"${srr_number}\": \"${sample_name}\".";

    # Create output FASTQ name.
    fastq_with_sample_name="${output_dir%/}/${sample_name}___${srr_number}___S1_L001_R1_001.fastq.gz";

    # Extract sequence name and CB and UMI and convert to bgzipped FASTQ format
    # (with fake quality scores).
    vdb-dump \
        --format tab \
        -C 'NAME,LINKAGE_GROUP' \
        "${sra_file}" \
      | mawk -F '\t|:|-' '
            {
                read_name = $1;
                # Cell barcode without "-1".
                cb = $3;
                umi = $5;

                if (cb == "") {
                    # No cell barcode defined for current read, set a fake one.
                    cb = "GGGGGGGGGGGGGGGG";
                    umi = "GGGGGGGGGG";
                }

                print "@" read_name "\n" cb umi "\n+\nEEEEEEEEEEEEEEEEEEEEEEEEEE" }' \
      | bgzip -@ 4 \
      > "${fastq_with_sample_name}";
}



srr_number_to_sample_name () {
    local srr_number="${1}";

    if [ ${#@} -ne 1 ] ; then
        printf '\nUsage:\n';
        printf '    srr_number_to_sample_name srr_number\n\n';
        printf 'Purpose:\n';
        printf '    Lookup human readable sample name that belongs to a certain SRR number.\n\n';
        printf 'Parameters:\n';
        printf '    srr_number: SRR number for which to retrieve the sample name.\n\n'
        return 1;
    fi

    if ! type efetch > /dev/null 2>&1 ; then
         printf 'Error: "efetch" (from EDirect) is not installed.\n' >&2;
         return 1;
    fi

    if ! type esearch > /dev/null 2>&1 ; then
         printf 'Error: "esearch" (from EDirect) is not installed.\n' >&2;
         return 1;
    fi

    if ! type xtract > /dev/null 2>&1 ; then
         printf 'Error: "xtract" (from EDirect) is not installed.\n' >&2;
         return 1;
    fi

    local sample_name=$(
        esearch -db sra -query "${srr_number}" \
          | efetch -format xml \
          | xtract -pattern SAMPLE -element TITLE
    );

    # Replace some "dangerous" characters with "_".
    sample_name="${sample_name// /_}";
    sample_name="${sample_name//,/_}";
    sample_name="${sample_name//\?/_}";
    sample_name="${sample_name//\//_}";
    sample_name="${sample_name//\\/_}";
    sample_name="${sample_name//[/_}";
    sample_name="${sample_name//]/_}";
    sample_name="${sample_name//{/_}";
    sample_name="${sample_name//\}/_}";

    echo "${sample_name}";
}



extract_cb_umi_as_fastq_from_sra_file "${@}";
