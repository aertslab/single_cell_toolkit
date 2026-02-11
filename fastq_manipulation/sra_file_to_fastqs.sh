#!/bin/bash
#
# Copyright (C) 2025 - Gert Hulselmans
#
# Purpose:
#   Convert SRA file to multiple FASTQ files:
#     - Include all useful reads:
#       - biological and technical reads (e.g. scATAC CB read).
#       - exclude index reads which are not needed for mapping.
#     - Use original read names (Illumnia read name) in FASTQ
#       file and use minimal quality comment line.
#     - Add human readable sample name to output FASTQ filenames.
#     - Add suffix to FASTQ filenames so CellRanger(/ATAC/ARC)
#       can work with those FASTQ files.



# Define temp dir location for fasterq-dump.
# To set a different tmp dir than "/tmp", set ${TMPDIR} environment variable
# with export before running this script:
#
#   export TMPDIR="/your/tmpdir";
FASTERQ_TMP_DIR="${TMPDIR:-/tmp}";



sra_file_to_fastqs () {
    # SRA file with SRR number.
    local sra_file="${1}";
    local output_dir="${2}";

    if [ ${#@} -ne 2 ] ; then
        printf '\nUsage:\n';
        printf '    sra_file_to_fastqs sra_file output_dir\n\n';
        printf 'Purpose:\n';
        printf '    Convert SRA file to multiple FASTQ files:\n';
        printf '      - Include all useful reads:\n';
        printf '        - biological and technical reads (e.g. scATAC CB read).\n';
        printf '        - exclude index reads which are not needed for mapping.\n';
        printf '      - Use original read names (Illumnia read name) in FASTQ file\n';
        printf '        and use minimal quality comment line.\n';
        printf '      - Add human readable sample name to output FASTQ filenames.\n';
        printf '      - Add suffix to FASTQ filenames so CellRanger(/ATAC/ARC)\n';
        printf '        can work with those FASTQ files.\n\n';
        printf 'Parameters:\n';
        printf '    sra_file:    SRA input filename.\n'
        printf '                 e.g.: downloaded with "prefetch -O ./ --max-size 100G ${SRR_number}"\n';
        printf '    output_dir:  Output dir to which the final FASTQ files will be written.\n\n';
        return 1;
    fi

    echo "Convert SRA file \"${sra_file}\" to multiple FASTQ files...";

    # Convert SRA file to multiple FASTQ files:
    #   - Write only Illumina read name for line 1 of a FASTQ record (seq-defline)
    #     instead of "@SRR_seq_id original_read_name length=X".
    #   - Write only "+" for line 3 of a FASTQ record (qual_defline)
    #     instead of "+SRR_seq_id original_read_name length=X".
    #   - Write all (biological and technical) reads (e.g. R1, R2, R3, R4) to
    #     different files so also reads with scATAC cell barcodes are written,
    #     but avoid writing pure index reads that are never needed (10 bp or less).
    fasterq-dump \
        --use-name \
        --seq-defline '@$sn' \
        --qual-defline '+' \
        --min-read-len 11 \
        --include-technical \
        --split-files \
        --outdir "${output_dir}" \
        -t "${FASTERQ_TMP_DIR}" \
        "${sra_file}";

    # Get SRR number from SRA filename:
    #   - Remove ".sra" extension.
    #   - Remove parent dirs.
    local srr_number="${sra_file%.sra}";
    srr_number="${srr_number##*/}";

    echo "Get sample name for \"${srr_number}\"..."

    # Get sample name for SRR number.
    local sample_name=$(srr_number_to_sample_name "${srr_number}");

    echo "Sample name for \"${srr_number}\": \"${sample_name}\".";

    local fastq='';
    local fastq_with_sample_name='';
    local -i read_number=1;

    for fastq in "${output_dir%/}/${srr_number}_"[0-4]'.fastq' ; do
        fastq_with_sample_name="${output_dir%/}/${sample_name}___${srr_number}___S1_R${read_number}_001.fastq";
        mv "${fastq}" "${fastq_with_sample_name}"
        echo "Compressing \"${fastq_with_sample_name}\"...";
        bgzip -@ 4 "${fastq_with_sample_name}";

        read_number+=1;
    done
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



sra_file_to_fastqs "${@}";
