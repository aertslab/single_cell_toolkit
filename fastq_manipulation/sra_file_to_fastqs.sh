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

    if ! type bgzip > /dev/null 2>&1 ; then
         printf 'Error: "bgzip" (from HTSlib) is not installed.\n' >&2;
         return 1;
    fi

    if ! type mawk > /dev/null 2>&1 ; then
         printf 'Error: "mawk" is not installed.\n' >&2;
         return 1;
    fi

    if ! type fasterq-dump > /dev/null 2>&1 ; then
         printf 'Error: "fasterq-dump" (from SRA-Toolkit) is not installed.\n' >&2;
         return 1;
    fi

    if ! type vdb-dump > /dev/null 2>&1 ; then
         printf 'Error: "vdb-dump" (from SRA-Toolkit) is not installed.\n' >&2;
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

    local fastq='';
    local fastq_with_sample_name='';
    local -i read_number=1;

    if vdb-dump -R 1 -C LINKAGE_GROUP "${sra_file}" 2>&1 \
            | grep -q 'column not found while updating cursor within virtual database module' ; then
        true;
    else
        # Get script dir.
        local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

        # Old 10X sc/snRNA dataset in GEO sometimes have no read with the CB and UMI, but
        # store the CB and UMI in the "LINKAGE_GROUP" field of each read in the SRA file.
        printf 'Running "%s/extract_cb_umi_as_fastq_from_sra_file.sh %s %s" to extract the CB/UMI read as R1 read from "LINKAGE_GROUP" field.\n' \
            "${script_dir}" \
            "${sra_file}" \
            "${output_dir}";

        "${script_dir}/extract_cb_umi_as_fastq_from_sra_file.sh" "${sra_file}" "${output_dir}";

        # As the CB and UMI read is the R1 read, set the starting read number
        # for renaming FASTQ files produced by fasterq-dump to 2.
        read_number+=1;
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

    for fastq in "${output_dir%/}/${srr_number}_"[0-4]'.fastq' ; do
        fastq_with_sample_name="${output_dir%/}/${sample_name}___${srr_number}___S1_L001_R${read_number}_001.fastq";
        mv "${fastq}" "${fastq_with_sample_name}"
        echo "Compressing \"${fastq_with_sample_name}\"...";
        bgzip -@ 4 "${fastq_with_sample_name}";

        read_number+=1;
    done

    printf '\n';
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

    # Extract sample title or experiment title as fallback.
    local sample_name=$(
        esearch -db sra -query "${srr_number}" \
          | efetch -format xml \
          | xtract \
                -pattern EXPERIMENT_PACKAGE \
                -if SAMPLE/TITLE \
                  -element SAMPLE/TITLE \
                -else \
                  -element EXPERIMENT/TITLE
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
