#!/bin/bash
#
# Copyright (C) 2024 - Gert Hulselmans
#
# Purpose:
#   Correct ScaleBio ATAC tagmentation barcode and use it to split read 1
#   and read 2 in different FASTQ files based on the tagmentation barcode.

set -e
set -o pipefail


SEQC_RUN="seqc run -release";
CODON_RUN="codon run -plugin seq -release";
CODON_OR_SEQ_RUN="${CODON_OR_SEQ_RUN:-${CODON_RUN}}";

cellranger_atac_737k_cratac_v1_barcode_list_filename='/staging/leuven/res_00001/barcodes/cellranger_atac.737K-cratac-v1.txt.gz';
hydrop_atac_v2_barcode_list_filename='/staging/leuven/res_00001/barcodes/HyDrop_v2.txt';



extract_and_correct_and_split_scalebio_atac_barcode_from_fastq () {
    local fastq_R1="${1}";
    local fastq_R2="${2}";
    local fastq_R3="${3}";
    local fastq_output_dir="${4}";
    local scalebio_bc_type="${5:-scalebioih}";
    local scatac_bc_type="${6:-10x}";

    if [ ${#@} -lt 4 ] ; then
        printf '\nUsage:\n';
        printf '    extract_and_correct_and_split_scalebio_atac_barcode_from_fastq \\\n';
        printf '        fastq_R1 \\\n';
        printf '        fastq_R2 \\\n';
        printf '        fastq_R3 \\\n';
        printf '        fastq_output_dir \\\n';
        printf '        [scalebio_bc_type] \\\n';
        printf '        [scatac_bc_type]\n\n';
        printf 'Purpose: Correct ScaleBio ATAC tagmentation barcode and use it to split read 1\n';
        printf '         and read 2 in different FASTQ files based on the tagmentation barcode.\n\n';
        printf 'Parameters:\n';
        printf '  - fastq_R1:   FASTQ R1 filename.\n';
        printf '  - fastq_R2:   FASTQ R2 filename.\n';
        printf '  - fastq_R3:   FASTQ R3 filename.\n';
        printf '  - fastq_output_dir: Output dir for FASTQ output file(s) (splitted per tagmentation barcode).\n';
        printf '  - scalebio_bc_type:\n';
        printf '        ScaleBio tagmentation barcode type set to look for:\n';
        printf '        "scalebio" (original ScaleBio) or "scalebioih" (ScaleBio inhouse: default).\n';
        printf '  - scatac_bc_type:\n';
        printf '        Single-cell ATAC barcode type set to look for:\n';
        printf '        "10x" (10x ATAC: default) or "hydrop" (HyDrop ATAC v2).\n';
        return 1;
    fi

    # Create correct tagmentation barcode filename.
    local corrected_bc_filename="${fastq_output_dir%/}/$(basename "${fastq_R2%.fastq.gz}.corrected_bc.zst")";

    scatac_bc_type="${scatac_bc_type,,}";

    if [ "${scatac_bc_type}" == "10x" ] ; then
        scatac_barcode_list_filename="${cellranger_atac_737k_cratac_v1_barcode_list_filename}";
    elif [ "${scatac_bc_type}" == "hydrop" ] ; then
        scatac_barcode_list_filename="${hydrop_atac_v2_barcode_list_filename}";
    else
        printf 'Error: Invalid scatac_bc_type "%s". Choose: "10x" or "hydrop".\n' "${scatac_bc_type}";
        return 1;
    fi


    if ! type "${CODON_OR_SEQ_RUN%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${CODON_OR_SEQ_RUN%% *}";
        return 1;
    fi

    if ! type zstd > /dev/null 2>&1 ; then
        printf 'Error: "zstd" not found or executable.\n';
        return 1;
    fi

    if ! type mawk > /dev/null 2>&1 ; then
        printf 'Error: "mawk" not found or executable.\n';
        return 1;
    fi

    if ! type igzip > /dev/null 2>&1 ; then
        printf 'Error: "igzip" not found or executable.\n';
        return 1;
    fi

    if ! type bgzip > /dev/null 2>&1 ; then
        printf 'Error: "bgzip" not found or executable.\n';
        return 1;
    fi


    # Get script dir.
    local script_dir="$(cd $(dirname "${BASH_SOURCE}") && pwd)";

    # Correct ScaleBio ATAC tagmentation barcodes.
    "${script_dir}/extract_and_correct_scalebio_atac_barcode_from_fastq.sh" \
        "${scatac_barcode_list_filename}" \
        "${fastq_R2}" \
        "${corrected_bc_filename}" \
        "${scalebio_bc_type}" \
        true \
        1 \
        1 \
        0.5;

    # Split original ScaleBio ATAC FASTQ files per tagmentation barcode.
    "${script_dir}/split_scalebio_atac_fastqs.sh" \
        "${fastq_R1}" \
        "${fastq_R3}" \
        "${corrected_bc_filename}" \
        "${fastq_output_dir}";

    return $?;
}



extract_and_correct_and_split_scalebio_atac_barcode_from_fastq "${@}";
