#!/bin/bash
#
# Copyright (C) 2021 - Gert Hulselmans
#
# Purpose:
#   Mark or remove duplicated reads per cell barcode.



mark_or_remove_duplicates_per_barcode () {
    local input_bam_file="${1}";
    local output_bam_file="${2}";
    local barcode_tag="${4}";
    local optical_duplicate_pixel_distance_arg="${5}";

    local remove_duplicates="false";
    local optical_duplicate_pixel_distance="100";

    if [ ${#@} -ne 5 ] ; then
        printf '\nUsage: mark_or_remove_duplicates_per_barcode input_bam output_bam mark_or_remove barcode_tag optical_duplicate_pixel_distance\n\n';
        printf 'Arguments:\n'
        printf '    - input_bam:       Input BAM file.\n';
        printf '    - output_bam:      Output BAM file.\n';
        printf '    - mark_or_remove:  "mark" or "remove" duplicated reads per cell barcode.\n';
        printf '    - barcode_tag:     Barcode tag to use ("CB", "CR", ...).\n';
        printf '    - optical_duplicate_pixel_distance:\n';
        printf '        The maximum offset between two duplicate clusters in order to\n';
        printf '        consider them optical duplicates.\n';
        printf '          - set actual optical duplicate pixel distance or use the\n';
        printf '            following keywords to set a reasonable default value.\n';
        printf '          - "unpatterned", "miseq", "nextseq500" or "hiseq2500": 100\n';
        printf '          - "patterned", "hiseq3000", "hiseq4000" or "novaseq6000": 2500\n\n';
        printf 'Purpose: Mark or remove duplicated reads per cell barcode.\n\n'

        return 1;
    fi

    case "${3}" in
        mark*|false)
            remove_duplicates="false";;
        remove*|true)
            remove_duplicates="true";;
        *)
            printf 'Error: Third parameter should be "mark" or "remove".\n' > /dev/stderr;
            return 1;;
    esac

    if [ "${#barcode_tag}" -ne 2 ] ; then
       printf 'Error: Barcode tag ("%s") should be 2 characters long.\n' "${barcode_tag}" > /dev/stderr;
       return 1;
    fi

    # Remove all underscores.
    optical_duplicate_pixel_distance_arg="${optical_duplicate_pixel_distance_arg//_/}";

    # Convert to lowercase.
    optical_duplicate_pixel_distance_arg="${optical_duplicate_pixel_distance_arg,,}";

    case "${optical_duplicate_pixel_distance_arg}" in
        [0-9]*)
            optical_duplicate_pixel_distance="${optical_duplicate_pixel_distance_arg}";;
        unpatterned|miseq*|nextseq*|hiseq2500)
            optical_duplicate_pixel_distance="100";;
        patterned|hiseq3000|hiseq4000|novaseq6000)
            optical_duplicate_pixel_distance="2500";;
        *)
            printf 'Error: Invalue optical duplicate pixel distance parameter.\n' > /dev/stderr;
            return 1;;
    esac

    module load 'picard/2.23.4-Java-1.8.0';

    java -Dpicard.useLegacyParser=false -jar "${EBROOTPICARD}/picard.jar" MarkDuplicates \
        -I "${input_bam_file}" \
        -O "${output_bam_file}" \
        -M "${output_bam_file%.bam}.mark_duplicates_metrics.txt" \
        -BARCODE_TAG "${barcode_tag}" \
        -REMOVE_DUPLICATES "${remove_duplicates}" \
        -OPTICAL_DUPLICATE_PIXEL_DISTANCE "${optical_duplicate_pixel_distance}";

    return $?;
}



mark_or_remove_duplicates_per_barcode "${@}";
