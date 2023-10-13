#!/usr/bin/env bash



mapping_stats () {
    local sample_id="${1}";
    local bam="${2}";
    local output_dir="${3}";

    if [ ${#@} -ne 3 ] ; then
        printf 'Usage:   mapping_stats.sh sample_id bam_file output_dir\n' >&2;
        printf 'Purpose: Get mapping statistics from BAM file\n' >&2;
        return 1;
    fi

    # Get mapping statistics from BAM file:
    #   - Read BAM file and write uncompressed BAM.
    #   - Uncompressed BAM file is written to each samtools command with tee (writes to each specified file and stdout).
    #   - samtools commands:
    #       - Get samtools statistics with:
    #           samtools stat "${bam}" > "${sample_id}.stat"
    #       - Uniquely mapped reads (BWA):
    #           samtools view -c -F 0x4 -F 0x100 -F 0x800 -e '! [XA] && ! [SA]' "${bam}"
    #       - Fraction of total read pairs mapped confidently to genome (>30 mapq):
    #           samtools view -c -F 0x4 -F 0x100 -F 0x800 -q 30 "${bam}"
    #   - Only use threads for "samtools stat". Using it with any of the other samtools commands
    #     makes everything slower than not using any threads at all.
    samtools view -u "${bam}" \
      | tee \
            >(samtools view -c -F 0x4 -F 0x100 -F 0x800 -e '! [XA] && ! [SA]' - > "${output_dir}/${sample_id}.uniquely_mapped_reads.txt") \
            >(samtools view -c -F 0x4 -F 0x100 -F 0x800 -q 30 - > "${output_dir}/${sample_id}.fraction_total_read_pairs.txt") \
      | samtools stat -@ 2 - > "${output_dir}/${sample_id}.stat"

    # Create one output file with combined mapping statistics from 3 samtools invocations:
    (
        printf "sample_id\t${sample_id}\n";

        awk \
            -F ':?\t' \
            '
            {
                if ( $1 == "SN" ) {
                   print $2 "\t" $3;
                }
            }
            ' \
            "${output_dir}/${sample_id}.stat";

        printf "uniquely mapped reads\t";
        cat "${output_dir}/${sample_id}.uniquely_mapped_reads.txt";

        printf "reads mapped with MAPQ>30\t";
        cat "${output_dir}/${sample_id}.fraction_total_read_pairs.txt";
    ) > "${output_dir}/${sample_id}.mapping_stats.tsv";

    rm \
        "${output_dir}/${sample_id}.stat" \
        "${output_dir}/${sample_id}.uniquely_mapped_reads.txt" \
        "${output_dir}/${sample_id}.fraction_total_read_pairs.txt";
}



mapping_stats "${@}";
