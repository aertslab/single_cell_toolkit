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
    #       - Filter BAM file similar as create_fragments_file but with SAMtools new expression syntax.
    #           samtools view -c --expr 'flag.proper_pair && !flag.secondary && !flag.supplementary && refid == mrefid && (tlen >= 10 || tlen <= -10) && mapq >= 30 && [MQ] >= 30 && [CB]' "${bam}"
    #
    #           Count read if and only if:
    #             - read is properly paired.
    #             - read is primary alignment.
    #             - read and its pair are located on the same chromosome.
    #             - read and its pair have a mapping quality of 30 or higher.
    #             - insert size is at least 10 in absolute value.
    #             - mapping quality is >=30 and mate mapping quality is >=30.
    #             - read has an associated CB tag.
    #   - Only use threads for "samtools stat". Using it with any of the other samtools commands
    #     makes everything slower than not using any threads at all.

    samtools view -u "${bam}" \
      | tee \
            >(samtools view -c -F 0x4 -F 0x100 -F 0x800 -e '! [XA] && ! [SA]' - > "${output_dir}/${sample_id}.uniquely_mapped_reads.txt") \
            >(samtools view -c -F 0x4 -F 0x100 -F 0x800 -q 30 - > "${output_dir}/${sample_id}.fraction_total_read_pairs.txt") \
            >(samtools view -c --expr 'flag.proper_pair && !flag.secondary && !flag.supplementary && refid == mrefid && (tlen >= 10 || tlen <= -10) && mapq >= 30 && [MQ] >= 30 && [CB]' - > "${output_dir}/${sample_id}.reads_used_for_fragments_file.txt") \
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

        printf "reads used to create fragments file\t";
        cat "${output_dir}/${sample_id}.reads_used_for_fragments_file.txt";
    ) > "${output_dir}/${sample_id}.mapping_stats.tsv";

    rm \
        "${output_dir}/${sample_id}.stat" \
        "${output_dir}/${sample_id}.uniquely_mapped_reads.txt" \
        "${output_dir}/${sample_id}.fraction_total_read_pairs.txt" \
        "${output_dir}/${sample_id}.reads_used_for_fragments_file.txt";
}



mapping_stats "${@}";
