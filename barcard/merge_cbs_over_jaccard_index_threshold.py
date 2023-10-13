#!/usr/bin/env python3

# Copyright (C) 2022-2023 - Gert Hulselmans


from __future__ import annotations

import argparse

import polars as pl


def merge_cbs_over_jaccard_index_threshold(
    CB1_vs_CB2_jaccard_filtered_tsv_filename: str,
    input_fragments_tsv_filename: str,
    output_fragments_tsv_filename: str,
) -> None:
    # Read CB1 vs CB2 Jaccard CSV file filtered by Otsu threshold.
    CB1_vs_CB2_jaccard_filtered_df = pl.read_csv(
        CB1_vs_CB2_jaccard_filtered_tsv_filename,
        separator="\t",
    )

    original_CB_to_merged_CBs_df = (
        # Get only "CB1" and "CB2" column and also append the same data,
        # but with the columns switched, so grouping on "CB1" will have
        # all possible CBs that potentially could get merged.
        pl.concat(
            [
                CB1_vs_CB2_jaccard_filtered_df.select(
                    [
                        pl.col("CB1"),
                        pl.col("CB2"),
                    ]
                ),
                CB1_vs_CB2_jaccard_filtered_df.select(
                    [
                        pl.col("CB2").alias("CB1"),
                        pl.col("CB1").alias("CB2"),
                    ]
                ),
            ]
        )
        .groupby("CB1")
        .agg(
            [
                # Get all CB2s that are linked to CB1.
                pl.col("CB2")
            ]
        )
        .select(
            [
                # Get CB1.
                pl.col("CB1").alias("original_CB"),
                # Get all CB2s liked to CB1, add CB1 to the list,
                # make unique and sort (so merged CB will be
                # independend of the order of the CBs) and convert
                # to string where merged CBs are separated by "_".
                pl.col("CB2")
                .list.concat([pl.col("CB1")])
                .list.unique()
                .list.sort()
                .list.join("_")
                .alias("merged_CBs"),
            ]
        )
    )

    # Write original CB to merged CBs TSV file.
    original_CB_to_merged_CBs_df.write_csv(
        f"{output_fragments_tsv_filename}.original_CB_to_merged_CBs.tsv",
        separator="\t",
    )

    fragments_df = (
        # Read original fragments file.
        pl.read_csv(
            input_fragments_tsv_filename,
            has_header=False,
            separator="\t",
            new_columns=["chrom", "start", "end", "CB", "CB_count"],
        )
        .join(
            # Add merged CBs as extra column or null if CB does not need to be merged.
            original_CB_to_merged_CBs_df,
            left_on="CB",
            right_on="original_CB",
            how="left",
        )
        .select(
            [
                pl.col(["chrom", "start", "end"]),
                # If the barcode needs to be merged, get "merged_CBs" column as CB,
                # else keep the original CB.
                pl.when(pl.col("merged_CBs").is_not_null())
                .then(pl.col("merged_CBs"))
                .otherwise(pl.col("CB"))
                .alias("CB"),
                pl.col("CB_count"),
            ]
        )
        .groupby(["chrom", "start", "end", "CB"], maintain_order=True)
        .agg(
            [
                # Combine counts for the same fragment for merged CBs.
                pl.col("CB_count")
                .sum()
                .alias("CB_count")
            ]
        )
        .write_csv(
            # Write output fragments TSV file with merged CBs.
            output_fragments_tsv_filename,
            separator="\t",
            has_header=False,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge CBs if they are over a precalculated Jaccard index threshold."
    )

    parser.add_argument(
        "-j",
        "--jaccard",
        dest="CB1_vs_CB2_jaccard_filtered_tsv_filename",
        action="store",
        type=str,
        required=True,
        help="Input TSV filename with CB1 vs CB2 with Jaccard index values filtered "
        "at a threshold. Barcodes above this threshold will be merged.",
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input_fragments_tsv_filename",
        action="store",
        type=str,
        required=True,
        help="Input fragments BED/TSV filename.",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_fragments_tsv_filename",
        action="store",
        type=str,
        required=True,
        help="Output fragments BED/TSV filename with CBs above Jaccard index merged.",
    )

    args = parser.parse_args()

    merge_cbs_over_jaccard_index_threshold(
        args.CB1_vs_CB2_jaccard_filtered_tsv_filename,
        args.input_fragments_tsv_filename,
        args.output_fragments_tsv_filename,
    )


if __name__ == "__main__":
    main()
