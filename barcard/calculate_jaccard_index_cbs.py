#!/usr/bin/env python3

# Copyright (C) 2022 - Gert Hulselmans


import argparse

import polars as pl


def calculate_jaccard_index_cbs(
    fragments_tsv_filename: str,
    CB1_vs_CB2_jaccard_tsv_filename: str,
    min_frags_per_CB: int = 1000,
):
    # Read fragments file, count number of fragments per CB and
    # keep only those fragments which have a CB above or equal to min_frags_per_CB threshold.
    fragments_df = (
        pl.read_csv(
            fragments_tsv_filename,
            has_header=False,
            sep="\t",
            new_columns=["chrom", "start", "end", "CB", "CB_count"],
        )
        .with_columns(
            [
                pl.col("chrom").cast(pl.Categorical),
                pl.col("CB").cast(pl.Categorical),
            ]
        )
        .with_column(pl.col("CB").count().over("CB").alias("per_CB_count"))
        .filter(pl.col("per_CB_count") >= min_frags_per_CB)
    )

    CB1_vs_CB2_jaccard_df = (
        fragments_df.lazy()
        .join(
            # Get all combinations of CB1 with CB2 for the same fragments.
            fragments_df.lazy(),
            how="inner",
            on=["chrom", "start", "end"],
        )
        .select(
            # Rename columns after join.
            [
                pl.col("chrom"),
                pl.col("start"),
                pl.col("end"),
                pl.col("CB").alias("CB1"),
                pl.col("CB_count").alias("CB1_count"),
                pl.col("per_CB_count").alias("per_CB1_count"),
                pl.col("CB_right").alias("CB2"),
                pl.col("CB_count_right").alias("CB2_count"),
                pl.col("per_CB_count_right").alias("per_CB2_count"),
            ]
        )
        .filter(
            # Only keep fragment rows for which we have 2 different CBs.
            pl.col("CB1")
            != pl.col("CB2")
        )
        .groupby(["CB1", "CB2"])
        .agg(
            # Get for each CB1:
            #   - total number of fragments with CB1s.
            #   - total number of fragments with CB2s.
            #   - total number of fragments with CB1s and CBs2 (intersection).
            [
                pl.col("per_CB1_count").first(),
                pl.col("per_CB2_count").first(),
                pl.count().alias("CB1_vs_CB2_intersection"),
            ]
        )
        .with_columns(
            # Convert categorical CB columns to strings so "<" will sort based
            # on lexicographic order (instead of numerical categorical order).
            [
                pl.col("CB1").cast(pl.Utf8).alias("CB1_str"),
                pl.col("CB2").cast(pl.Utf8).alias("CB2_str"),
            ]
        )
        .with_column(
            # Make new "CB1_CB2" column with CB1 and CB2 sorted
            # so CB1 vs CB2 and CB2 vs CB1 give the same value,
            # so we can group on that value to not count the
            # same fragments overlap twice.
            pl.when(pl.col("CB1_str") < pl.col("CB2_str"))
            .then(pl.col("CB1_str") + "_" + pl.col("CB2_str"))
            .otherwise(pl.col("CB2_str") + "_" + pl.col("CB1_str"))
            .alias("CB1_CB2")
        )
        .sort(
            # Sort on CB1_CB2 and then on CB1 so the output is deterministic.
            by=["CB1_CB2", "CB1"],
        )
        .groupby(["CB1_CB2"], maintain_order=True)
        .agg(
            # Only keep the first CB1 vs CB2 row (to avoid counting
            # the same overlap twice).
            pl.col(
                [
                    "CB1",
                    "CB2",
                    "CB1_vs_CB2_intersection",
                    "per_CB1_count",
                    "per_CB2_count",
                ]
            ).first()
        )
        .drop(
            # Remove unneeded column.
            ["CB1_CB2"],
        )
        .with_column(
            # Calculate Jaccard index for CB1 vs CB2.
            (
                pl.col("CB1_vs_CB2_intersection")
                / (
                    pl.col("per_CB1_count")
                    + pl.col("per_CB2_count")
                    - pl.col("CB1_vs_CB2_intersection")
                )
            ).alias("jaccard")
        )
        .sort(
            # Sort on Jaccard index value.
            by=["jaccard"],
            reverse=True,
        )
        .with_column(
            # Include rank based on Jaccard index value.
            pl.col("jaccard")
            .rank(method="ordinal", reverse=True)
            .alias("jaccard_rank")
        )
        .collect()
    )

    # Write output to TSV file.
    CB1_vs_CB2_jaccard_df.write_csv(
        CB1_vs_CB2_jaccard_tsv_filename,
        sep="\t",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Calculate Jaccard index for CB pairs based on their common fragments."
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="fragments_tsv_filename",
        action="store",
        type=str,
        required=True,
        help="Input fragments BED/TSV filename.",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="CB1_vs_CB2_jaccard_tsv_filename",
        action="store",
        type=str,
        required=True,
        help="Output filename with CB1 vs CB2 with Jaccard index values.",
    )

    parser.add_argument(
        "-t",
        "--threshold",
        dest="min_frags_per_CB",
        action="store",
        type=int,
        required=False,
        default=1000,
        help="Minimum number of fragments needed per CB. Default: 1000.",
    )

    args = parser.parse_args()

    calculate_jaccard_index_cbs(
        args.fragments_tsv_filename,
        args.CB1_vs_CB2_jaccard_tsv_filename,
        min_frags_per_CB=args.min_frags_per_CB,
    )


if __name__ == "__main__":
    main()
