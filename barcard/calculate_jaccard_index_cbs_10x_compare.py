#!/usr/bin/env python3

# Copyright (C) 2022-2023 - Gert Hulselmans


from __future__ import annotations

import argparse
from typing import Sequence

import polars as pl


def get_mo_rna_atac_cb_mapping() -> pl.DataFrame:
    cellranger_mo_rna_cb_filename = (
        "/staging/leuven/res_00001/barcodes/cellranger_arc_rna.737K-arc-v1.txt"
    )
    cellranger_mo_atac_cb_filename = (
        "/staging/leuven/res_00001/barcodes/cellranger_arc_atac.737K-arc-v1.txt"
    )
    cellranger_mo_atac_rev_comp_cb_filename = "/staging/leuven/res_00001/barcodes/cellranger_arc_atac.737K-arc-v1.REV_COMP.txt"

    mo_rna_atac_cb_df = pl.concat(
        [
            pl.read_csv(
                cellranger_mo_rna_cb_filename,
                separator="\t",
                has_header=False,
                new_columns=["CB_mo_rna"],
            ),
            pl.read_csv(
                cellranger_mo_atac_cb_filename,
                separator="\t",
                has_header=False,
                new_columns=["CB_mo_atac"],
            ),
            pl.read_csv(
                cellranger_mo_atac_rev_comp_cb_filename,
                separator="\t",
                has_header=False,
                new_columns=["CB_mo_atac_rev_comp"],
            ),
        ],
        how="horizontal",
    )
    return mo_rna_atac_cb_df


def get_atac_cb_mapping() -> pl.DataFrame:
    cellranger_atac_cb_filename = (
        "/staging/leuven/res_00001/barcodes/cellranger_atac.737K-cratac-v1.txt"
    )
    cellranger_atac_rev_comp_cb_filename = (
        "/staging/leuven/res_00001/barcodes/cellranger_atac.737K-cratac-v1.REV_COMP.txt"
    )

    cellranger_atac_cb_df = pl.concat(
        [
            pl.read_csv(
                cellranger_atac_cb_filename,
                separator="\t",
                has_header=False,
                new_columns=["CB_mo_rna"],
            ),
            pl.read_csv(
                cellranger_atac_cb_filename,
                separator="\t",
                has_header=False,
                new_columns=["CB_mo_atac"],
            ),
            pl.read_csv(
                cellranger_atac_rev_comp_cb_filename,
                separator="\t",
                has_header=False,
                new_columns=["CB_mo_atac_rev_comp"],
            ),
        ],
        how="horizontal",
    )
    return cellranger_atac_cb_df


def read_fragments_file(
    fragments_tsv_filename: str,
    pumatac_barcodes_tsv_filename: str,
    mo_rna_atac_cb_df: pl.DataFrame | None = None,
    chromosomes: str | Sequence | None = None,
) -> pl.DataFrame:
    print(f'Reading fragments file "{fragments_tsv_filename}" ...')

    vsc_atac_barcodes_df = pl.read_csv(
        pumatac_barcodes_tsv_filename,
        has_header=False,
        separator="\t",
        comment_prefix="#",
        new_columns=["CB"],
    )

    if isinstance(mo_rna_atac_cb_df, pl.DataFrame):
        cellranger_arc_rev_comp_barcodes_df = vsc_atac_barcodes_df.join(
            mo_rna_atac_cb_df,
            left_on="CB",
            right_on="CB_mo_atac_rev_comp",
            how="inner",
        )

        cellranger_arc_barcodes_df = vsc_atac_barcodes_df.join(
            mo_rna_atac_cb_df,
            left_on="CB",
            right_on="CB_mo_atac",
            how="inner",
        )

        if (
            cellranger_arc_rev_comp_barcodes_df.height
            > cellranger_arc_barcodes_df.height
        ):
            selected_barcodes = cellranger_arc_rev_comp_barcodes_df.select(
                [pl.col("CB_mo_rna").alias("CB_select"), pl.col("CB").alias("CB_final")]
            )
        else:
            selected_barcodes = cellranger_arc_barcodes_df.select(
                [pl.col("CB_mo_rna").alias("CB_select"), pl.col("CB").alias("CB_final")]
            )
    else:
        selected_barcodes = vsc_atac_barcodes_df.select(
            [pl.col("CB").alias("CB_select"), pl.col("CB").alias("CB_final")]
        )

    # Read fragments file, count number of fragments per CB and
    # keep only those fragments which have a CB above or equal to min_frags_per_CB threshold.
    fragments_df = (
        pl.read_csv(
            fragments_tsv_filename,
            has_header=False,
            separator="\t",
            comment_prefix="#",
            new_columns=["chrom", "start", "end", "CB", "CB_count"],
        )
        .with_columns(
            [
                pl.col("CB").str.replace("-[0-9]$", ""),
            ]
        )
        .join(
            selected_barcodes,
            left_on="CB",
            right_on="CB_select",
            how="inner",
        )
        .select(
            [
                pl.col(["chrom", "start", "end"]),
                pl.col("CB_final").alias("CB"),
                pl.col("CB_count"),
            ]
        )
    )

    if chromosomes:
        if isinstance(chromosomes, str):
            if chromosomes.find(",") != -1:
                # Chromosomes names to keep were specified as a comma separated list.
                chromosomes = chromosomes.split(",")
            elif chromosomes.find(" ") != -1:
                # Chromosomes names to keep were specified as a space separated list.
                chromosomes = chromosomes.split(" ")

        if isinstance(chromosomes, str):
            # Assume chromosome names to keep is a regex.
            chroms_to_keep = fragments_df.select(
                # Get all chromosome names available in the fragments file.
                pl.col("chrom").unique(),
            ).filter(
                pl.col("chrom").cast(pl.Utf8).str.contains(chromosomes),
            )
        elif isinstance(chromosomes, Sequence):
            chroms_to_keep = fragments_df.select(
                # Get all chromosome names available in the fragments file.
                pl.col("chrom").unique(),
            ).filter(
                # Only keep chromosomes available in the list.
                pl.col("chrom")
                .cast(pl.Utf8)
                .is_in(chromosomes),
            )

        print(
            "Keeping chromosomes: "
            + ", ".join(
                [
                    str(chr)
                    for chr in chroms_to_keep.get_column("chrom").sort().to_list()
                ]
            )
        )

        fragments_df = fragments_df.join(chroms_to_keep, how="inner", on="chrom")

    return fragments_df


def calculate_jaccard_index_cbs(
    pumatac_fragments_df: pl.DataFrame,
    cellranger_arc_fragments_df: pl.DataFrame,
    CB1_vs_CB2_jaccard_tsv_filename: str,
    chromosomes: str | Sequence | None = None,
) -> None:
    pumatac_fragments_df = pumatac_fragments_df.with_columns(
        pl.col("CB").count().over("CB").alias("per_CB_count")
    )
    cellranger_arc_fragments_df = cellranger_arc_fragments_df.with_columns(
        pl.col("CB").count().over("CB").alias("per_CB_count")
    )

    vsn_cellranger_outer_join_fragments_df = pumatac_fragments_df.join(
        cellranger_arc_fragments_df,
        on=["chrom", "start", "end", "CB"],
        how="outer",
    )

    print("Calculate Jaccard index for CB pairs based on their common fragments ...")

    # fe
    CB1_vs_CB2_jaccard_df = (
        # Get all fragments from VSN atac fragments file and CellRanger-ARC
        # fragments file and annotate if the fragment was found in them.
        pumatac_fragments_df.join(
            cellranger_arc_fragments_df,
            on=["chrom", "start", "end", "CB"],
            how="outer",
            suffix="2",
        )
        .group_by(["CB"])
        .agg(
            [
                #   - total number of fragments with CB1s and CBs2 (intersection).
                #   - total number of fragments with CB1s.
                #   - total number of fragments with CB2s.
                (
                    pl.col("per_CB_count").is_not_null()
                    & pl.col("per_CB_count2").is_not_null()
                )
                .sum()
                .alias("CB1_vs_CB2_intersection"),
                pl.col("per_CB_count").max().alias("per_CB1_count"),
                pl.col("per_CB_count2").max().alias("per_CB2_count"),
            ]
        )
        .with_columns(
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
            descending=True,
        )
        .with_columns(
            # Include rank based on Jaccard index value.
            pl.col("jaccard")
            .rank(method="ordinal", descending=True)
            .alias("jaccard_rank")
        )
    )

    print(
        "Write Jaccard index for CB pairs based on their common fragments to "
        f'"{CB1_vs_CB2_jaccard_tsv_filename}" ...'
    )

    # Write output to TSV file.
    CB1_vs_CB2_jaccard_df.write_csv(
        CB1_vs_CB2_jaccard_tsv_filename,
        separator="\t",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate Jaccard index between PUMATAC and CellRanger-ARC fragment files for selected CBs based on their common fragments."
    )

    parser.add_argument(
        "--vsn",
        dest="pumatac_fragments_tsv_filename",
        action="store",
        type=str,
        required=True,
        help="Input PUMATAC fragments BED/TSV filename.",
    )

    parser.add_argument(
        "--cellranger-atac",
        dest="cellranger_atac_fragments_tsv_filename",
        action="store",
        type=str,
        required=False,
        help="Input CellRanger-ATAC fragments BED/TSV filename.",
    )

    parser.add_argument(
        "--cellranger-arc",
        dest="cellranger_arc_fragments_tsv_filename",
        action="store",
        type=str,
        required=False,
        help="Input CellRanger-ARC fragments BED/TSV filename.",
    )

    parser.add_argument(
        "--barcodes",
        dest="pumatac_barcodes_tsv_filename",
        action="store",
        type=str,
        required=True,
        help="Input PUMATAC cell barcodes filename.",
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
        "-c",
        "--chromosomes",
        dest="chromosomes",
        action="store",
        type=str,
        required=False,
        default=None,
        help="Only use specified chromosome names. "
        "Specify chromosomes to keep as a comma or space separated list "
        '(e.g.: "chr1,chr2,chr3" or "chr1 chr2 chr3") or a regular expression '
        '(e.g.: "^(chr)?([0-9]+|[XY])$"), or None to keep all chromosomes. '
        "Default: None.",
    )

    args = parser.parse_args()

    if args.cellranger_arc_fragments_tsv_filename:
        mo_rna_atac_cb_df = get_mo_rna_atac_cb_mapping()

    if args.cellranger_atac_fragments_tsv_filename:
        mo_rna_atac_cb_df = get_atac_cb_mapping()

    pumatac_fragments_df = read_fragments_file(
        fragments_tsv_filename=args.pumatac_fragments_tsv_filename,
        pumatac_barcodes_tsv_filename=args.pumatac_barcodes_tsv_filename,
        mo_rna_atac_cb_df=None,
        chromosomes=args.chromosomes,
    )

    cellranger_arc_fragments_df = read_fragments_file(
        fragments_tsv_filename=args.cellranger_arc_fragments_tsv_filename
        if args.cellranger_arc_fragments_tsv_filename
        else args.cellranger_atac_fragments_tsv_filename,
        pumatac_barcodes_tsv_filename=args.pumatac_barcodes_tsv_filename,
        mo_rna_atac_cb_df=mo_rna_atac_cb_df,
        chromosomes=args.chromosomes,
    )

    calculate_jaccard_index_cbs(
        pumatac_fragments_df=pumatac_fragments_df,
        cellranger_arc_fragments_df=cellranger_arc_fragments_df,
        CB1_vs_CB2_jaccard_tsv_filename=args.CB1_vs_CB2_jaccard_tsv_filename,
    )


if __name__ == "__main__":
    main()
