#!/usr/bin/env python3

import argparse
import os
import shutil
import sys

import xopen

import polars as pl


def get_orig_cb_idx_to_cb(barcodes_tsv_filename):
    """Create a Polars DataFrame with original CB indices and CB names from barcodes TSV file."""
    cb_idx_orig_to_cb_df = pl.read_csv(
        barcodes_tsv_filename,
        sep="\t",
        has_header=False,
        comment_char="#",
        columns=["column_1"],
        new_columns=["CB"],
    ).with_row_count("CB_idx_orig", offset=1)

    return cb_idx_orig_to_cb_df


def read_mtx(matrix_mtx_filename):
    """Read MatrixMarket matrix file as Polars DataFrame and a header (first 3 lines)."""
    matrix_mtx_header = b""

    with xopen.xopen(matrix_mtx_filename, "rb") as fh:
        for idx, line in enumerate(fh):
            if idx == 0:
                if line != b"%%MatrixMarket matrix coordinate integer general\n":
                    print(line)
                    raise ValueError(
                        f'"{matrix_mtx_filename}" is not in "MatrixMarket matrix coordinate integer general" format.'
                    )
            if idx <= 2:
                matrix_mtx_header += line
            else:
                break

    matrix_mtx_df = pl.read_csv(
        matrix_mtx_filename,
        sep=" ",
        has_header=False,
        skip_rows=3,
        columns=["column_1", "column_2", "column_3"],
        new_columns=["feature_idx", "CB_idx_orig", "value"],
        dtypes={"feature_idx": pl.UInt32, "CB_idx_orig": pl.UInt32, "value": pl.UInt32},
    )

    return matrix_mtx_df, matrix_mtx_header


def write_filtered_barcodes_and_matrix_mtx(
    cb_idx_orig_to_cb_df,
    matrix_mtx_df,
    matrix_mtx_header,
    barcodes_tsv_out_filename,
    mtx_out_filename,
):
    # Get CB idx mapping from original CB indices to new CB indices for CBs that appear
    # in the MatrixMarket matrix dataframe at least once.
    filtered_cb_idx_orig_to_cb_idx_new_df = (
        # Get all unique CB index values from the MatrixMarket matrix dataframe.
        matrix_mtx_df.lazy()
        .select(pl.col("CB_idx_orig").unique())
        # And keep only those CB names from the barcode dataframe, which appear in the
        # MatrixMarket matrix dataframe at least once.
        .join(
            cb_idx_orig_to_cb_df.lazy(),
            how="left",
            on="CB_idx_orig",
        )
        # And create a new CB idx value for those selected CBs.
        .with_row_count("CB_idx_new", offset=1)
        .select(pl.col(["CB", "CB_idx_orig", "CB_idx_new"]))
    ).collect(streaming=True)

    # Write filtered list of CBs to barcodes filename.
    with xopen.xopen(barcodes_tsv_out_filename, "wb") as fh_barcodes_tsv:
        filtered_cb_idx_orig_to_cb_idx_new_df.select(pl.col("CB")).write_csv(
            fh_barcodes_tsv,
            has_header=False,
        )

    # Correct CB idx column in MatrixMarket matrix dataframe by mapping the original
    # CB indices to the filtered list of CB indices.
    matrix_mtx_out_df = (
        matrix_mtx_df.lazy()
        .join(
            filtered_cb_idx_orig_to_cb_idx_new_df.lazy(),
            how="left",
            on="CB_idx_orig",
        )
        .select(
            pl.col("feature_idx"),
            pl.col("CB_idx_new"),
            pl.col("value"),
        )
    ).collect(streaming=True)

    # Write output MatrixMarket matrix file.
    with xopen.xopen(mtx_out_filename, "wb") as fh_matrix_mtx:
        # Write MatrixMarket matrix header from original file.
        fh_matrix_mtx.write(matrix_mtx_header)

        # Write MatrixMarket matrix dataframe with corrected CB indices.
        matrix_mtx_out_df.write_csv(
            fh_matrix_mtx,
            sep=" ",
            has_header=False,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Filter STAR solo raw count matrix files and retain only those CBs which appear at least once.",
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input_dir",
        action="store",
        type=str,
        required=True,
        help="Input dir with STAR solo count matrix files (barcodes.tsv, features.tsv, matrix.mtx).",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        action="store",
        type=str,
        required=True,
        help="Output dir to which filtered STAR solo count matrix files (barcodes.tsv, features.tsv, matrix.mtx) will be written.",
    )

    args = parser.parse_args()

    if os.path.realpath(args.input_dir) == os.path.realpath(args.output_dir):
        print(
            "Error: Input and output dir can not be the same.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.exists(args.input_dir):
        print(
            f'Error: Input directory "{args.input_dir}" does not exist.',
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.exists(args.output_dir):
        print(
            f'Error: Output directory "{args.output_dir}" does not exist.',
            file=sys.stderr,
        )
        sys.exit(1)

    barcodes_tsv_filename = os.path.join(args.input_dir, "barcodes.tsv")
    features_tsv_filename = os.path.join(args.input_dir, "features.tsv")
    matrix_mtx_filename = os.path.join(args.input_dir, "matrix.mtx")

    if not os.path.exists(barcodes_tsv_filename):
        if os.path.exists(barcodes_tsv_filename + ".gz"):
            barcodes_tsv_filename += ".gz"
        else:
            print(
                f'Error: Barcode TSV file "{barcodes_tsv_filename}" or "{barcodes_tsv_filename + ".gz"}" does not exist.',
                file=sys.stderr,
            )
            sys.exit(1)

    if not os.path.exists(features_tsv_filename):
        if os.path.exists(features_tsv_filename + ".gz"):
            features_tsv_filename += ".gz"
        else:
            print(
                f'Error: Features TSV file "{features_tsv_filename}" or "{features_tsv_filename + ".gz"}" does not exist.',
                file=sys.stderr,
            )
            sys.exit(1)

    if not os.path.exists(matrix_mtx_filename):
        if os.path.exists(matrix_mtx_filename + ".gz"):
            matrix_mtx_filename += ".gz"
        else:
            print(
                f'Error: Matrix mtx file "{matrix_mtx_filename}" or "{matrix_mtx_filename + ".gz"}" does not exist.',
                file=sys.stderr,
            )
            sys.exit(1)

    barcodes_tsv_out_filename = os.path.join(
        args.output_dir, os.path.basename(barcodes_tsv_filename)
    )
    features_tsv_out_filename = os.path.join(
        args.output_dir, os.path.basename(features_tsv_filename)
    )
    matrix_mtx_out_filename = os.path.join(
        args.output_dir, os.path.basename(matrix_mtx_filename)
    )

    # Read barcodes TSV file and create a Polars DataFrame with original CB indices
    # and CB names.
    cb_idx_orig_to_cb_df = get_orig_cb_idx_to_cb(barcodes_tsv_filename)

    # Read MatrixMarket matrix file as Polars DataFrame and a header (first 3 lines).
    matrix_mtx_df, matrix_mtx_header = read_mtx(matrix_mtx_filename)

    #   - Filter CBs from barcodes TSV file, so only those that are available in
    #     at least once in the MatrixMarket matrix file are retained.
    #   - Write filtered CBs to output barcodes TSV file.
    #   - Write output MatrixMarket matrix file with corrected CB indices.
    write_filtered_barcodes_and_matrix_mtx(
        cb_idx_orig_to_cb_df,
        matrix_mtx_df,
        matrix_mtx_header,
        barcodes_tsv_out_filename,
        matrix_mtx_out_filename,
    )

    # Copy features.tsv file to output dir without modificatons.
    shutil.copy(features_tsv_filename, features_tsv_out_filename)


if __name__ == "__main__":
    main()
