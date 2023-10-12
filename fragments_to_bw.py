#!/usr/bin/env python3

import argparse
import numpy as np
import pyBigWig
import polars as pl
import os
import gzip
import pyarrow as pa
import pyarrow.csv
import numba

from typing import Literal
from pathlib import Path


def get_chromosome_sizes(chrom_sizes_filename: str):
    chrom_sizes = {}

    with open(chrom_sizes_filename, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            chrom, chrom_size = line.rstrip("\n").split("\t")[0:2]

            chrom_sizes[chrom] = int(chrom_size)

    return chrom_sizes


def normalise_filepath(path: str | Path, check_not_directory: bool = True) -> str:
    """
    Create a string path, expanding the home directory if present.

    """
    path = os.path.expanduser(path)
    if check_not_directory and os.path.exists(path) and os.path.isdir(path):
        raise IsADirectoryError(f"Expected a file path; {path!r} is a directory")
    return path


def read_bed_to_polars_df(
    bed_filename: str,
    engine: str | Literal["polars"] | Literal["pyarrow"] = "pyarrow",
    min_column_count: int = 3,
) -> pl.DataFrame:
    """
    Read BED file to a Polars DataFrame.

    Parameters
    ----------
    bed_filename
        BED filename.
    engine
        Use Polars or pyarrow to read the BED file (default: `pyarrow`).
    min_column_count
        Minimum number of required columns needed in BED file.

    Returns
    -------
    Polars DataFrame with BED entries.

    See Also
    --------
    pycisTopic.fragments.read_fragments_to_polars_df

    Examples
    --------
    Read BED file to Polars DataFrame with pyarrow engine.
    >>> bed_df_pl = read_bed_to_polars_df("test.bed", engine="pyarrow")

    Read BED file to Polars DataFrame with pyarrow engine and require that the BED
    file has at least 4 columns.
    >>> bed_with_at_least_4_columns_df_pl = read_bed_to_polars_df(
    ...     "test.bed",
    ...     engine="pyarrow",
    ...     min_column_count=4,
    ... )

    """
    bed_column_names = (
        "Chromosome",
        "Start",
        "End",
        "Name",
        "Score",
        "Strand",
        "ThickStart",
        "ThickEnd",
        "ItemRGB",
        "BlockCount",
        "BlockSizes",
        "BlockStarts",
    )

    bed_filename = normalise_filepath(bed_filename)

    # Set the correct open function, depending upon if the fragments BED file is gzip
    # compressed or not.
    open_fn = gzip.open if bed_filename.endswith(".gz") else open

    skip_rows = 0
    column_count = 0
    with open_fn(bed_filename, "rt") as bed_fh:
        for line in bed_fh:
            # Remove newlines and spaces.
            line = line.strip()

            if not line or line.startswith("#"):
                # Count number of empty lines and lines which start with a comment
                # before the actual data.
                skip_rows += 1
            else:
                # Get number of columns from the first real BED entry.
                column_count = len(line.split("\t"))

                # Stop reading the BED file.
                break

    if column_count < min_column_count:
        raise ValueError(
            f"BED file needs to have at least {min_column_count} columns. "
            f'"{bed_filename}" contains only {column_count} columns.'
        )

    # Enable global string cache so categorical columns from multiple Polars DataFrames
    # can be joined later, if necessary.
    pl.enable_string_cache(True)

    if engine == "polars":
        # Read BED file with Polars.
        bed_df_pl = pl.read_csv(
            bed_filename,
            has_header=False,
            skip_rows=skip_rows,
            separator="\t",
            use_pyarrow=False,
            new_columns=bed_column_names[:column_count],
            dtypes={
                bed_column: dtype
                for bed_column, dtype in {
                    "Chromosome": pl.Categorical,
                    "Start": pl.Int32,
                    "End": pl.Int32,
                    "Name": pl.Categorical,
                    "Strand": pl.Categorical,
                }.items()
                if bed_column in bed_column_names[:column_count]
            },
        )
    elif engine == "pyarrow":
        # Read BED file with pyarrow.
        bed_df_pl = pl.from_arrow(
            pa.csv.read_csv(
                bed_filename,
                read_options=pa.csv.ReadOptions(
                    use_threads=True,
                    skip_rows=skip_rows,
                    column_names=bed_column_names[:column_count],
                ),
                parse_options=pa.csv.ParseOptions(
                    delimiter="\t",
                    quote_char=False,
                    escape_char=False,
                    newlines_in_values=False,
                ),
                convert_options=pa.csv.ConvertOptions(
                    column_types={
                        "Chromosome": pa.dictionary(pa.int32(), pa.large_string()),
                        "Start": pa.int32(),
                        "End": pa.int32(),
                        "Name": pa.dictionary(pa.int32(), pa.large_string()),
                        "Strand": pa.dictionary(pa.int32(), pa.large_string()),
                    },
                ),
            )
        )
    else:
        raise ValueError(
            f'Unsupported engine value "{engine}" (allowed: ["polars", "pyarrow"]).'
        )

    return bed_df_pl


def read_fragments_to_polars_df(
    fragments_bed_filename: str,
    engine: str | Literal["polars"] | Literal["pyarrow"] = "pyarrow",
) -> pl.DataFrame:
    """
    Read fragments BED file to a Polars DataFrame.

    If fragments don't have a Score column, a Score columns is created by counting
    the number of fragments with the same chromosome, start, end and CB.

    Parameters
    ----------
    fragments_bed_filename
        Fragments BED filename.
    engine
        Use Polars or pyarrow to read the fragments BED file (default: pyarrow).

    Returns
    -------
    Polars DataFrame with fragments.

    See Also
    --------
    pycisTopic.fragments.read_bed_to_polars_df

    Examples
    --------
    Read gzipped fragments BED file to a Polars DataFrame.
    >>> fragments_df_pl = read_fragments_to_polars_df(
    ...     fragments_bed_filename="fragments.tsv.gz",
    ... )

    Read uncompressed fragments BED file to a Polars DataFrame.
    >>> fragments_df_pl = read_fragments_to_polars_df(
    ...     fragments_bed_filename="fragments.tsv",
    ... )

    """
    fragments_df_pl = read_bed_to_polars_df(
        bed_filename=fragments_bed_filename,
        engine=engine,
        min_column_count=4,
    )

    # If no score is provided or score column is ".", generate a score column with the
    # number of fragments which have the same chromosome, start, end and CB.
    if (
        "Score" not in fragments_df_pl.columns
        or fragments_df_pl.schema["Score"] == pl.Utf8
    ):
        fragments_df_pl = fragments_df_pl.groupby(
            ["Chromosome", "Start", "End", "Name"]
        ).agg(pl.count().cast(pl.Int32()).alias("Score"))
    else:
        fragments_df_pl = fragments_df_pl.with_columns(pl.col("Score").cast(pl.Int32()))

    return fragments_df_pl


@numba.njit
def calculate_depth(chrom_size, starts, ends):
    # Initialize array for current chromosome to store the depth per basepair.
    chrom_depth = np.zeros(chrom_size, dtype=np.uint32)

    # Require same number of start and end positions.
    assert starts.shape[0] == ends.shape[0]

    for i in range(starts.shape[0]):
        # Add 1 depth for each basepair in the current fragment.
        chrom_depth[starts[i] : ends[i]] += numba.uint32(1)

    return chrom_depth


@numba.njit
def collapse_consecutive_values(X):
    # Length.
    n = X.shape[0]

    # Create idx array with enough space to store all possible indices
    # in case there are no consecutive values in the input that are
    # the same.
    idx = np.empty(n + 1, dtype=np.uint32)

    # First index position will always be zero.
    idx[0] = 0
    # Next index postion to fill in in idx.
    j = 1

    # Loop over the whole input array and store indices for only those
    # positions for which the previous value was different.
    for i in numba.prange(1, n):
        if X[i - 1] == X[i]:
            continue

        # Current value is different from previous value, so store the index
        # position of the current index i in idx[j].
        idx[j] = i

        # Increase next index position to fill in in idx.
        j += 1

    # Store length of input as last value. Needed to calculate the number of times
    # the last consecutive value gets repeated.
    idx[j] = n

    # Get all consecutive different values from the input.
    values = X[idx[:j]].astype(np.float32)

    # Calculate the number of times each value in values gets consecutively
    # repeated in the input.
    lengths = idx[1 : j + 1] - idx[:j]

    # Restrict indices array to same number of element than the values and lentgths
    # arrays.
    idx = idx[:j].copy()

    # To reconstruct the original input: X == values.repeat(lenghts)
    return idx, values, lengths


def fragments_to_bw(
    fragments_df: pl.DataFrame,
    chrom_sizes: dict[str, int],
    bw_filename: str,
    normalize: bool = True,
    scaling_factor: float = 1.0,
):
    chrom_arrays = {}

    with pyBigWig.open(bw_filename, "wb") as bw:
        print("Add chromosome sizes to bigWig header")
        bw.addHeader([(chrom, chrom_size) for chrom, chrom_size in chrom_sizes.items()])

        for chrom, chrom_size in chrom_sizes.items():
            chrom_arrays[chrom] = np.zeros(chrom_size, dtype=np.uint32)

        no_fragments = 0

        print(f"Number of fragments: {fragments_df.height}")
        print("Split fragments df by chromosome")
        per_chrom_fragments_dfs = fragments_df.partition_by("Chromosome", as_dict=True)

        print("Calculate depth per chromosome:")
        for chrom in per_chrom_fragments_dfs:
            print(f"  - {chrom} ...")
            if chrom not in chrom_sizes:
                print(f"    Skipping {chrom} as it is not in chrom sizes file.")
                continue
            starts, ends = (
                per_chrom_fragments_dfs[chrom].select(["Start", "End"]).to_numpy().T
            )
            chrom_arrays[chrom] = calculate_depth(chrom_sizes[chrom], starts, ends)
            no_fragments += per_chrom_fragments_dfs[chrom].height

        # Calculate RPM scaling factor.
        rpm_scaling_factor = no_fragments / 1_000_000.0

        print(
            "Compact depth array per chromosome (make ranges for consecutive the same values and remove zeros):"
        )
        for chrom, chrom_size in chrom_sizes.items():
            print(f"  - Compact {chrom} ...")
            idx, values, lengths = collapse_consecutive_values(chrom_arrays[chrom])
            non_zero_idx = np.flatnonzero(values)

            if non_zero_idx.shape[0] == 0:
                # Skip chromosomes with no depth > 0.
                continue

            # Select only consecutive different values and calculate start and end coordinates
            # (in BED format) for each of those ranges.
            chroms = np.repeat(chrom, len(non_zero_idx))
            starts = idx[non_zero_idx]
            ends = idx[non_zero_idx] + lengths[non_zero_idx]
            values = values[non_zero_idx]

            if normalize:
                values = values / rpm_scaling_factor * scaling_factor
            elif scaling_factor != 1.0:
                values *= scaling_factor

            print(f"  - Write {chrom} to bigWig ...")
            bw.addEntries(chroms=chroms, starts=starts, ends=ends, values=values)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate genome coverage for fragments and write result to bigWig file."
    )

    parser.add_argument(
        "-c",
        "--chrom",
        dest="chrom_sizes_filename",
        action="store",
        type=str,
        required=True,
        help="Filename with chromosome sizes (*.chrom.sizes, *.fa.fai).",
    )

    parser.add_argument(
        "-i",
        "--frag",
        dest="fragments_filename",
        action="store",
        type=str,
        required=True,
        help="Fragments TSV file for which to calculate genome coverage.",
    )

    parser.add_argument(
        "-o",
        "--bw",
        dest="bigwig_filename",
        action="store",
        type=str,
        required=True,
        help="BigWig filename to which the genome coverage data will be written.",
    )

    parser.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store",
        type=str,
        required=False,
        choices={"yes", "no"},
        default="yes",
        help='Normalize genome coverage data. Default: "yes".',
    )

    parser.add_argument(
        "-s",
        "--scaling",
        dest="scaling_factor",
        action="store",
        type=float,
        required=False,
        default=1.0,
        help='Scaling factor for genome coverage data. Default: "1.0".',
    )

    args = parser.parse_args()

    chrom_sizes = get_chromosome_sizes(args.chrom_sizes_filename)
    fragments_df = read_fragments_to_polars_df(args.fragments_filename)

    fragments_to_bw(
        fragments_df, chrom_sizes, args.bigwig_filename, args.normalize == "yes", args.scaling_factor
    )


if __name__ == "__main__":
    main()
