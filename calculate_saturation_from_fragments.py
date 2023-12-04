#!/usr/bin/env python

import argparse
import gzip
import logging
import os

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import polars as pl

__author__ = "Swann Floc Hlay"
__contributors__ = "Gert Hulselmans"
__version__ = "v0.2.0"


FORMAT = "%(asctime)-15s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)
logger = logging.getLogger("calculate_saturation_from_fragments")


# Fractions for which to sample from the fragments file:
#     [
#         float(f'{sampling_fraction:0.2f}')
#         for sampling_fraction in list(np.concatenate((np.arange(0.0, 0.5, 0.1), np.arange(0.5, 0.9, 0.05), np.arange(0.9, 1.01, 0.02))))
#     ]
sampling_fractions_default = [
    0.0, 0.1, 0.2, 0.3, 0.4,
    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
    0.9, 0.92, 0.94, 0.96, 0.98, 1.0
]  # fmt: skip


### initialise function and classes


def read_bc_and_counts_from_fragments_file(fragments_bed_filename: str) -> pl.DataFrame:
    """
    Read cell barcode (column 4) and counts per fragment (column 5) from fragments BED file.

    Cell barcodes will appear more than once as they have counts per fragment, but as
    the fragment locations are not needed, they are not returned.

    Parameters
    ----------
    fragments_bed_filename
        Fragments BED filename.

    Returns
    -------
    Polars dataframe with cell barcode and count per fragment (column 4 and 5 of BED file).

    """
    # Set the correct open function depending if the fragments BED file is gzip compressed or not.
    open_fn = gzip.open if fragments_bed_filename.endswith(".gz") else open

    skip_rows = 0
    nbr_columns = 0

    with open_fn(fragments_bed_filename, "rt") as fragments_bed_fh:
        for line in fragments_bed_fh:
            # Remove newlines and spaces.
            line = line.strip()

            if not line or line.startswith("#"):
                # Count number of empty lines and lines which start with a comment before the actual data.
                skip_rows += 1
            else:
                # Get number of columns from the first real BED entry.
                nbr_columns = len(line.split("\t"))

                # Stop reading the BED file.
                break

    if nbr_columns < 5:
        raise ValueError(
            "Fragments BED file needs to have at least 5 columns. "
            f'"{fragments_bed_filename}" contains only {nbr_columns} columns.'
        )

    # Read cell barcode (column 4) and counts (column 5) per fragment from fragments BED file.
    fragments_df = pl.read_csv(
        fragments_bed_filename,
        has_header=False,
        skip_rows=skip_rows,
        separator="\t",
        use_pyarrow=False,
        columns=["column_1", "column_2", "column_3", "column_4", "column_5"],
        new_columns=["Chromosome", "Start", "End", "CellBarcode", "FragmentCount"],
        dtypes=[pl.Categorical, pl.UInt32, pl.UInt32, pl.Categorical, pl.UInt32],
    )

    return fragments_df


def MM(x, Vmax, Km):
    """Define the Michaelis-Menten Kinetics model that will be used for the model fitting."""
    y = Vmax * x / (Km + x) if Vmax > 0 and Km > 0 else 10000000000.0
    return y


# sub-sampling function
def sub_sample_fragments(
    fragments_df,
    n_reads,
    cbs=None,
    min_uniq_frag=200,
    sampling_fractions=sampling_fractions_default,
    stats_tsv_filename="sampling_stats.tsv",
):
    sampling_fractions_length = len(sampling_fractions)

    # Initialize dataframe for storing all statistics results.
    stats_df = pd.DataFrame(
        {
            "total_unique_frag_count": np.zeros(sampling_fractions_length, np.uint32),
            "total_frag_count": np.zeros(sampling_fractions_length, np.uint32),
            "mean_frag_per_bc": np.zeros(sampling_fractions_length, np.float64),
            "median_uniq_frag_per_bc": np.zeros(sampling_fractions_length, np.float64),
            "cell_barcode_count": np.zeros(sampling_fractions_length, np.uint32),
            "total_reads": np.array(sampling_fractions, np.float64) * n_reads,
        },
        index=pd.Index(data=np.array(sampling_fractions), name="sampling_fraction"),
    )

    if cbs is not None:
        if isinstance(cbs, pl.DataFrame) and cbs.columns == ["CellBarcode"]:
            # Selected cell barcodes are provided in the correct format.
            selected_cbs_df = cbs
        else:
            # Convert selected cell barcodes provided as input to a Polars DataFrame.
            selected_cbs_df = pl.DataFrame(
                [
                    pl.Series("CellBarcode", cbs, dtype=pl.Categorical),
                ]
            )
    else:
        # Get all cell barcodes which have more than min_uniq_frag fragments from the
        # fragments file.
        selected_cbs_df = (
            fragments_df.group_by("CellBarcode")
            .agg(pl.col("FragmentCount").count().alias("nbr_frags_per_CBs"))
            .filter(pl.col("nbr_frags_per_CBs") > min_uniq_frag)
        )

    # Count all selected cell barcodes.
    nbr_selected_cbs = selected_cbs_df.height

    if 1.0 in sampling_fractions:
        # As there is no need to sample when sampling fraction is 100%,
        # the median number of unique fragments per barcode can be
        # calculated much more efficiently on the original fragments
        # file dataframe with counts than the expanded one, which is
        # needed when sampling is required.

        logger.info("Calculate statistics for sampling fraction 100.0%.")

        logger.info("Keep fragments with selected cell barcodes.")
        fragments_for_selected_bcs_df = selected_cbs_df.join(
            fragments_df,
            on="CellBarcode",
            how="left",
        )

        logger.info("Calculate total number of fragments.")

        stats_df.loc[1.0, "total_frag_count"] = fragments_for_selected_bcs_df.select(
            pl.col("FragmentCount").sum()
        ).item()

        logger.info("Calculate total unique number of fragments.")

        stats_df.loc[1.0, "total_unique_frag_count"] = (
            fragments_for_selected_bcs_df.group_by(
                ["CellBarcode", "Chromosome", "Start", "End"]
            )
            .agg([pl.first("Start").alias("Start_tmp")])
            .select(pl.count())
            .item()
        )

        logger.info(
            "Calculate mean number of fragments per cell barcode and median number of "
            "unique fragments per cell barcode."
        )
        stats_df_pl = (
            fragments_for_selected_bcs_df.group_by("CellBarcode")
            .agg(
                [
                    pl.col("FragmentCount").sum().alias("MeanFragmentsPerCB"),
                    pl.count().alias("UniqueFragmentsPerCB"),
                ]
            )
            .select(
                [
                    pl.col("MeanFragmentsPerCB").mean(),
                    pl.col("UniqueFragmentsPerCB").median(),
                ]
            )
        )

        stats_df.loc[1.0, "mean_frag_per_bc"] = stats_df_pl["MeanFragmentsPerCB"][0]
        stats_df.loc[1.0, "median_uniq_frag_per_bc"] = stats_df_pl[
            "UniqueFragmentsPerCB"
        ][0]
        stats_df.loc[1.0, "cell_barcode_count"] = nbr_selected_cbs

        # Delete dataframe to free memory.
        del fragments_for_selected_bcs_df

    # Create dataframe where each row contains one fragment:
    #   - Original dataframe has a count per fragment with the same cell barcode.
    #   - Create a row for each count, so we can sample fairly afterward.
    logger.info("Create dataframe with all fragments (for sampling).")
    fragments_all_df = fragments_df.with_columns(
        pl.col("FragmentCount").repeat_by(pl.col("FragmentCount"))
    ).explode("FragmentCount")

    # Delete input dataframe to free memory.
    del fragments_df

    for sampling_fraction in sampling_fractions:
        if sampling_fraction == 0.0:
            # All statistics are zero and already set when the stats_df dataframe is
            # created.
            continue
        elif sampling_fraction == 1.0:
            # Statistics for 100% sampling are already calculated as there is no need
            # to have the fragments_all_df dataframe as no sampling is needed.
            # This avoids the need to use the expensive group_by operations for the
            # calculations of the median number of unique fragments per barcode.
            continue

        logger.info(
            "Calculate statistics for sampling fraction "
            f"{round(sampling_fraction * 100, 1)}%."
        )

        # Sample x% from all fragments (with duplicates) and keep fragments which have
        # selected cell barcodes.
        logger.info(
            f"Sample {round(sampling_fraction * 100, 1)}% from all fragments and keep "
            "fragments with selected cell barcodes."
        )
        fragments_sampled_for_selected_cb_df = selected_cbs_df.join(
            fragments_all_df.sample(fraction=sampling_fraction),
            on="CellBarcode",
            how="inner",
        )

        logger.info("Calculate total number of fragments.")
        stats_df.loc[
            sampling_fraction, "total_frag_count"
        ] = fragments_sampled_for_selected_cb_df.height

        logger.info("Calculate total unique number of fragments.")
        stats_df.loc[sampling_fraction, "total_unique_frag_count"] = (
            fragments_sampled_for_selected_cb_df.group_by(
                ["CellBarcode", "Chromosome", "Start", "End"]
            )
            .agg([pl.first("Start").alias("Start_tmp")])
            .select(pl.count())
            .item()
        )

        logger.info("Calculate mean number of fragments per cell barcode.")
        stats_df.loc[sampling_fraction, "mean_frag_per_bc"] = (
            fragments_sampled_for_selected_cb_df.select(
                [pl.col("CellBarcode"), pl.col("FragmentCount")]
            )
            .group_by("CellBarcode")
            .agg([pl.count("FragmentCount").alias("FragmentsPerCB")])
            .select([pl.col("FragmentsPerCB").mean().alias("MeanFragmentsPerCB")])
            .item()
        )

        logger.info("Calculate median number of unique fragments per cell barcode.")
        stats_df.loc[sampling_fraction, "median_uniq_frag_per_bc"] = (
            fragments_sampled_for_selected_cb_df.group_by(
                ["CellBarcode", "Chromosome", "Start", "End"]
            )
            .agg([pl.col("FragmentCount").first().alias("FragmentCount")])
            .select([pl.col("CellBarcode"), pl.col("FragmentCount")])
            .group_by("CellBarcode")
            .agg(pl.col("FragmentCount").count().alias("UniqueFragmentsPerCB"))
            .select(pl.col("UniqueFragmentsPerCB").median())
            .item()
        )

        logger.info("Calculate number of cell barcodes.")
        stats_df.loc[sampling_fraction, "cell_barcode_count"] = (
            fragments_sampled_for_selected_cb_df.group_by(["CellBarcode"])
            .agg([pl.first("Start").alias("Start_tmp")])
            .select(pl.count())
            .item()
        )

        # Delete dataframe to free memory.
        del fragments_sampled_for_selected_cb_df

    logger.info("Add extra statistics.")
    stats_df["mean_reads_per_barcode"] = (
        stats_df["total_reads"] / stats_df["cell_barcode_count"]
    )

    stats_df["duplication_rate"] = (
        stats_df["total_frag_count"] - stats_df["total_unique_frag_count"]
    ) / stats_df["total_frag_count"]

    if 0.0 in sampling_fractions:
        stats_df.loc[0.0, "mean_reads_per_barcode"] = 0.0
        stats_df.loc[0.0, "duplication_rate"] = 0.0

    logger.info(f'Saving statistics in "{stats_tsv_filename}".')
    stats_df.to_csv(stats_tsv_filename, sep="\t")

    return stats_df


# Format axis labels
def format_axis(text):
    format_dict = {
        "mean_frag_per_bc": "Mean fragments / barcode",
        "median_uniq_frag_per_bc": "Median unique fragments / barcode",
        "total_frag_count": "Total fragment count",
    }
    for i, j in format_dict.items():
        text = text.replace(i, j)
    return text


# MM-fit function
def fit_MM(
    stats_df,
    percentages=[0.3, 0.6, 0.9],
    saturation_plot_filename="saturation.png",
    x_axis="total_frag_count",
    y_axis="median_uniq_frag_per_bc",
):
    # select x/y data from MM fit from subsampling stats
    x_data = np.array(stats_df.loc[:, x_axis]) / 10**6
    y_data = np.array(stats_df.loc[:, y_axis])
    # fit to MM function
    best_fit_ab, covar = curve_fit(MM, x_data, y_data, bounds=(0, +np.inf))
    # expand fit space
    x_fit = np.linspace(0, int(np.max(x_data) * 100), num=500)
    y_fit = MM(x_fit, *(best_fit_ab))
    # impute maximum saturation to plot as 95% of y_max
    y_val = best_fit_ab[0] * 0.95
    # subset x_fit space if bigger then y_val
    if y_val < max(y_fit):
        x_coef = np.where(y_fit >= y_val)[0][0]
        x_fit = x_fit[0:x_coef]
        y_fit = y_fit[0:x_coef]
    # plot model
    plt.plot(x_fit, MM(x_fit, *best_fit_ab), label="fitted", c="black", linewidth=1)
    # plot raw data
    plt.scatter(x=x_data, y=y_data, c="crimson", s=10)
    # mark curent saturation
    x_idx = np.where(y_fit >= max(y_data))[0][0]
    x_coef = x_fit[x_idx]
    y_coef = y_fit[x_idx]
    plt.plot([x_coef, x_coef], [0, y_coef], linestyle="--", c="crimson")
    plt.plot([0, x_coef], [y_coef, y_coef], linestyle="--", c="crimson")
    plt.text(
        x=x_fit[-1],
        y=y_coef,
        s=str(round(100 * max(y_data) / best_fit_ab[0])) + f"% {x_coef:.2f}",
        c="crimson",
        ha="right",
        va="bottom",
    )
    # plot percentaged values
    for perc in percentages:
        # Find read count for percent saturation
        y_val = best_fit_ab[0] * perc
        # Find closest match in fit
        if max(y_fit) > y_val:
            x_idx = np.where(y_fit >= y_val)[0][0]
            x_coef = x_fit[x_idx]
            y_coef = y_fit[x_idx]
            # Draw vline
            plt.plot([x_coef, x_coef], [0, y_coef], linestyle="--", c="grey")
            # Draw hline
            plt.plot([0, x_coef], [y_coef, y_coef], linestyle="--", c="grey")
            # Plot imputed read count
            plt.text(
                x=x_fit[-1],
                y=y_coef,
                s=str(round(100 * perc)) + f"% {x_coef:.2f}",
                c="grey",
                ha="right",
                va="bottom",
            )
    # save figure
    plt.xlabel(format_axis(x_axis) + " (millions)")
    plt.ylabel(format_axis(y_axis))
    plt.title(os.path.basename(saturation_plot_filename).replace(".saturation.png", ""))
    plt.savefig(saturation_plot_filename)
    plt.close()


def main():
    sampling_fractions_default_str = ",".join(
        [str(x) for x in sampling_fractions_default]
    )

    parser = argparse.ArgumentParser(
        description="Infer saturation of scATAC from fragments file."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="fragments_input_bed_filename",
        action="store",
        type=str,
        required=True,
        help="Fragment input BED filename.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_prefix",
        action="store",
        type=str,
        required=True,
        help=(
            "Output prefix, which will contain PNG file with saturation curve and TSV"
            " file with summary of reads and additional reads needed to reach"
            " saturation specified by percentages."
        ),
    )
    parser.add_argument(
        "-n",
        "--raw-reads",
        dest="n_reads",
        action="store",
        type=str,
        required=True,
        help="Number of raw sequenced reads for this fragments file.",
    )
    parser.add_argument(
        "-c",
        "--cbs",
        dest="cbs_filename",
        type=str,
        help=(
            "Filename with list of selected cell barcodes. "
            "If not specified --min_frags_per_cb is used to get the list of "
            "selected cell barcodes"
        ),
        default=None,
    )
    parser.add_argument(
        "-p",
        "--percentages",
        dest="percentages",
        type=str,
        help=(
            "Comma separated list of decimal percentages to predict. Default: "
            '"0.3,0.6,0.9"'
        ),
        default="0.3,0.6,0.9",
    )
    parser.add_argument(
        "-m",
        "--min_frags_per_cb",
        dest="min_frags_per_cb",
        type=int,
        help="Minimum number of unique fragments per cell barcodes. Default: 200",
        default=200,
    )
    parser.add_argument(
        "-s",
        "--sampling_fractions",
        dest="sampling_fractions",
        type=str,
        help=(
            "Fractions at which to perform the sub-samplings. Default: "
            f'"{sampling_fractions_default_str}"'
        ),
        default=sampling_fractions_default_str,
    )

    parser.add_argument("-V", "--version", action="version", version=f"{__version__}")

    args = parser.parse_args()

    sampling_fractions = [float(x) for x in args.sampling_fractions.split(",")]
    percentages = [float(x) for x in args.percentages.split(",")]

    # Enable global string cache.
    pl.enable_string_cache()

    # Load fragments BED file.
    logger.info("Loading fragments BED file started.")
    fragments_df = read_bc_and_counts_from_fragments_file(
        args.fragments_input_bed_filename
    )
    logger.info("Loading fragments BED file finished.")

    if args.cbs_filename:
        logger.info("Loading selected cell barcodes.")
        cbs_df = pl.read_csv(
            args.cbs_filename,
            has_header=False,
            columns=["column_1"],
            new_columns=["CellBarcode"],
            dtypes=[pl.Categorical],
        )
    else:
        cbs_df = None

    # Sub-sample.
    stats_df = sub_sample_fragments(
        fragments_df=fragments_df,
        cbs=cbs_df,
        min_uniq_frag=args.min_frags_per_cb,
        sampling_fractions=sampling_fractions,
        stats_tsv_filename=f"{args.output_prefix}.sampling_stats.tsv",
    )

    logger.info("fit_MM.")
    # Fit'n'plot for total count.
    fit_MM(
        stats_df,
        percentages=percentages,
        saturation_plot_filename=f"{args.output_prefix}.saturation.png",
        x_axis="total_frag_count",
        y_axis="median_uniq_frag_per_bc",
    )
    logger.info("Finished.")


if __name__ == "__main__":
    main()
