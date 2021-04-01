#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import json
import pandas as pd
from scipy.optimize import curve_fit
from uncertainties import ufloat
from typing import Tuple
from pathlib import Path
import matplotlib.pylab as plt
import numpy as np
import argparse

__author__ = "Jasper Janssens"
__contributors__ = "Swan Floc’Hlay, Maxime De Waegeneer"
__version__ = "v0.2.0"
__contact__ = "jasper.janssens@kuleuven.be"

f"""
calculate_saturation.py ({__version__}): 
    Infer saturation of 10x scATAC/scRNA sample

    Requirements:
        - Python 3
        - Packages: pandas, numpy, scipy, matplotlib, uncertainties, argparse

    How to run:
        - Example:
            python calculate_saturation.py \
                --dir 10X_dir \
                --type RNA \
                --o output_dir 
    Input:
        - Output folder from CellRanger (tested with 3.1.0)
    
    Output:
        - 1 png with saturation/complexity curve
        - 1 tsv with summary of reads and additional reads needed for 50/75% 
"""

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="10x folder", required=True)
parser.add_argument("-t", "--type", help="RNA/ATAC", required=True)
parser.add_argument("-o", "--output", help="path to output", required=True)
parser.add_argument(
    "--percentages",
    metavar="percentages",
    type=str,
    help="Comma separated list of decimal percentages to predict",
    default="0.5,0.75",
)

args = parser.parse_args()

percentages = [float(x) for x in args.percentages.split(",")]

print(f"Calculating saturation curves for 10x {args.type} data...")
print(f"10x folder used: {args.dir}")

if not os.path.isdir(args.output):
    # print(f"{args.output} does not exist. Creating...", file=sys.stderr)
    os.mkdir(args.output)


if args.type == "ATAC":
    complexity_info_dir = (
        Path(args.dir)
        / "SC_ATAC_COUNTER_CS"
        / "SC_ATAC_COUNTER"
        / "_SC_ATAC_METRIC_COLLECTOR"
        / "ESTIMATE_LIBRARY_COMPLEXITY"
        / "fork0"
    )
    if not complexity_info_dir.exists():
        raise FileNotFoundError(
            f"The given 10x {args.type} folder {args.dir} does not exits."
        )

    a = os.scandir(path=complexity_info_dir)
    file_path = [x.name for x in a if x.name.startswith("join-")]
    complexity_info_path = (
        Path(args.dir)
        / "SC_ATAC_COUNTER_CS"
        / "SC_ATAC_COUNTER"
        / "_SC_ATAC_METRIC_COLLECTOR"
        / "ESTIMATE_LIBRARY_COMPLEXITY"
        / "fork0"
        / file_path[0]
        / "files"
        / "singlecell_complexity.json"
    )
    if not complexity_info_path.exists():
        raise FileNotFoundError(
            f"The complexity JSON file of the given 10x {args.type} folder {args.dir} does not exist."
        )

    summary_info_path = Path(args.dir) / "outs" / "summary.csv"
    if not summary_info_path.exists():
        raise FileNotFoundError(
            f"The suymmary info file of the given 10x {args.type} folder {args.dir} does not exist."
        )
    summary = pd.read_csv(summary_info_path)
    num_cells = int(summary["annotated_cells"])

elif args.type == "RNA":
    path = args.dir
    summary_info_path = (
        Path(args.dir)
        / "SC_RNA_COUNTER_CS"
        / "SC_RNA_COUNTER"
        / "SUMMARIZE_REPORTS"
        / "fork0"
        / "files"
        / "metrics_summary_csv.csv"
    )
    if not summary_info_path.exists():
        raise FileNotFoundError(
            f"The suymmary info file of the given 10x {args.type} folder {args.dir} does not exist."
        )
    summary = pd.read_csv(summary_info_path)
    num_cells = summary["Estimated Number of Cells"][0]
    num_cells = int(re.sub(",", "", num_cells))
    complexity_info_path = (
        Path(args.dir)
        / "SC_RNA_COUNTER_CS"
        / "SC_RNA_COUNTER"
        / "SUMMARIZE_REPORTS"
        / "fork0"
        / "files"
        / "metrics_summary_json.json"
    )


# Create output path
if args.type == "ATAC":
    complexity_info_path_stripping = re.sub(
        "/SC_ATAC_COUNTER_.*", "", str(complexity_info_path)
    )
    complexity_info_path_stripping = re.split("/", complexity_info_path_stripping)
    project_name = complexity_info_path_stripping[-1]
    output_path = Path(args.output) / f"{project_name}_complexity.png"
elif args.type == "RNA":
    project_name = "RNA"  # needs to be updated
    output_path = (Path(args.output) / f"{project_name}_saturation.png",)


def prepare_data(complexity_info_path: Path, assay_type):
    """
    Prepare the data for the model fit.
    Returns a Tuple2 with X coordinates and Y coordinates for the model fitting
    """
    # Open complexity file
    with open(complexity_info_path) as complexity_info_fh:
        complexity_info = json.load(complexity_info_fh)
    complexity_info_df = pd.DataFrame(complexity_info)

    if assay_type == "ATAC":
        # get x data
        x_data = np.array(complexity_info_df["total_depth"])
        # get y data
        y_data = np.array(complexity_info_df["unique"])
        return x_data, y_data
    elif assay_type == "RNA":
        subsampled_columns = [
            x for x in complexity_info_df.columns if "multi_raw_rpc_" in x
        ]
        saturation_data = complexity_info_df[subsampled_columns].copy()
        saturation_data = pd.DataFrame(saturation_data.max())
        saturation_data.index = [
            np.float(
                re.sub(
                    "multi_raw_rpc_", "", re.sub("_subsampled_duplication_frac", "", x)
                )
            )
            for x in subsampled_columns
        ]
        saturation_data = saturation_data.loc[saturation_data[0] != 0].copy()
        saturation_data = saturation_data.sort_values(by=0)
        saturation_data = saturation_data.reset_index()
        saturation_data.columns = [0, 1]
        # get x data
        x_data = np.array(saturation_data[0])
        # get y data
        y_data = np.array(saturation_data[1])
        return x_data, y_data

    else:
        raise Exception(f"Not a valid assay type: {assay_type}")


# Define the Michaelis-Menten Kinetics model that will be used for the model fitting
def MM(x, Vmax, Km):
    if Vmax > 0 and Km > 0:
        y = (Vmax * x) / (Km + x)
    else:
        y = 1e10
    return y


def fit_model(model, x_data, y_data):
    """
    Fit the given model through the data given by x_data and y_data
    Returns:
    - The best model fit
    - The best model fit params: Tuple(a, b, sigma_ab) (along with their standard deviations)
    - R2 (R aquared)
    """
    best_fit_ab, covar = curve_fit(model, x_data, y_data)
    sigma_ab = np.sqrt(np.diagonal(covar))

    # get parameter value
    a = ufloat(best_fit_ab[0], sigma_ab[0])
    b = ufloat(best_fit_ab[1], sigma_ab[1])
    best_fit_params = a, b, sigma_ab
    text_res = "Best fit parameters:\na = {}\nb = {}".format(a, b)
    print(text_res)

    # get rsquared value
    residuals = y_data - model(x_data, *(best_fit_ab))
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r_sq = 1 - (ss_res / ss_tot)
    # Return the best model fit, the model parameters and the rsquared value
    return best_fit_ab, best_fit_params, r_sq


def drawline(
    perc,
    coef,
    y_fit,
    x_fit,
    assay_type,
    x_max=None,
    color="k",
    linestyle="--",
):

    # get max value possible
    y_max = coef[0]
    if assay_type == "ATAC":
        y_val = y_max * perc
        y_val = int(y_val)
    elif assay_type == "RNA":
        y_val = perc
        y_val = np.float(y_val)

    # find x value at which fitted line surpasses y_val
    x_coef = np.where(y_fit >= y_val)[0][0]
    x_coef = x_fit[x_coef]

    # draw a vertical line at x_coef
    if assay_type == "ATAC":
        plt.axvline(x=x_coef, ymin=0, ymax=perc, color=color, linestyle=linestyle)
    elif assay_type == "RNA":
        plt.axvline(
            x=x_coef, ymin=0, ymax=perc / 1.05, color=color, linestyle=linestyle
        )

    # draw a horizontal line at y_val
    if not x_max:
        x_max = np.max(x_fit)
    plt.axhline(y=y_val, xmin=0, xmax=x_coef / x_max, color=color, linestyle=linestyle)
    return x_coef


def plot_saturation_curve(
    model,
    model_fit,
    model_fit_params,
    model_fit_r_sq: float,
    x_data: np.array,
    assay_type: str,
    output_path: Path,
):
    """
    Plot saturation curve and returns pd.DataFrame, reeds needed and saturation percentage along with mean reads per cell (X coefficient)
    """

    model_param_a, model_param_b, model_sigma_ab = model_fit_params

    plt.figure(figsize=(10, 7))

    # Plot the model
    if assay_type == "ATAC":
        x_fit = np.linspace(0, int(np.max(x_data) * 100), 100000)
    elif assay_type == "RNA":
        x_fit = np.linspace(0, int(np.max(x_data) * 2000), 100000)
    y_fit = model(x_fit, *(model_fit))

    # Get max value possible and highest we want to plot
    y_max = model_fit[0]
    y_val = y_max * 0.9

    # Find x value at which fitted line surpasses y_val
    x_coef = np.where(y_fit >= y_val)[0][0]
    x_fit = x_fit[0:x_coef]
    y_fit = y_fit[0:x_coef]
    y_max, x_max = np.max(y_fit), np.max(x_fit)

    # Plot model
    plt.plot(
        x_fit,
        model(x_fit, *model_fit),
        label="a = " + str(model_param_a),
        c="k",
        linewidth=3,
    )
    plt.plot(
        x_fit,
        model(x_fit, *model_fit),
        label="b = " + str(model_param_b),
        c="k",
        linewidth=3,
    )

    bound_upper = model(x_fit, *(model_fit + model_sigma_ab))
    bound_lower = model(x_fit, *(model_fit - model_sigma_ab))

    # Plotting the confidence intervals
    plt.fill_between(x_fit, bound_lower, bound_upper, color="gray", alpha=0.4)

    # Plot raw data
    plt.scatter(x=x_data, y=y_data, c="red")

    # Plot different levels
    # Loop over different percentages + save results in dataframe
    reads_needed = {}
    for perc in percentages:
        x_coef = drawline(perc, model_fit, y_fit, x_fit, assay_type=assay_type)
        # update summary
        reads_needed[perc] = x_coef
        plt.text(
            x=0.6 * x_max,
            y=perc * model_fit[0],
            s=str(int(perc * 100)) + "% = " + str(int(x_coef)) + " reads",
            ha="left",
            fontsize=18,
            wrap=True,
        )

    saturation_pct = np.max(y_data) / model_fit[0]
    mean_reads_per_cell = drawline(
        perc=saturation_pct,
        coef=model_fit,
        y_fit=y_data,
        x_fit=x_data,
        assay_type=args.type,
        x_max=x_max,
        color="r",
    )
    # update summary
    reads_needed["now"] = mean_reads_per_cell
    plt.text(
        x=0.6 * x_max,
        y=saturation_pct * model_fit[0],
        s=str(int(saturation_pct * 100))
        + "% = "
        + str(int(mean_reads_per_cell))
        + " reads",
        ha="left",
        fontsize=18,
        wrap=True,
    )

    plt.text(
        x=0.6 * x_max,
        y=saturation_pct * model_fit[0],
        s=str(int(saturation_pct * 100))
        + "% = "
        + str(int(mean_reads_per_cell))
        + " reads",
        ha="left",
        fontsize=18,
        wrap=True,
    )

    # Add miscellaneous (labels, legends, ...)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    lg = plt.legend(title=r"$R^2$ = " + "{0:.2f}".format(model_fit_r_sq), fontsize=18)
    title = lg.get_title()
    title.set_fontsize(18)

    if assay_type == "ATAC":
        plt.ylabel("Median Unique Fragments per Cell", fontsize=22)
        plt.xlabel("Mean Reads per Cell", fontsize=22)
        plt.ylim([0, model_fit[0]])
        plt.xlim([0, x_max])
    elif assay_type == "RNA":
        plt.ylabel("Saturation", fontsize=22)
        plt.xlabel("Mean Reads per Cell", fontsize=22)
        plt.ylim([0, 1.05])
        plt.xlim([0, x_max])

    plt.savefig(
        Path(args.output) / f"{project_name}_complexity.png",
        dpi=350,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.close()
    return reads_needed, saturation_pct, mean_reads_per_cell


x_data, y_data = prepare_data(
    complexity_info_path=complexity_info_path, assay_type=args.type
)
best_model_fit, best_model_fit_params, r_sq = fit_model(
    model=MM, x_data=x_data, y_data=y_data
)

reads_needed, saturation_pct, mean_reads_per_cell = plot_saturation_curve(
    model=MM,
    model_fit=best_model_fit,
    model_fit_params=best_model_fit_params,
    model_fit_r_sq=r_sq,
    x_data=x_data,
    assay_type=args.type,
    output_path=output_path,
)


# Initialize summary
reads_summary = []


# Update summary
reads_needed["num_cells"] = num_cells
reads_needed["saturation"] = int(saturation_pct * 100)
reads_needed["med_uniq_frag_per_cell"] = saturation_pct * best_model_fit[0]
reads_needed["mean_read_per_cell"] = int(mean_reads_per_cell)
reads_summary.append(reads_needed)

# Write summary of results
reads_summary = pd.DataFrame(reads_summary)
# from mean reads/cell to total reads
for perc in percentages:
    reads_summary[perc] = reads_summary[perc] * reads_summary["num_cells"]

reads_summary["now"] = reads_summary["now"] * reads_summary["num_cells"]
# Calculate additional reads needed
for perc in percentages:
    reads_summary[f"extra_needed_for_{perc * 100:.2f}"] = (
        reads_summary[perc] - reads_summary["now"]
    )
    reads_summary.drop(perc, axis=1, inplace=True)


# Clean up
reads_summary = reads_summary.astype("int")
reads_summary = reads_summary.sort_values(
    by=f"extra_needed_for_{min(percentages) * 100:.2f}", ascending=False
)
# Save summary
if args.type == "ATAC":
    reads_summary.to_csv(Path(args.output) / f"{project_name}_complexity.tsv", sep="\t")
elif args.type == "RNA":
    reads_summary.to_csv(Path(args.output) / f"{project_name}_saturation.tsv", sep="\t")

print(f"Done.")
