#!/usr/bin/env python

"""
calculate_saturation.py: 
    Infer saturation of 10x scATAC/scRNA sample

    Requirements:
        - Python 3
        - Packages: pandas, numpy, scipy, json, matplotlib, re, os, uncertainties, argparse

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

__author__ = "Jasper Janssens"
__contributors__ = "Swan Floc’Hlay"
__version__ = "0.2.0"
__contact__ = "jasper.janssens@kuleuven.vib.be"

import re
import os
import sys
import json
import pandas as pd
from scipy.optimize import curve_fit
from uncertainties import ufloat
from pathlib import Path
import matplotlib.pylab as plt
import numpy as np

# get arguments
import argparse

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

# Open complexity file
with open(complexity_info_path) as complexity_info_fh:
    complexity_info = json.load(complexity_info_fh)
complexity_info_df = pd.DataFrame(complexity_info)

if args.type == "RNA":
    subsampled_columns = [
        x for x in complexity_info_df.columns if "multi_raw_rpc_" in x
    ]
    saturation_data = complexity_info_df[subsampled_columns].copy()
    saturation_data = pd.DataFrame(saturation_data.max())
    saturation_data.index = [
        np.float(
            re.sub("multi_raw_rpc_", "", re.sub("_subsampled_duplication_frac", "", x))
        )
        for x in subsampled_columns
    ]
    saturation_data = saturation_data.loc[saturation_data[0] != 0].copy()
    saturation_data = saturation_data.sort_values(by=0)
    saturation_data = saturation_data.reset_index()
    saturation_data.columns = [0, 1]

# initialize summary
reads_summary = []

# define Michaelis-Menten Kinetics
def MM(x, Vmax, Km):
    if Vmax > 0 and Km > 0:
        y = (Vmax * x) / (Km + x)
    else:
        y = 1e10
    return y


# get x/y data
if args.type == "ATAC":
    # get x data
    x_data = np.array(complexity_info_df["total_depth"])
    # get y data
    y_data = np.array(complexity_info_df["unique"])
elif args.type == "RNA":
    # get x data
    x_data = np.array(saturation_data[0])
    # get y data
    y_data = np.array(saturation_data[1])

######################################################################
# fit function
best_fit_ab, covar = curve_fit(MM, x_data, y_data)
sigma_ab = np.sqrt(np.diagonal(covar))

# get parameter value
a = ufloat(best_fit_ab[0], sigma_ab[0])
b = ufloat(best_fit_ab[1], sigma_ab[1])
text_res = "Best fit parameters:\na = {}\nb = {}".format(a, b)

# get rsquared value
residuals = y_data - MM(x_data, *(best_fit_ab))
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
r_sq = 1 - (ss_res / ss_tot)


plt.figure(figsize=(10, 7))
######################################################################
# plotting the model

if args.type == "ATAC":
    x_fit = np.linspace(0, int(np.max(x_data) * 100), 100000)
elif args.type == "RNA":
    x_fit = np.linspace(0, int(np.max(x_data) * 2000), 100000)
y_fit = MM(x_fit, *(best_fit_ab))

# get max value possible and highest we want to plot
y_max = best_fit_ab[0]
y_val = y_max * 0.9

# find x value at which fitted line surpasses y_val
x_coef = np.where(y_fit >= y_val)[0][0]
x_fit = x_fit[0:x_coef]
y_fit = y_fit[0:x_coef]
y_max, x_max = np.max(y_fit), np.max(x_fit)

# plot model
plt.plot(x_fit, MM(x_fit, *best_fit_ab), label="a = " + str(a), c="k", linewidth=3)
plt.plot(x_fit, MM(x_fit, *best_fit_ab), label="b = " + str(b), c="k", linewidth=3)

bound_upper = MM(x_fit, *(best_fit_ab + sigma_ab))
bound_lower = MM(x_fit, *(best_fit_ab - sigma_ab))

# plotting the confidence intervals
plt.fill_between(x_fit, bound_lower, bound_upper, color="gray", alpha=0.4)
######################################################################

# plot raw data
plt.scatter(x=x_data, y=y_data, c="red")

######################################################################
# plot different levels


def drawline(
    perc, coef, y_fit, x_fit, assay_type, x_max=None, color="k", linestyle="--"
):

    # get max value possible
    y_max = best_fit_ab[0]
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


# loop over different percentages + save results in dataframe
reads_needed = {}
for perc in percentages:
    x_coef = drawline(perc, best_fit_ab, y_fit, x_fit, assay_type=args.type)
    # update summary
    reads_needed[perc] = x_coef
    plt.text(
        x=0.6 * x_max,
        y=perc * best_fit_ab[0],
        s=str(int(perc * 100)) + "% = " + str(int(x_coef)) + " reads",
        ha="left",
        fontsize=18,
        wrap=True,
    )

perc = np.max(y_data) / best_fit_ab[0]
x_coef = drawline(
    perc=perc,
    coef=best_fit_ab,
    y_fit=y_data,
    x_fit=x_data,
    assay_type=args.type,
    x_max=x_max,
    color="r",
)
# update summary
reads_needed["now"] = x_coef
plt.text(
    x=0.6 * x_max,
    y=perc * best_fit_ab[0],
    s=str(int(perc * 100)) + "% = " + str(int(x_coef)) + " reads",
    ha="left",
    fontsize=18,
    wrap=True,
)

plt.text(
    x=0.6 * x_max,
    y=perc * best_fit_ab[0],
    s=str(int(perc * 100)) + "% = " + str(int(x_coef)) + " reads",
    ha="left",
    fontsize=18,
    wrap=True,
)

# update summary
reads_needed["num_cells"] = num_cells
reads_needed["saturation"] = int(perc * 100)
reads_needed["med_uniq_frag_per_cell"] = perc * best_fit_ab[0]
reads_needed["mean_read_per_cell"] = int(x_coef)
reads_summary.append(reads_needed)
######################################################################

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
lg = plt.legend(title=r"$R^2$ = " + "{0:.2f}".format(r_sq), fontsize=18)
title = lg.get_title()
title.set_fontsize(18)
if args.type == "ATAC":
    plt.ylabel("Median Unique Fragments per Cell", fontsize=22)
    plt.xlabel("Mean Reads per Cell", fontsize=22)
    plt.ylim([0, best_fit_ab[0]])
    plt.xlim([0, x_max])
elif args.type == "RNA":
    plt.ylabel("Saturation", fontsize=22)
    plt.xlabel("Mean Reads per Cell", fontsize=22)
    plt.ylim([0, 1.05])
    plt.xlim([0, x_max])

if args.type == "ATAC":
    complexity_info_path_stripping = re.sub(
        "/SC_ATAC_COUNTER_.*", "", str(complexity_info_path)
    )
    complexity_info_path_stripping = re.split("/", complexity_info_path_stripping)
    project_name = complexity_info_path_stripping[-1]

    plt.savefig(
        args.output + "/" + project_name + "_complexity.png",
        dpi=350,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.close()
elif args.type == "RNA":
    project_name = "RNA"  # needs to be updated
    plt.savefig(
        args.output + "/" + project_name + "_saturation.png",
        dpi=350,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.close()

# write summary of results
reads_summary = pd.DataFrame(reads_summary)
# from mean reads/cell to total reads
for perc in percentages:
    reads_summary[perc] = reads_summary[perc] * reads_summary["num_cells"]

reads_summary["now"] = reads_summary["now"] * reads_summary["num_cells"]
# calculate additional reads needed
for perc in percentages:
    reads_summary[f"extra_needed_for_{perc * 100:.2f}"] = (
        reads_summary[perc] - reads_summary["now"]
    )
    reads_summary.drop(perc, axis=1, inplace=True)


# clean up
reads_summary = reads_summary.astype("int")
reads_summary = reads_summary.sort_values(
    by=f"extra_needed_for_{min(percentages) * 100:.2f}", ascending=False
)
# save summary
if args.type == "ATAC":
    reads_summary.to_csv(args.output + "/" + project_name + "_complexity.tsv", sep="\t")
elif args.type == "RNA":
    reads_summary.to_csv(args.output + "/" + project_name + "_saturation.tsv", sep="\t")
