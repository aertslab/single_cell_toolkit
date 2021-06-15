#!/usr/bin/env python

### load libs
import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pylab as plt
import bisect
import random
from collections import Counter, Sequence
import pyranges as pr

### get arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help='path to fragment file', required=True)
parser.add_argument('-o','--output', help='path to output dir', default = "./")
parser.add_argument('-p','--percentages', metavar='percentages', type=str, help='Comma separated list of decimal percentages to predict', default='0.3,0.6,0.9')
parser.add_argument('-m','--minimum_frag_for_bc', metavar='integer', type=int, help='Minimum number of uniq fragments to select real cell barcodes', default=200)
parser.add_argument('-n','--n_chunk', metavar='integer', type=int, help='Number of sub-sampling to perform', default=10)
parser.add_argument('-w','--whitelist', metavar='string', type=str, help='Path to file with barcode whitelist', default=None)

# Parse arguments
args = parser.parse_args()
percentages = [float(x) for x in args.percentages.split(',')]

### initialise function and classes

# class of fragment data weighted by duplicate count for sampling
class WeightedPopulation(Sequence) :
    def __init__(self, population, weights):
        assert len(population) == len(weights) > 0
        self.population = population
        self.cumweights = []
        cumsum = 0 # compute cumulative weight
        for w in weights:
            cumsum += w   
            self.cumweights.append(cumsum)  
    def __len__(self):
        return self.cumweights[-1]
    def __getitem__(self, i):
        if not 0 <= i < len(self):
            raise IndexError(i)
        return self.population[bisect.bisect(self.cumweights, i)]


# fragment file loading function
def read_fragment(path_to_frag) :
    fragments_df = pr.read_bed(path_to_frag).df
    if 'Score' not in fragments_df or all(fragments_df['Score'] == '.'):
        raise ValueError('Duplicate information is missing')
    # Convert to category/float32 for memory efficiency
    fragments_df["Name"] = fragments_df["Name"].astype("category")
    fragments_df["Score"] = fragments_df["Score"].astype("int32")
    return(fragments_df)


def MM(x, Vmax, Km):
    """
    Define the Michaelis-Menten Kinetics model that will be used for the model fitting.
    """
    if Vmax > 0 and Km > 0:
        y = (Vmax * x) / (Km + x)
    else:
        y = 1e10
    return y


# sub-sampling function
def sub_sample_fragment(fragments_df,
                        min_uniq_frag = 200,
                        n_chunk = 10,
                        outfile = "sampling_stats.tab",
                        whitelist = None) :
    # sample on total size with weight
    pop = WeightedPopulation(list(fragments_df.index), fragments_df['Score'])
    lp = len(pop)
    sublist = random.sample(pop , lp)
    # init stats buket
    stat_buket = {'mean_frag_per_bc' : {'chunk 0' : 0} , # mean read per cell
        'median_uniq_frag_per_bc' : {'chunk 0' : 0} , # median uniq fragment per cell
        'total_frag_count' : {'chunk 0' : 0} , # total read count (all barcodes)
        'cell_barcode_count' : {'chunk 0' : 0} # number of barcodes with n_reads > min_uniq_frag
        }
    # Find real cell barcode
    chunk = 'chunk ' + str(n_chunk)
    pop_count = pd.DataFrame.from_dict(Counter(sublist), orient='index', columns = [chunk])
    pop_count = pop_count.join(fragments_df)
    tmp_mufpc = pop_count.astype({chunk : 'bool'}).groupby("Name")[chunk].sum()
    bc = list(tmp_mufpc[tmp_mufpc > min_uniq_frag].index)
    print("Keeping " + str(len(bc)) + " barcodes with sufficient fragment counts")
    # Filter barcodes with whitelist
    if whitelist is not None :
        whitelist_bc =  list(pd.read_csv(whitelist, names = ["valid_bc"])["valid_bc"])
        bc = [x for x in bc if x in whitelist_bc]
        print("Keeping " + str(len(bc)) + " barcodes in whitelist")
    # compute stats per chunk of increasing size
    for i in range(1,n_chunk + 1) :
        chunk = 'chunk ' + str(i)
        print('Processing ' + chunk)
        pop_count = [sublist[index] for index in range(0,round(i * lp/n_chunk))]
        pop_count = pd.DataFrame.from_dict(Counter(pop_count), orient='index', columns = [chunk])
        # Reattribute cell barcode
        pop_count = pop_count.join(fragments_df)
        tmp_mrpc = pop_count.groupby("Name")[chunk].sum()
        tmp_mufpc = pop_count.astype({chunk : 'bool'}).groupby("Name")[chunk].sum()
        # Count fragment irrespective of bc threshold
        stat_buket['total_frag_count'][chunk] = np.sum(tmp_mrpc)
        # select barcode and compute stats
        stat_buket['mean_frag_per_bc'][chunk] = np.mean(tmp_mrpc.loc[bc])
        stat_buket['median_uniq_frag_per_bc'][chunk] = np.median(tmp_mufpc.loc[bc])
        stat_buket['cell_barcode_count'][chunk] = len(bc)
    # Save data as tab file
    stat_buket = pd.DataFrame(stat_buket)
    stat_buket.to_csv(outfile, sep = "\t")
    return(stat_buket)


# MM-fit function
def fit_MM(stat_bucket, percentages = [0.3,0.6,0.9], path_to_fig = "./", x_axis = "total_frag_count", y_axis = "median_uniq_frag_per_bc"):
    # select x/y data fro MM fit from subsampling stats
    x_data = np.array(stat_bucket.loc[:,x_axis])
    y_data = np.array(stat_bucket.loc[:,y_axis])
    # fit to MM function
    best_fit_ab, covar = curve_fit(MM, x_data, y_data, bounds=(0, +np.inf))
    # expand fit space
    x_fit = np.linspace(0, int(np.max(x_data) * 10), num = 200)
    y_fit = MM(x_fit, *(best_fit_ab))
    # impute maximum saturation to plot as 95% of y_max
    y_val = best_fit_ab[0] * 0.95
    # subset x_fit space if bigger then y_val
    if y_val < max(y_fit) :
        x_coef = np.where(y_fit>=y_val)[0][0]
        x_fit = x_fit[0:x_coef]
        y_fit = y_fit[0:x_coef]
    # plot model
    plt.plot(x_fit, MM(x_fit, *best_fit_ab),label="fitted", c='black',linewidth=1)
    # plot raw data
    plt.scatter(x=x_data,y=y_data,c='crimson', s = 10)
    # mark curent saturation
    x_idx = np.where(y_fit>=max(y_data))[0][0]
    x_coef = x_fit[x_idx]
    y_coef = y_fit[x_idx]
    plt.plot([x_coef, x_coef], [0, y_coef], linestyle = '--', c = 'crimson')
    plt.plot([0, x_coef], [y_coef, y_coef], linestyle = '--', c = 'crimson')
    plt.text(x = x_fit[-1], y = y_coef, s = str(round(100 * max(y_data) / best_fit_ab[0])) + "% {:.2e}".format(round(x_coef)), c = "crimson", ha = 'right', va = 'bottom')
    # plot percentaged values
    for perc in percentages:
        # Find read count for percent saturation
        y_val = best_fit_ab[0]*perc
        # Find closest match in fit
        if max(y_fit) > y_val :
            x_idx = np.where(y_fit>=y_val)[0][0]
            x_coef = x_fit[x_idx]
            y_coef = y_fit[x_idx]
            # Draw vline
            plt.plot([x_coef, x_coef], [0, y_coef], linestyle = '--', c = 'grey')
            # Draw hline
            plt.plot([0, x_coef], [y_coef, y_coef], linestyle = '--', c = 'grey')
            # Plot imputed read count
            plt.text(x = x_fit[-1], y = y_coef, s = str(round(100*perc)) + "% {:.2e}".format(round(x_coef)), c = "grey", ha = 'right', va = 'bottom')
    # save figure
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.savefig(path_to_fig)
    plt.close()


### Run script for argparsed agruments
# load frag files
fragments_df = read_fragment(args.input)
# sub-sample
stat_buket = sub_sample_fragment(fragments_df,
                        min_uniq_frag = args.minimum_frag_for_bc,
                        n_chunk = args.n_chunk,
                        outfile = args.output + "sampling_stats.tab",
                        whitelist = args.whitelist)
# fit'n'plot for total count
fit_MM(stat_buket, percentages = percentages, path_to_fig = args.output + "saturation.png", x_axis = "total_frag_count", y_axis = "median_uniq_frag_per_bc")
