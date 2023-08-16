import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import os
from sys import argv
plt.rc('font', family='DejaVu Sans')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')


def map_MCMC_chains(input_file):
    """Function to plot maps of the (logs .vs. logq) parameter space mapped
    out by a Markov Chain Monte Carlo (MCMC) algorithm search.

    The first plot is a 2D histogram of the number of MCMC chains
    for each pixel in (logs .vs. logq) space, as a proxy for a map of the
    population of solutions calculated during a
    search of parameter space.

    The second plot is a 2D histogram of the log likelihood of each pixel
    in the parameter space.
    """

    # Load the output data from MCMC chains
    results = np.loadtxt(input_file)

    # Plot ln likelihood data for all chains as a 2D function of binary mass ratio and projected separation
    fontsize = 16
    ticksize = 14
    fig = plt.figure(1, (10, 5))
    plt.scatter(results[:, :, 4].ravel(), results[:, :, 5].ravel(),
                c=results[:, :, -1].ravel(), alpha=0.25)
    plt.xlabel('$log_{10}$(Binary projected separation)', fontsize=fontsize)
    plt.ylabel('$log_{10}$(Binary mass ratio)', fontsize=fontsize)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    cbar = plt.colorbar(label="-ln likelihood")
    cbar.ax.tick_params(labelsize=ticksize)

    plt.savefig('MCMC_likelihood_map.png')

    plt.close(1)


if __name__ == '__main__':

    if len(argv) == 1:
        input_file = raw_input('Please enter the path to the DE population output: ')

    else:
        input_file = argv[1]

    map_MCMC_chains(input_file)
