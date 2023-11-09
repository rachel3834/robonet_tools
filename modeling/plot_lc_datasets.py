from os import path
import argparse
import matplotlib.pyplot as plt
import numpy as np
import plot_phase_curve
import json

def plot_lightcurve_datasets(config):
    dt = 2450000.0
    colours = ['b', 'g', 'm', 'r', 'k']
    markers = ['*', 'v', 'x']     # One per bandpass

    fig = plt.figure(1, (15, 10))
    plt.rcParams.update({'font.size': 18})

    for j,f in enumerate(['g', 'r', 'i']):
        for i,dataset in enumerate(config[f]):
            data = plot_phase_curve.load_data(dataset['file'])

            plt.errorbar(data[:, 0] - dt, data[:, 1], yerr=data[:, 2],
                 fmt=colours[i]+markers[j], label=dataset['label'])

    plt.xlabel('HJD - ' + str(dt))
    plt.ylabel('Magnitude')
    plt.legend()
    plt.grid('top-right')

    [xmin, xmax, ymin, ymax] = plt.axis()
    plt.axis([xmin-100, xmax+300, ymax, ymin])

    plt.annotate('SDSS-g', (7700, 19.8))
    plt.annotate('SDSS-r', (7700, 18.0))
    plt.annotate('SDSS-i', (7700, 15.5))
    plt.title(config['title'])

    plt.savefig(config['output_file'], bbox_inches='tight')

    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help="Path to plot configuration file")
    args = parser.parse_args()

    config = json.loads(open(args.config_file,'r').read())

    plot_lightcurve_datasets(config)