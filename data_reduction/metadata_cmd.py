from os import path
from pyDANDIA import  metadata
import argparse
import matplotlib.pyplot as plt
import numpy as np

def plot_vphas_cmd(args):
    default_marker_colour = '#8c6931'
    meta = metadata.MetaData()
    meta.load_all_metadata(args.red_dir, 'pyDANDIA_metadata.fits')

    # Select all stars with valid VPHAS+ data
    gdx = np.where(meta.star_catalog[1]['gmag'] > 0.0)[0]
    rdx = np.where(meta.star_catalog[1]['rmag'] > 0.0)[0]
    idx = np.where(meta.star_catalog[1]['imag'] > 0.0)[0]
    stars = set(gdx).intersection(set(rdx))
    stars = list(stars.intersection(set(idx)))

    # Plot CMD
    fig = plt.figure(1, (10, 10))
    ax = plt.subplot(111)
    plt.rcParams.update({'font.size': 25})
    gi = meta.star_catalog[1]['gmag'] - meta.star_catalog[1]['imag']
    plt.scatter(gi[stars], meta.star_catalog[1]['imag'][stars],
                c=default_marker_colour, marker='*', s=1)
    plt.xlabel('SDSS (g-i) [mag]')
    plt.ylabel('SDSS i [mag]')
    plt.title(args.title)
    (xmin,xmax,ymin,ymax) = plt.axis()
    plt.axis([xmin,xmax,ymax,ymin])
    plt.grid()

    plt.savefig(path.join(args.red_dir, 'vphas_CMD.png'))


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('red_dir', help='Path to reduction directory')
    parser.add_argument('title', help='Title for plot')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    plot_vphas_cmd(args)