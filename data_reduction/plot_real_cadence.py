from os import path
import argparse
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def plot_cadence_data(args):

    data = load_cadence_data(args)

    (fig, axlist) = plt.subplots(3)
    cmap = plt.get_cmap('PuBuGn')
    norm = mpl.colors.Normalize(vmin=0.0, vmax=25.0)

    for ip, filter in enumerate(['g', 'r', 'i']):
        dataset = data[filter]
        img = axlist[ip].imshow(np.array(dataset['cadence']), cmap=cmap, norm=norm)
        inityticks = np.array(axlist[ip].get_yticks())
        yincr = (inityticks.max() - inityticks.min())/3.0
        yticks = np.arange(inityticks.min()+yincr/2.0, inityticks.max(), yincr)
        yticklabels = [dataset['rows'][x][0] + ' to ' + dataset['rows'][x][1] for x in [0,1,2]]
        xticks = np.arange(0,20,1, dtype='int')
        xticklabels = [str(x+1) for x in xticks]
        axlist[ip].set_yticks(yticks)
        axlist[ip].set_yticklabels(yticklabels)
        axlist[ip].set_xticks(xticks)
        axlist[ip].set_xticklabels(xticklabels)
        axlist[ip].set_title('Cadence in ' + filter + '-band')
        if ip==2:
            axlist[ip].set_xlabel('Field number')
    plt.subplots_adjust(left=0.32, right=0.85, top=0.98, bottom=0.1, hspace=-0.4)

    cbar_ax = fig.add_axes([0.87, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(img, cax=cbar_ax)
    cbar.ax.set_ylabel('Median observation interval [hrs]')

    #plt.show()
    plt.savefig(path.join(args.output_dir, 'realized_cadence.png'), bbox_inches='tight')

def load_cadence_data(args):

    with open(args.input_file, 'r') as f:
        data = json.load(f)

    return data

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to input data file')
    parser.add_argument('output_dir', help='Path to output directory')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    plot_cadence_data(args)