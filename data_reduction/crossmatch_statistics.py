from os import path
import argparse
from pyDANDIA import crossmatch
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

def run_statistics(args):

    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.crossmatch_file,log=None)

    filter_list = ['gp', 'rp', 'ip']

    calc_image_totals(xmatch, filter_list)
    calc_star_totals(xmatch)
    calc_cadence(xmatch, filter_list)

def calc_image_totals(xmatch, filter_list):
    """Function to calculate the total number of images across all datasets"""

    filter_count = {'gp': 0, 'rp': 0, 'ip': 0}
    for i,dset in enumerate(xmatch.datasets['dataset_code']):
        idx = np.where(xmatch.images['dataset_code'] == dset)[0]
        print(dset+' contains '+str(len(idx)))
        f = xmatch.datasets['dataset_filter'][i]
        filter_count[f] += len(idx)

    for f in filter_list:
        print('Total '+f+' '+str(filter_count[f]))

def calc_star_totals(xmatch):
    print('N stars = '+str(len(xmatch.stars)))

def calc_cadence(xmatch, filter_list):
    survey_seasons = [(datetime.strptime('2017-01-01',"%Y-%m-%d"),
                        datetime.strptime('2017-10-30',"%Y-%m-%d")),
                        (datetime.strptime('2018-01-01',"%Y-%m-%d"),
                        datetime.strptime('2018-10-30',"%Y-%m-%d")),
                        (datetime.strptime('2019-01-01',"%Y-%m-%d"),
                        datetime.strptime('2019-10-30',"%Y-%m-%d"))]

    image_dates = []
    for ts in xmatch.images['datetime']:
        image_dates.append(datetime.fromisoformat(ts)) # "%Y-%m-%dT%H:%M:%D"
    image_dates = np.array(image_dates)

    (fig, ax) = plt.subplots(1, 1, figsize=(10, 10))

    for f in filter_list:
        idx1 = np.where(xmatch.images['filter'] == f)[0]
        idx2 = np.where(xmatch.images['hjd'] > 0.0)[0]
        cadences = np.array([])
        for (season_start, season_end) in survey_seasons:
            idx3 = np.where(image_dates > season_start)[0]
            idx4 = np.where(image_dates <= season_end)[0]
            idx = set(idx1).intersection(set(idx2))
            idx = idx.intersection(set(idx3))
            idx = list(idx.intersection(set(idx4)))

            hjds = xmatch.images['hjd'][idx]
            hjds = np.sort(hjds)
            dhjds = hjds[1:] - hjds[0:-1]
            cadence = np.median(dhjds) * 24.0
            sigma = dhjds.std() * 24.0
            cadences = np.concatenate((cadences, dhjds*24.0))

            print('Median cadence in '+f+' for season '
                    +season_start.strftime("%Y-%m-%d")
                    +' to '+season_end.strftime("%Y-%m-%d")+' = '+str(round(cadence,1))+'hrs')

        ax.hist(cadences, bins='auto', label=f)

    ax.set_xlabel('Cadence [hrs]')
    ax.set_ylabel('Frequency')
    plt.legend()
    params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large'}
    plt.rcParams.update(params)
    plt.xlim(0,100)
    dirpath = path.dirname(args.crossmatch_file)
    plt.tight_layout()
    plt.savefig(path.join(dirpath, 'logs', 'cadence_histogram.png'), )
    plt.close(1)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('crossmatch_file', type=str,
                    help='Path to crossmatch file')

    return parser.parse_args()

if __name__ == '__main__':
  args = get_args()
  run_statistics(args)
