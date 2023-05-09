from os import path
import argparse
from pyDANDIA import crossmatch
import numpy as np

def run_statistics(args):

    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.crossmatch_file,log=None)

    calc_image_totals(xmatch)

def calc_image_totals(xmatch):
    """Function to calculate the total number of images across all datasets"""
    filter_list = ['gp', 'rp', 'ip']

    filter_count = {'gp': 0, 'rp': 0, 'ip': 0}
    for i,dset in enumerate(xmatch.datasets['dataset_code']):
        idx = np.where(xmatch.images['dataset_code'] == dset)[0]
        print(dset+' contains '+str(len(idx)))
        f = xmatch.datasets['dataset_filter'][i]
        filter_count[f] += len(idx)

    for f in filter_list:
        print('Total '+f+' '+str(filter_count[f]))

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('crossmatch_file', type=str,
                    help='Path to crossmatch file')

    return parser.parse_args()

if __name__ == '__main__':
  args = get_args()
  run_statistics(args)
