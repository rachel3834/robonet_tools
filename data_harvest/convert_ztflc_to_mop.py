from os import path
import argparse
import numpy as np
from astropy.table import Table, Column
import csv
import pandas as pd

def convert_lc(args):
    """
    Function to convert ZTF lightcurve data from the CSV file format they provide to the CSV format
    required for MOP upload.
    """

    # Load ZTF file format
    if not path.isfile(args.input_file):
        raise IOError('Cannot find input ZTF lightcurve file at ' + args.input_file)

    ztf_data = load_ztf_lc_csv(args.input_file)

    # Output lightcurves, one per filter, in MOP format
    output_to_mop_format(ztf_data, args.output_file)

def output_to_mop_format(lc_data, output_file):
    """
    Function to output Table format lightcurve data, potentially in multiple filters,
    to CSV format files suitable for upload to the MOP system.  If data from multiple
    filters is included, multiple lightcurve files will be created.
    """

    filter_set = np.unique(lc_data['filter'].data)

    file_root_name = (path.basename(output_file)).split('.')[0]
    output_file_root = path.join(path.dirname(output_file), file_root_name)

    for f in filter_set:
        file_path = output_file_root + '_' + str(f) + '.csv'
        idx = np.where(lc_data['filter'] == f)[0]

        with open(file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['time', 'filter', 'magnitude', 'error'])
            for i in idx:
                writer.writerow([
                    lc_data['JD'][i],
                    lc_data['filter'][i],
                    lc_data['mag'][i],
                    lc_data['mag_error'][i]
                ])

        print('Output ' + f + '-band lightcurve to ' + file_path)

def load_ztf_lc_csv(input_file):
    """
    Function to load a ZTF lightcurve from a CSV-format input file downloaded from the FINK portal
    https://fink-portal.org/api
    """

    # Mapping of filter indices, taken from
    # https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/ztf_release_notes_latest
    filter_map = {
        1: 'ZTF-g',
        2: 'ZTF-r',
        3: 'ZTF-i'
    }

    with open(input_file) as f:
        raw_data = pd.read_csv(f)

    # Remove invalid measurements
    valid = np.where(np.isnan(raw_data['i:magpsf']) == False)[0]

    lc_data = Table([
        Column(name='JD', data=raw_data['i:jd'][valid]),
        Column(name='mag', data=raw_data['i:magpsf'][valid]),
        Column(name='mag_error', data=raw_data['i:sigmapsf'][valid]),
        Column(name='filter', data=raw_data['i:fid'][valid].apply(lambda x: filter_map[x]))
    ])

    return lc_data

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to input ZTF lightcurve CSV file')
    parser.add_argument('output_file', help='Path to output MOP-format CSV file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    convert_lc(args)
