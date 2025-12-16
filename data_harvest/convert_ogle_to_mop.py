from os import path
import argparse
import numpy as np
from astropy.table import Table, Column
import csv

def convert_lc(args):
    """
    Function to convert ZTF lightcurve data from the CSV file format they provide to the CSV format
    required for MOP upload.
    """

    # Load ZTF file format
    if not path.isfile(args.input_file):
        raise IOError('Cannot find input ZTF lightcurve file at ' + args.input_file)

    ogle_data = load_ogle_lc(args)

    # Output lightcurves, one per filter, in MOP format
    output_to_mop_format(ogle_data, args.output_file)

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
                    lc_data['HJD'][i],
                    lc_data['filter'][i],
                    lc_data['mag'][i],
                    lc_data['mag_error'][i]
                ])

        print('Output ' + f + '-band lightcurve to ' + file_path)

def load_ogle_lc(args):
    """
    Function to load an OGLE lightcurve from a dat-format input file downloaded from the project website
    """

    with open(args.input_file) as f:
        raw_data = np.loadtxt(f)

    lc_data = Table([
        Column(name='HJD', data=raw_data[:,0]),
        Column(name='mag', data=raw_data[:,1]),
        Column(name='mag_error', data=raw_data[:,2]),
        Column(name='filter', data=[args.filter]*len(raw_data))
    ])

    return lc_data

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to input OGLE lightcurve dat file')
    parser.add_argument('filter', help='Name of OGLE filter for label')
    parser.add_argument('output_file', help='Path to output MOP-format CSV file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    convert_lc(args)
