from os import path
import argparse
import csv
import pandas as pd
from astropy.table import Table, Column

def convert_csv_file(args):
    """
    Function to convert the CSV of photometry downloaded from the BHTOM system
    to a set of lightcurve files in the format that can be uploaded to MOP
    """

    # Read the BHTOM-format CSV photometry file
    bhtom_data = load_bhtom_photometry(args)

    # Convert the data columns to the MOP-standard format
    mop_data = bhtom_to_mop(bhtom_data)

    # Output data in MOP-compatible format
    output_mop_csv(args, mop_data)

def output_mop_csv(args, mop_data):
    """
    Function to output a Table of photometry data to the MOP-format CSV file
    """

    with open(args.output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerow(['time','filter','magnitude','error'])
        for i in range(0, len(mop_data), 1):
            csvwriter.writerow(mop_data[i])

def bhtom_to_mop(bhtom_data):
    """
    Function to repackage the BHTOM photometry data into the columns accepted by MOP's
    photometry processor
    """

    label = ['BH-'+bhtom_data['facility'][i]+'-'+bhtom_data['filter'][i]+'-'+bhtom_data['observer'][i]
             for i in range(0,len(bhtom_data),1)]

    mop_data = Table([
        Column(name='time', data=bhtom_data['MJD']),
        Column(name='filter', data=label),
        Column(name='magnitude', data=bhtom_data['mag']),
        Column(name='error', data=bhtom_data['mag_error']),
    ])

    return mop_data

def load_bhtom_photometry(args):
    """
    Function to load a CSV file of photometry in the BHTOM-standard output format
    This is expected to be a semi-colon-separated CSV file with columns:
    MJD;Magnitude;Error;Facility;Filter;Observer
    """
    with open(args.input_file) as f:
        raw_data = pd.read_csv(f, delimiter=';')

    lc_data = Table([
        Column(name='MJD', data=raw_data['MJD']),
        Column(name='mag', data=raw_data['Magnitude']),
        Column(name='mag_error', data=raw_data['Error']),
        Column(name='facility', data=raw_data['Facility']),
        Column(name='filter', data=raw_data['Filter']),
        Column(name='observer', data=raw_data['Observer']),
    ])

    return lc_data

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="Path to input CSV file of photometry from BHTOM")
    parser.add_argument('output_file', help="Path to output lightcurve")
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    convert_csv_file(args)
