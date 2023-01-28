from os import path
from sys import argv
import numpy as np
from pyDANDIA import calibrate_photometry
from pyDANDIA import pipeline_setup
from pyDANDIA import logs

def rerun_phot_calibration(datasets):
    """Function to re-run the photometric calibration for a dataset,
    based on the parameters of extracted from the log of the previous run of
    calibrate_photometry"""

    for data_dir in datasets:

        # Use the data_dir to create a pipeline_setup object and
        # params dictionary with red_dir, metadata, log_dir
        params = {'red_dir': data_dir, 'log_dir': data_dir,
                  'metadata': 'pyDANDIA_metadata.fits',
                  'use_gaia_phot': False, 'set_phot_calib': False,
                  'cat_mags_min': None, 'a0': None, 'a1': None}
        setup = pipeline_setup.pipeline_setup(params)

        # Parse the log from the previous run of calibration_photometry and
        # extract the values for parameters det_mags_max, det_mags_min,
        # cat_mags_max and cat_merr_max
        (status, params) = parse_calib_phot_log(params)

        if status:
            (status, report) = calibrate_photometry.calibrate_photometry_catalog(setup, **params)

def parse_calib_phot_log(params):
    log_file = path.join(params['red_dir'], 'phot_calib.log')
    par_list = ['det_mags_max', 'det_mags_min', 'cat_mags_max', 'cat_merr_max']

    if not path.isfile(log_file):
        print('WARNING: cannot find calibrate_photometry log for '+params['red_dir']+', skipping dataset')
        status = False

    else:
        status = True
        file_lines = open(log_file).readlines()
        for line in file_lines:
            for par in par_list:
                # Note colon is required to avoid getting multiple entries in file
                if par+':' in line:
                    entries = line.replace('\n','').split()
                    params[par] = float(entries[2])

    return status, params

def get_list_of_datasets():
    datasets = []

    if len(argv) > 1:
        option = argv[1]
    else:
        option = input('Please enter the path to a reduction directory, or path to list of directories with @ prefix: ')

    if '@' in option[0:1]:
        if not path.isfile(option[1:]):
            raise IOError('Cannot find input list of datasets at '+option[1:])

        file_lines = open(option[1:]).readlines()
        for line in file_lines:
            if len(line.replace('\n','')) > 0:
                entries = line.replace('\n','').split()
                datasets.append(entries[0])

    else:
        datasets.append(option)

    return datasets

if __name__ == '__main__':
    datasets = get_list_of_datasets()
    rerun_phot_calibration(datasets)
