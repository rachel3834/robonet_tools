from os import path, rename, getcwd
from sys import argv, executable
import numpy as np
from pyDANDIA import calibrate_photometry
from pyDANDIA import pipeline_setup
from pyDANDIA import logs
import subprocess

def rerun_phot_calibration(datasets, software_dir):
    """Function to re-run the photometric calibration for a dataset,
    based on the parameters of extracted from the log of the previous run of
    calibrate_photometry"""

    for data_dir in datasets:
        print('-> Calibrating '+path.basename(data_dir))

        # Use the data_dir to create a pipeline_setup object and
        # params dictionary with red_dir, metadata, log_dir
        params = {'red_dir': data_dir, 'log_dir': data_dir,
                  'metadata': 'pyDANDIA_metadata.fits',
                  'use_gaia_phot': False, 'set_phot_calib': False,
                  'software_dir': software_dir,
                  'cat_mags_min': None, 'a0': None, 'a1': None}
        setup = pipeline_setup.pipeline_setup(params)

        # Parse the log from the previous run of calibration_photometry and
        # extract the values for parameters det_mags_max, det_mags_min,
        # cat_mags_max and cat_merr_max
        (status, params) = parse_calib_phot_log(params)
        backup_phot_calib_log(setup)

        if status:
            command = path.join(software_dir,'calibrate_photometry.py')
            args = [executable, command, data_dir,
                        'pyDANDIA_metadata.fits', data_dir]
            for key in ['det_mags_max', 'det_mags_min', 'cat_merr_max', \
                        'cat_mags_max', 'cat_mags_min']:
                if key in params.keys() and 'none' not in str(params[key]).lower():
                    args.append(key+'='+str(params[key]))

            # Calling out as a separate process so that the logs are properly
            # handled and not concatenated into a single file when multiple
            # datasets are calibrated. 
            #(status, report) = calibrate_photometry.calibrate_photometry_catalog(setup, **params)
            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            p.wait()

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

def backup_phot_calib_log(setup):

    phot_file_name = path.join(setup.red_dir,'phot_calib.log')
    idx = 1
    bkup_file_name = path.join(setup.red_dir,'phot_calib.log.'+str(idx))
    if path.isfile(phot_file_name):
        while path.isfile(bkup_file_name):
            idx +=1
            bkup_file_name = path.join(setup.red_dir,'phot_calib.log.'+str(idx))
        rename(phot_file_name, bkup_file_name)

def get_args():
    datasets = []

    if len(argv) > 1:
        option = argv[1]
        software_dir = argv[2]
    else:
        option = input('Please enter the path to a reduction directory, or path to list of directories with @ prefix: ')
        software_dir = input('Please enter the path to the pyDANDIA software directory: ')

    if '@' in option[0:1]:
        if not path.isfile(option[1:]):
            raise IOError('Cannot find input list of datasets at '+option[1:])

        print('Reading list of datasets from '+option[1:])
        file_lines = open(option[1:]).readlines()
        for line in file_lines:
            if len(line.replace('\n','')) > 0:
                entries = line.replace('\n','').split()
                datasets.append(entries[0])
        print('Found '+str(len(datasets))+' datasets to process')

    else:
        datasets.append(option)
        print('Received dataset '+datasets[0]+' to process')

    return datasets, software_dir

if __name__ == '__main__':
    (datasets, software_dir) = get_args()
    rerun_phot_calibration(datasets, software_dir)
