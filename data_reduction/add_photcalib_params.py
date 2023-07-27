from sys import argv
from os import path
import numpy as np
from pyDANDIA import metadata

def add_params_metadata(red_dir):
    """Function to add the photometric covarience parameters from the
    phot_calib log to the corresponding metadata table"""

    reduction_metadata = metadata.MetaData()
    reduction_metadata.load_all_metadata(red_dir, 'pyDANDIA_metadata.fits')

    covar_fit = parse_phot_calib_log(red_dir)

    reduction_metadata.phot_calib[1]['c0'] = covar_fit[0,0]
    reduction_metadata.phot_calib[1]['c1'] = covar_fit[0,1]
    reduction_metadata.phot_calib[1]['c2'] = covar_fit[1,0]
    reduction_metadata.phot_calib[1]['c3'] = covar_fit[1,1]

    reduction_metadata.save_a_layer_to_file(red_dir,
                                            'pyDANDIA_metadata.fits',
                                            'phot_calib')
    print('Added photometric calibration parameters to metadata for '+path.basename(red_dir))
    
def parse_phot_calib_log(red_dir):

    log_file = path.join(red_dir, 'phot_calib.log')
    if not path.isfile(log_file):
        raise IOError('ERROR: Cannot find photometric log file '+log_file)

    covar_fit = np.zeros((2,2))
    logdata = open(log_file,'r').readlines()
    for i,line in enumerate(logdata):
        if 'Fit covarience' in line:
            entries = line.replace('array([[','').replace(']',' ').replace(',','').split()
            covar_fit[0,0] = float(entries[-2])
            covar_fit[0,1] = float(entries[-1])
            break
    entries = logdata[i+1].replace(']])\n','').replace(',','').replace('[','').split()
    covar_fit[1,0] = float(entries[-2])
    covar_fit[1,1] = float(entries[-1])

    return covar_fit

if __name__ == '__main__':
    if len(argv) == 1:
        red_dir = input('Please enter the path to the reduction directory: ')
    else:
        red_dir = argv[1]
    add_params_metadata(red_dir)
