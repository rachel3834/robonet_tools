from os import path
from sys import argv
from pyDANDIA import metadata
import numpy as np

def set_psf_dimensions():

    params = get_args()

    reduction_metadata = metadata.MetaData()
    reduction_metadata.load_all_metadata(metadata_directory=params['red_dir'],
                                         metadata_name='pyDANDIA_metadata.fits')

    reduction_metadata.calc_psf_radii()

    reduction_metadata.save_updated_metadata(params['red_dir'],'pyDANDIA_metadata.fits')

def get_args():

    params = {}
    if len(argv) == 1:
        params['red_dir'] = input('Please enter the path to the reduction directory: ')
    else:
        params['red_dir'] = argv[1]

    return params


if __name__ == '__main__':
    set_psf_dimensions()
