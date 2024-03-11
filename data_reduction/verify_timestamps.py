from os import path
import argparse
from pyDANDIA import time_utils
from pyDANDIA import crossmatch
from pyDANDIA import field_photometry
from pyDANDIA import logs
from pyDANDIA import hd5_utils

def verify_field_timestamps(args):
    """
    Function to verify the Heliocentric Julian Date timestamps for a given field by
    re-calculating them from scratch.

    This function then updates the HJD timestamps stored in the crossmatch Image table
    and the photometry HDF5 files for each quadrant.
    """

    log = logs.start_stage_log(args.data_dir, 'verify_timestamps')

    # Load images table from the crossmatch file:
    xmatch = crossmatch.CrossMatchTable()
    crossmatch_file = path.join(args.data_dir, args.field_name + '_field_crossmatch.fits')
    xmatch.load(crossmatch_file, log=None)
    log.info('Loaded crossmatch table')

    # Re-calculate HJD timestamps for all images in the dataset
    # This is a loop because the datasets include data from multiple sites,
    # the geographic locations of which are looked up separately
    for i in range(0, len(xmatch.images), 1):
        tel_code = '-'.join(
            xmatch.images['dataset_code'][i].split('_')[1].split('-')[0:3]
        ) + 'a'
        hjd, ltt_helio = time_utils.calc_hjd(
            xmatch.images['datetime'][i],
            xmatch.images['RA'][i], xmatch.images['Dec'][i],
            tel_code,
            xmatch.images['exposure'][i],
            debug=False
        )

        xmatch.images['hjd'][i] = hjd
    log.info('Verified HJD calculations')

    # Store the verified timestamps
    xmatch.save(crossmatch_file)
    log.info('Updated crossmatch table')

    # Verify the HJD timestamps in the timeseries photometry data
    for qid in range(1, 5, 1):     
        phot_file = path.join(args.data_dir,
                              args.field_name + '_quad' + str(qid) + '_photometry.hdf5')

        if not path.isfile(phot_file):
            raise IOError('Cannot find input timeseries photometry file ' + phot_file)

        quad_phot = hd5_utils.read_phot_from_hd5_file(phot_file,
                                                      return_type='array')
        quad_phot[:,:,0] = xmatch.images['hjd']

        # Output tables to quadrant HDF5 file
        params = {'crossmatch_file': crossmatch_file, 'field_name': args.field_name}
        field_photometry.output_quad_photometry(params, xmatch, quad_phot, qid, log)

    logs.close_log(log)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Path to data directory')
    parser.add_argument('field_name', help='Field name file prefix')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    verify_field_timestamps(args)