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

        # CATCH FOR DATA WITH BAD DATE-OBS
        # Due to a temporary operational issue at LCO sites, the DATE-OBS parameter recorded
        # in the image headers for some images on 2019-04-30 were set to 1970-01-01.  The DATE
        # parameter in the corresponding header appears to be correct but records the time of
        # file write rather than the exposure start time.  Using this to estimate the DATE-OBS would
        # be close but inexact, and I've decided to simply flag these frames as bad data instead.
        # Affected frames have a zero-length exposure, so no valid photometry is produced.
        # They can be identified because the calculated HJD will be prior to 2017-01-01 = JD 2457754.5,
        # which is the earliest possible date for the ROME survey.
        if hjd > 2457754.5:
            xmatch.images['hjd'][i] = hjd
        else:
            log.info('Caught bad DATE-OBS: '
             + xmatch.images['filename'][i] + ' DATE-OBS=' + xmatch.images['datetime'][i]
             + ' RA, Dec=' + str(xmatch.images['RA'][i]) + ', ' + str(xmatch.images['Dec'][i])
             + ' ExpTime=' + str(xmatch.images['exposure'][i]))
            xmatch.images['hjd'][i] = 0.0

    log.info('Verified HJD calculations')

    # Store the verified timestamps
    xmatch.save(crossmatch_file)
    log.info('Updated crossmatch table')

    # Verify the HJD timestamps in the timeseries photometry data.  This section can be switched off
    # for faster testing
    apply_to_phot = True
    if apply_to_phot:
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