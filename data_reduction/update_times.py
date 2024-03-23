from os import path
import argparse
from astropy.io import fits
from astropy.time import Time, TimeDelta
import glob
import copy

# On 2024-03-20, the clocks on site at LSC-1m04 were early by 7107.053965sec according to
# LCO Tel-Ops.  This script is designed to amend the header information to the correct time.
# The affected images from this program are as follows:
AFFECTED_IMAGES = [
    'lsc1m004-fa03-20240320-0352-e91.fits',
    'lsc1m004-fa03-20240320-0353-e91.fits',
    'lsc1m004-fa03-20240320-0354-e91.fits',
]

# The clocks were measured to be early by 7107.053965sec
dt = TimeDelta(7107.053965, format='sec')

def run(args):

    # Make a list of images in the data directory
    image_list = glob.glob(path.join(args.data_dir, '*fits'))

    # If one of the affected images is present, edit the affected timestamps in its header
    # DATE-OBS, UTSTART, UTEND, MJD-OBS
    for image in image_list:
        if path.basename(image) in AFFECTED_IMAGES:

            hdu = fits.open(image)
            hdu2 = copy.deepcopy(hdu)

            dateobs = hdu[0].header['DATE-OBS']
            utstart = hdu[0].header['UTSTART']
            utend = hdu[0].header['UTSTOP']
            exptime = float(hdu[0].header['EXPTIME'])

            # Corrected DATE-OBS and UTSTART is the existing values plus the offset
            dateobs1 = Time(dateobs, format='isot', scale='utc')
            dateobs1 += dt
            utstart1 = str(dateobs1.value).split('T')[1]

            # Corrected UTEND is the UTSTART + exposure time
            utend1 = dateobs1 + Time(exptime, format='sec')
            utend1 = str(utend1.value).split('T')[1]

            # Corrected MJD-OBS is the corrected DATE-OBS converted to MJD
            mjdobs1 = dateobs1.jd - 2400000.5

            # Update the header:
            hdu2[0].header['DATE-OBS'] = str(dateobs1.value)
            hdu2[0].header['UTSTART'] = utstart1
            hdu2[0].header['UTSTOP'] = utend1
            hdu2[0].header['MJD-OBS'] = str(mjdobs1)

            # Save the corrected image file
            hdu2.writeto(image.replace('.fits', '_c.fits'), overwrite=True)
            print('Corrected timestatmp for image ' + image)


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Path to input image data directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    run(args)