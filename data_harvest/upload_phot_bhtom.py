from os import path
import argparse
import requests
import log_utils
import config_utils
from astropy.io import fits
import glob

DRY_RUN = 'False'

def upload_photometry(args):
    """
    Function to upload the photometry files for a set of images reduced with CCDphot

    This code expects to be pointed at a reduction directory containing a set of
    FITS images, for which there are corresponding output photometry ASCII files.
    While it is the ASCII photometry that will be uploaded, information gathered
    from the FITS image headers is also required.
    """

    config = config_utils.get_config(args.config_file)

    log = log_utils.start_day_log(config, 'bhtom_upload')

    # Sanity check the given reduction directory for FITS images and photometric output
    data_set = review_data(args, log)

    # Looping over all FITS images with corresponding photometry, gather
    # the required information from the image header, then upload the ASCII data table
    for image, phot_file in data_set.items():
        hdr = fits.getheader(path.join(args.red_dir, 'data', image))
        upload_phot_file(args, config, hdr, phot_file, log)

    log_utils.close_log(log)

def user_lut(config, params, log):
    """
    Function to fetch the BHTOM credentials of the user who requested the observation,
    based on the LCO user ID listed in the image header, and the configured credential
    information.

    Parameters:
        params dict     Header dictionary from an LCO FITS image.  Requires parameters
                        SITEID, TELID, INSTRUME and ORIGIN

    Returns:
        user_token  str BHTOM credentials for the user who took the data
    """

    if params['USERID'] in config['bhtom_users'].keys():
        log.info('Identified user')
        return config['bhtom_users'][params['USERID']]
    else:
        raise IOError('Data obtained by user unknown to BHTOM; no credentials to upload')

def oname_lut(params, log):
    """
    Function to look-up the ONAME (observatory name) for a facility as defined by BHTOM2,
    given parameters in the image header.

    Parameters:
        params dict     Header dictionary from an LCO FITS image.  Requires parameters
                        SITEID, TELID, INSTRUME and ORIGIN

    Returns:
        oname   str     BHTOM's identifier for the instrument that took the image
    """

    sites = {
        'tfn': 'Teide',
        'lsc': 'CTIO',
        'ogg': 'HO',
        'elp': 'MCD',
        'cpt': 'SAAO',
        'coj': 'SS',
    }
    site = sites[params['SITEID']]

    if '1m0a' in params['TELID']:
        aperture_class = '1m'
    elif '2m0a' in params['TELID']:
        aperture_class = '2m',
    elif '0m4a' in params['TELID']:
        aperture_class = '40cm'

    if 'fa' in params['INSTRUME']:
        instrument = '4K'
    elif 'sb' in params['INSTRUME']:
        instrument = 'SBIG6303'
    elif 'sq' in params['INSTRUME']:
        instrument = 'QHY600'

    oname = params['ORIGIN'] + '-' + site + '-' + aperture_class + '_' + instrument

    log.info('Identified image from ' + oname)

    return oname

def parse_filter_id(lco_filter):
    """
    Function to interpret the LCO filter names to BHTOM's syntax.
    This can be 'GaiaSP/any' or filter specific, e.g. 'GaiaSP/i'
    """

    bh_filter = 'GaiaSP/any'
    if lco_filter == 'ip':
        bh_filter = 'GaiaSP/i'
    elif lco_filter == 'rp':
        bh_filter = 'GaiaSP/i'
    elif lco_filter == 'gp':
        bh_filter = 'GaiaSP/g'

    return bh_filter

def review_data(args, log):
    """
    Function to review the FITS images and output photometry tables available for a
    single reduction directory.
    """

    # Create lists of the names of FITS images and photometry files available
    image_list = [path.basename(x) for x in glob.glob(path.join(args.red_dir, 'data', '*.fits'))]
    phot_list = [path.basename(x) for x in glob.glob(path.join(args.red_dir, 'phot', '*.pphot'))]

    log.info('Found ' + str(len(image_list)) + ' images in ' + args.red_dir)
    log.info('Found ' + str(len(phot_list)) + ' photometry data files in ' + args.red_dir)

    # Create a dictionary pairing the images with their corresponding
    # photometry tables if available
    data_set = {}
    for image in image_list:
        phot_file = (path.basename(image)).replace('.fits', '.pphot')
        if phot_file in phot_list:
            data_set[image] = phot_file

    log.info('Paired ' + str(len(data_set)) + ' phot files with images for upload')

    return data_set

def upload_phot_file(args, config, hdr, phot_file, log):
    """
    Function to upload the photometry from a single image, reduced by CCDPhot
    """

    oname = oname_lut(hdr, log)
    user_token, observer_name = user_lut(config, hdr, log)
    bh_filter = parse_filter_id(hdr['FILTER'])

    ur = {
        'target': hdr['OBJECT'],
        'filter': bh_filter,
        'data_product_type': 'photometry',
        'dry_run': DRY_RUN,
        'observatory': oname,
        'mjd': hdr['MJD-OBS'],
        'observer': observer_name,
    }

    response = requests.post(
        url='https://uploadsvc2.astrolabs.pl/upload/',
        headers={
            'Authorization': "Token " + str(user_token)
        },
        data=ur,
        files={'files': phot_file}
    )

    if log != None:
        log.info('Uploaded photometry file ' + phot_file + ' to BHTOM')
        log.info('with response: ' + repr(response.text))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help='Path to the configuration file')
    parser.add_argument('red_dir', help='Path to reduction directory')
    args = parser.parse_args()

    upload_photometry(args)
