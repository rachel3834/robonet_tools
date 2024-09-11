from astropy.io import votable
from os import path
import argparse

def parse_gaia_votable(args):
    """
    Function to parse a VO table file downloaded from
    https://gaia.ari.uni-heidelberg.de/cone.html
    as a cone search selecting ALL COLUMNS
    """

    if not path.isfile(args.file_path):
        raise IOError('Cannot find the input VO Table file ' + args.file_path)

    results = votable.parse_single_table(args.file_path)

    table = results.to_table()
    print('Got ' + str(len(table)) + ' results, storing first entry')

    param_map = field_mapping()

    results = {}
    for field, key in param_map.items():
        results[key] = table[field].value[0]
    print(results)

def field_mapping():
    """
    Function to provide the mapping of the VO Table fields to MOP parameter names
    """

    mapping = {
        'source_id': 'gaia_source_id',
        'phot_g_mean_mag': 'gmag',
        'phot_rp_mean_mag': 'rpmag',
        'phot_bp_mean_mag': 'bpmag',
        'bp_rp': 'bprp',
        'ebpminrp_gspphot':'reddening_bprp',
        'azero_gspphot': 'extinction_g',
        'distance_gspphot': 'distance',
        'teff_gspphot': 'teff',
        'logg_gspphot': 'logg',
        'ruwe': 'ruwe'
    }

    return mapping
def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('file_path', help='Path to input VO Table file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    parse_gaia_votable(args)