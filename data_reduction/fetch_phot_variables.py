from os import path
import json
from pyDANDIA import crossmatch, search_crossmatch_bulk
import argparse

def fetch_phot_for_variables(args):

    # Load catalog of variable stars
    catalog = search_crossmatch_bulk.load_object_list(args.variable_cat)

    # Load the crossmatch table
    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.xmatch_file, log=None)

    # For each variable, using the closest match in the catalog to identify the most likely star
    # from ROME/REA, extract the reference magnitude in each band
    for name, var_data in catalog.items():

        # Find the closest candidate ROME star
        min_sep = 1e6
        j = -1
        for candidate in var_data['rome_stars']:
            if candidate['separation_deg'] < min_sep:
                min_sep = candidate['separation_deg']
                j = candidate['field_id']

        # Extract the reference image magnitudes if the candidate was measured in the
        # primary reference for the field.
        if j > 0:
            var_data['cal_g_mag_lsc_doma'] = xmatch.stars['cal_g_mag_lsc_doma'][j]
            var_data['cal_r_mag_lsc_doma'] = xmatch.stars['cal_r_mag_lsc_doma'][j]
            var_data['cal_i_mag_lsc_doma'] = xmatch.stars['cal_i_mag_lsc_doma'][j]
            catalog[name] = var_data
            print(name, var_data)

    # Output the amended catalog:
    search_crossmatch_bulk.output_search_results(args.output_file, catalog)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('variable_cat', help="Path to catalog of variable stars from a field")
    parser.add_argument('xmatch_file', help="Path to the crossmatch table for the field")
    parser.add_argument('output_file', help="Path to output file")
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    fetch_phot_for_variables(args)