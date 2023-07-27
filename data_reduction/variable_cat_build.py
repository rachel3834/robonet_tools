from os import path
import json
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
import event_cat_build

def build_variable_catalog(args):

    # Load the variable catalogs from the different surveys.  The OGLE catalog
    # is treated as the main catalog overall
    variables = event_cat_build.load_object_list(args.ogle_var_file,
                ['ra','dec','ogle_class','ogle_subclass','rome_field'])
    vvv_variables = event_cat_build.load_object_list(args.vvv_var_file,
                ['ra', 'dec', 'VVV_class', 'VVV_period', 'Gaia_EDR3_ID', 'rome_field'])
    gaia_alerts = event_cat_build.load_object_list(args.gaia_file,
                ['ra', 'dec', 'Gaia_alert_class', 'Gaia_alert_comment',
                'ATel', 'rome_field'])

    # Combine the information from the VVV catalog and Gaia alerts
    variables = add_vvv_variables(variables, vvv_variables)
    variables = event_cat_build.mark_gaia_alerts(variables, gaia_alerts)

    # Output the revised main_catalog
    event_cat_build.output_catalog(args.output_file, variables)

def add_vvv_variables(main_catalog, vvv_catalog):
    """Add variables to the main catalog from the VVV list.  If the same star
    is already present in the main catalog, then the data from VVV is combined
    with the existing entry.  Otherwise, a new entry is made.
    The format of the VVV catalog is:
    [radeg, decdeg, classification, period, gaiaEDR3, romefield]
    """

    # Build a coordinate index of the main catalog:
    (targets, target_index) = event_cat_build.make_coordinate_index(main_catalog)

    # Crossmatch all of the VVV variables against this target index,
    # combining or adding entries as appropriate
    tol = 2.0/3600.0 * u.deg
    for vvv_name, vvv_data in vvv_catalog.items():
        v = SkyCoord(vvv_data['ra'], vvv_data['dec'], frame='icrs', unit=(u.deg, u.deg))
        (idx, d2d, d3d) = v.match_to_catalog_sky(targets)
        if d2d[0] <= tol:
            target_name = target_index[int(idx)]
            target_data = main_catalog[target_name]
            target_data['VVV_name'] = vvv_name
            target_data['VVV_ra'] = vvv_data['ra']
            target_data['VVV_dec'] = vvv_data['dec']
            target_data['VVV_class'] = vvv_data['VVV_class']
            target_data['VVV_period'] = vvv_data['VVV_period']
            target_data['Gaia_EDR3_ID'] = vvv_data['Gaia_EDR3_ID']
            main_catalog[target_name] = target_data

        else:
            target_data = {'ra': vvv_data['ra'],
                           'dec': vvv_data['dec'],
                           'ogle_class': None,
                           'ogle_subclass': None,
                           'VVV_name': vvv_name,
                           'VVV_ra': vvv_data['ra'],
                           'VVV_dec': vvv_data['dec'],
                           'VVV_class': vvv_data['VVV_class'],
                           'VVV_period': vvv_data['VVV_period'],
                           'Gaia_EDR3_ID': vvv_data['Gaia_EDR3_ID']}
            main_catalog[vvv_name] = target_data

    # Review the main catalog and back-fill the entries for any unmatched
    # objects, to ensure consistent formatting
    for target_name, data in main_catalog.items():
        if len(data) == 5:
            data['VVV_name'] = None
            data['VVV_class'] = None
            data['VVV_period'] = None
            data['Gaia_EDR3_ID'] = None
            main_catalog[target_name] = data

    return main_catalog

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ogle_var_file", help="Path to OGLE variable file in JSON format", type=str)
    parser.add_argument("vvv_var_file", help="Path to VVV variables file in JSON format", type=str)
    parser.add_argument("gaia_file", help="Path to the table of Gaia events", type=str)
    parser.add_argument("output_file", help="Path to output file of matching objects, in JSON", type=str)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_variable_catalog(args)
