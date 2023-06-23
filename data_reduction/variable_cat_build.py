from os import path
import json
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
import event_cat_build

def build_variable_catalog(args):

    # Load the variable catalogs from the different surveys.  The OGLE catalog
    # is treated as the main catalog overall
    variables = event_cat_build.load_object_list(args.ogle_var_file)
    vvv_variables = event_cat_build.load_object_list(args.vvv_var_file)
    gaia_alerts = event_cat_build.load_object_list(args.gaia_file)

    # Combine the information from the VVV catalog and Gaia alerts
    variables = add_vvv_variables(variables, vvv_variables)
    variables = event_cat_build.mark_gaia_alerts(variables, gaia_alerts, main_entry_length=5)

    # Output the revised main_catalog
    event_cat_build.output_catalog(args.output_file, variables)

def add_vvv_variables(main_catalog, vvv_catalog, main_entry_length=5):
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
        v = SkyCoord(vvv_data[0], vvv_data[1], frame='icrs', unit=(u.deg, u.deg))
        (idx, d2d, d3d) = v.match_to_catalog_sky(targets)
        if d2d[0] <= tol:
            target_name = target_index[int(idx)]
            target_data = main_catalog[target_name]
            target_data.append(vvv_name)
            target_data.append(vvv_data[2])
            target_data.append(vvv_data[3])
            target_data.append(vvv_data[4])
            main_catalog[target_name] = target_data

        else:
            target_data = [vvv_data[0], vvv_data[1]]
            add_length = main_entry_length - 2
            target_data += [None]*add_length
            target_data += [vvv_name, vvv_data[2], vvv_data[3], vvv_data[4]]
            main_catalog[vvv_name] = target_data

    # Review the main catalog and back-fill the entries for any unmatched
    # objects, to ensure consistent formatting
    for target_name, data in main_catalog.items():
        if len(data) == 5:
            data.append('no_vvv_entry')
            data.append('no_vvv_class')
            data.append('no_vvv_period')
            data.append('no_gaiaEDR3')
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
