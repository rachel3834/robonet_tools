from os import path
import json
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u

def build_event_catalog(args):

    # Load the main event catalog, the table of Spitzer events and the list
    # of Gaia alerts
    events = load_object_list(args.main_file, ['ra', 'dec', 'baseline_mag'])
    spitzer_events = load_object_list(args.spitzer_file, ['ra', 'dec', 'rome_field'])
    gaia_alerts = load_object_list(args.gaia_file,
                ['ra', 'dec', 'Gaia_alert_class', 'Gaia_alert_comment',
                'ATel', 'rome_field'])

    # Combined the catalogs
    events = mark_spitzer_events(events, spitzer_events)
    events = mark_gaia_alerts(events, gaia_alerts)

    # Output revised catalog
    output_catalog(args.output_file, events)

def output_catalog(output_file, main_catalog):

    json_data = json.dumps(main_catalog, indent=4)
    with open(output_file, 'w') as write_file:
        write_file.write(json_data)
        write_file.close()

def make_coordinate_index(catalog):
    """Function to produce a coordinate-based index of an input catalog.
    This expects an input catalog in dictionary format, where:
    {target1: [ra, dec, data*]}
    """

    # Make a coordinate-based index of all variables in the main catalog
    ra = []
    dec = []
    target_index = {}
    i = -1
    for target, data in catalog.items():
        ra.append(data['ra'])
        dec.append(data['dec'])
        i += 1
        target_index[i] = target
    targets = SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg))

    return targets, target_index

def mark_gaia_alerts(main_catalog, gaia_catalog, main_entry_length=3):
    """Function to crossmatch the main catalog against the catalog of alerts.
    Parameter main_entry_length indicates the expected default number of items
    in the list of data for each main catalog entry, prior to the addition of
    other data.
    """

    # Build an index of the coordinates of all objects in tha main catalog.
    # Keep track of the index numbers of the targets in the list for future
    # reference
    (targets, target_index) = make_coordinate_index(main_catalog)

    # Search for each Gaia event in the main target catalog and append the
    # Gaia alert data to the corresponding entry in the main catalog, if there
    # is one
    tol = 2.0/3600.0 * u.deg
    for gaia_name, data in gaia_catalog.items():
        g = SkyCoord(data['ra'], data['dec'], frame='icrs', unit=(u.deg, u.deg))
        (idx, d2d, d3d) = g.match_to_catalog_sky(targets)
        if d2d[0] <= tol:
            target_name = target_index[int(idx)]
            target_data = main_catalog[target_name]
            target_data['Gaia_alert_ID'] = gaia_name
            target_data['Gaia_alert_class'] = data['Gaia_alert_class']
            target_data['Gaia_alert_comment'] = data['Gaia_alert_comment']
            target_data['ATel'] = data['ATel']
            main_catalog[target_name] = target_data

        # If the Gaia classification is ULENS, and there is no known event
        # at those coordinates, then add the event to the main catalog
        else:
            target_data = {'ra': data[0],
                           'dec': data[1],
                           'baseline_mag': None,
                           'Gaia_alert_ID': gaia_name,
                           'Gaia_alert_class': data['Gaia_alert_class'],
                           'Gaia_alert_comment': data['Gaia_alert_comment'],
                           'ATel': data['ATel']}
            main_catalog[gaia_name] = target_data

    # Review the main catalog.  For those event where no Gaia alert was
    # issued, fill in blanks for the Gaia entries to ensure consistency of
    # formatting
    for target_name, data in main_catalog.items():
        if 'Gaia_alert_ID' not in data.keys():
            data['Gaia_alert_ID'] = None
            data['Gaia_alert_class'] = None
            data['Gaia_alert_comment'] = None
            data['ATel'] = None
            main_catalog[target_name] = data

    return main_catalog

def mark_spitzer_events(main_catalog, add_catalog):
    """Function to crossmatch two catalogs on event name and combine the
    results"""

    combined_catalog = {}
    for target_name, data in main_catalog.items():
        if target_name in add_catalog.keys():
            add_data = add_catalog[target_name]
            data['spitzer_target'] = True
        else:
            data['spitzer_target'] = False
        combined_catalog[target_name] = data

    return combined_catalog

def load_object_list(input_file, key_list):

    if not path.isfile(input_file):
        raise IOError('Cannot find input list of target coordinates at '
                        + input_file)

    with open(input_file, "r") as read_file:
        data = json.load(read_file)

    # Reformat for readability
    catalog = {}
    for name, entry in data.items():
        target_data = {}
        for i,key in enumerate(key_list):
            target_data[key] = entry[i]
        catalog[name] = target_data

    return catalog

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("main_file", help="Path to the main JSON file of events", type=str)
    parser.add_argument("spitzer_file", help="Path to the table of Spitzer events", type=str)
    parser.add_argument("gaia_file", help="Path to the table of Gaia events", type=str)
    parser.add_argument("output_file", help="Path to output file of matching objects, in JSON", type=str)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_event_catalog(args)
