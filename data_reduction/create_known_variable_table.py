from os import path
import json
import argparse
import numpy as np

def create_variable_lookup_table(args):
    """
    Function to convert the JSON file of known variables and events from each field into
    a look-up table indexed by ROME field ID
    """

    # Load the tables of known events and variables for this field
    event_catalog = load_json_target_catalog(args.event_catalog_file)
    variable_catalog = load_json_target_catalog(args.variable_catalog_file)

    # Reviewing all entries of known events and variables, create a look-up table indexed by
    # ROME field ID.  Combine both events and variables together.
    LUT = {}
    LUT = extract_rome_star_matches(LUT, event_catalog, catalog_type='event')
    LUT = extract_rome_star_matches(LUT, variable_catalog, catalog_type='variable')

    # Output LUT
    json_data = json.dumps(LUT, indent=4)
    with open(args.lookup_file, 'w') as write_file:
        write_file.write(json_data)
        write_file.close()

def extract_rome_star_matches(LUT, catalog, catalog_type='event'):
    """
    Function to take a JSON catalog of events or variables matched against the ROME catalog.
    If a known variable has a match in the ROME catalog, its entry will include a 'rome_stars'
    dictionary in the catalog.

    The match radius typically used for cross-matching when the input catalogs were generated was 2arcsec.
    It is still possible that the target has multiple matches, in which case the closest match is marked
    as the variable.

    The catalog_type parameter indicates whether the input catalog flags events or variables,
    so that the corresponding identifiers can be stored with the correct keyword.
    """

    def extract_keyword_entry(key, entry):
        value = None
        if key in entry.keys():
            if 'null' not in str(entry[key]).lower() and 'false' not in str(entry[key]).lower():
                value = entry[key]
        return value

    for name, entry in catalog.items():
        if 'rome_stars' in entry.keys():
            separations = [x['separation_deg'] for x in entry['rome_stars']]
            if len(separations) > 0:
                jdx = np.argsort(separations)[0]
                # Check to see if we already have an entry for this object;
                # this can happen if multiple surveys discover the same event for instance
                if entry['rome_stars'][jdx]['field_id'] in LUT.keys():
                    lut_entry = LUT[entry['rome_stars'][jdx]['field_id']]
                else:
                    lut_entry = {
                        'ogle_event_id': None,
                        'ogle_variable_id': None,
                        'moa_event_id': None,
                        'kmtnet_event_id': None,
                        'spitzer_event': None,
                        'vvv_variable_id': None
                    }

                # Collate data for events, potentially detected by multiple surveys
                if catalog_type == 'event':
                    if 'OGLE' in name:
                        lut_entry['ogle_event_id'] = name
                        lut_entry['kmtnet_event_id'] = extract_keyword_entry('KMTNet_alert_ID', entry['target_data'])

                    elif 'KMT' in name:
                        lut_entry['kmtnet_event_id'] = name
                    lut_entry['spitzer_event'] = extract_keyword_entry('spitzer_target', entry['target_data'])
                    lut_entry['moa_event_id'] = extract_keyword_entry('MOA_alert_ID', entry['target_data'])

                # Handle different formating of variable star catalog
                else:
                    if 'OGLE' in name:
                        lut_entry['ogle_variable_id'] = name
                        lut_entry['vvv_variable_id'] = extract_keyword_entry('VVV_name', entry['target_data'])
                    else:
                        lut_entry['vvv_variable_id'] = name

                LUT[entry['rome_stars'][jdx]['field_id']] = lut_entry

    return LUT

def load_json_target_catalog(catalog_path):
    """Function to load a catalog of selected objects in JSON format"""

    with open(catalog_path, "r") as read_file:
        data = json.load(read_file)

    return data

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('event_catalog_file', help='Path to event catalog JSON file')
    parser.add_argument('variable_catalog_file', help='Path to variable catalog JSON file')
    parser.add_argument('lookup_file', help='Path to output look-up table file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    create_variable_lookup_table(args)