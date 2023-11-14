from os import path
import json
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
import requests
from datetime import datetime
import numpy as np

def build_event_catalog(args):

    # Load the main event catalog, the table of Spitzer events and the list
    # of Gaia alerts
    events = load_object_list(args.main_file, ['ra', 'dec', 'baseline_mag'])
    print('Loaded '+str(len(events))+' events from the main list')
    (targets, target_index) = make_coordinate_index(events)

    spitzer_events = load_object_list(args.spitzer_file, ['ra', 'dec', 'rome_field'])
    gaia_alerts = load_object_list(args.gaia_file,
                ['ra', 'dec', 'Gaia_alert_class', 'Gaia_alert_comment',
                'ATel', 'rome_field'])

    # Combined the catalogs
    events = mark_spitzer_events(events, spitzer_events)
    events = mark_gaia_alerts(events, gaia_alerts, targets, target_index)

    # Retrieve known event parameters from OGLE, MOA
    events = fetch_ogle_event_parameters(events, targets, target_index)
    events = fetch_moa_event_parameters(events, targets, target_index)
    events = fetch_kmtnet_event_parameters(events, targets, target_index)

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

def mark_gaia_alerts(main_catalog, gaia_catalog, targets, target_index):
    """Function to crossmatch the main catalog against the catalog of alerts.
    """

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
            target_data['Gaia_alert_ra'] = data['ra']
            target_data['Gaia_alert_dec'] = data['dec']
            target_data['Gaia_alert_class'] = data['Gaia_alert_class']
            target_data['Gaia_alert_comment'] = data['Gaia_alert_comment']
            target_data['ATel'] = data['ATel']
            main_catalog[target_name] = target_data

        # If the Gaia classification is ULENS, and there is no known event
        # at those coordinates, then add the event to the main catalog
        elif 'ULENS' in data['Gaia_alert_class'] and d2d[0] > tol:
            target_data = {'ra': data['ra'],
                           'dec': data['dec'],
                           'baseline_mag': None,
                           'Gaia_alert_ID': gaia_name,
                           'Gaia_alert_ra': data['ra'],
                           'Gaia_alert_dec': data['dec'],
                           'Gaia_alert_class': data['Gaia_alert_class'],
                           'Gaia_alert_comment': data['Gaia_alert_comment'],
                           'ATel': data['ATel']}
            main_catalog[gaia_name] = target_data

    # Review the main catalog.  For those events where no Gaia alert was
    # issued, fill in blanks for the Gaia entries to ensure consistency of
    # formatting
    for target_name, data in main_catalog.items():
        if 'Gaia_alert_ID' not in data.keys():
            data['Gaia_alert_ID'] = None
            data['Gaia_alert_ra'] = None
            data['Gaia_alert_dec'] = None
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
    """The event object list was compiled by Yiannis Tsapras, using software that loaded
    the list of OGLE events first, followed by MOA and then KMTNet.  If an event was
    identified in more than one survey, it is only entered into the list under the
    naming convention of the first survey in the sequence, so duplicated detections
    are not flagged.

    Tables of cross-identified events are available from the surveys:
    For MOA https://www.massey.ac.nz/~iabond/moa/alertXXXX/fetchtxt.php?path=moa/alertXXXX/moa2ogle.txt
    For KMTNet a webpage scrapper code was used:
    https://kmtnet.kasi.re.kr/ulens/event/2017/
    """
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

def fetch_ogle_event_parameters(catalog, targets, target_index):

    BROKER_URL = 'https://www.astrouw.edu.pl/ogle/ogle4/ews'
    tol = 2.0 / 3600.0 * u.deg

    # Fetch the online OGLE event catalogs for the years of the ROME/REA survey
    survey_years = [2017, 2018, 2019]
    ogle_events = {}
    for year in survey_years:
        par_file_url = path.join(BROKER_URL, str(year), 'lenses.par')
        response = requests.request('GET', par_file_url)
        if response.status_code == 200:
            for line in response.iter_lines():
                line = str(line)
                if 'StarNo' not in line and len(line) > 5:  # Skip the file header
                    entries = line.split()
                    name = 'OGLE-' + entries[0].replace("b'", "")
                    ra = entries[3]
                    dec = entries[4]
                    to = float(entries[5])
                    te = float(entries[7])
                    uo = float(entries[8])
                    s = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
                    ogle_events[name] = {'RA': ra, 'Dec': dec, 't0': to, 'tE': te, 'u0': uo}

                    if name in catalog.keys():
                        target_data = catalog[name]
                        target_data['t0'] = ogle_events[name]['t0']
                        target_data['tE'] = ogle_events[name]['tE']
                        target_data['u0'] = ogle_events[name]['u0']
                        catalog[name] = target_data

                    else:
                        (idx, d2d, d3d) = s.match_to_catalog_sky(targets)
                        if d2d[0] <= tol:
                            target_name = target_index[int(idx)]
                            target_data = catalog[target_name]
                            target_data['OGLE_alert_ID'] = name
                            target_data['OGLE_alert_ra'] = ogle_events[name]['RA']
                            target_data['OGLE_alert_dec'] = ogle_events[name]['Dec']
                            target_data['OGLE_t0'] = ogle_events[name]['t0']
                            target_data['OGLE_tE'] = ogle_events[name]['tE']
                            target_data['OGLE_u0'] = ogle_events[name]['u0']
                            catalog[target_name] = target_data

    return catalog

def fetch_moa_event_parameters(catalog, targets, target_index):

    BROKER_URL = 'https://www.massey.ac.nz/~iabond/moa/'
    tol = 2.0 / 3600.0 * u.deg

    # Fetch the online MOA event catalogs for the years of the ROME/REA survey
    survey_years = [2017, 2018, 2019]
    moa_events = {}
    for year in survey_years:
        par_file_url = path.join(BROKER_URL, 'alert'+str(year), 'alert.php')
        response = requests.request('GET', par_file_url)
        if response.status_code == 200:
            for line in response.iter_lines():
                line = str(line)
                if str(year) in line:
                    entries = line.replace('<',' ').replace('>',' ').replace('\n','').split()
                    if len(entries) == 38:
                        name = 'MOA-'+entries[5]
                        (date,dayfrac) = entries[18].split('.')
                        to = Time(datetime.strptime(date, "%Y-%b-%d")).jd + float(dayfrac)
                        te = float(entries[22])
                        amax = float(entries[26])
                        s = SkyCoord(entries[10], entries[14], frame='icrs', unit=(u.hourangle, u.deg))
                        moa_events[name] = {'t0': to, 'tE': te, 'amax': amax, 'ra': s.ra.deg, 'dec': s.dec.deg}

                        if name in catalog.keys():
                            target_data = catalog[name]
                            target_data['ra'] = moa_events[name]['ra']
                            target_data['dec'] = moa_events[name]['dec']
                            target_data['t0'] = moa_events[name]['t0']
                            target_data['tE'] = moa_events[name]['tE']
                            target_data['amax'] = moa_events[name]['amax']
                            catalog[name] = target_data

                        else:
                            (idx, d2d, d3d) = s.match_to_catalog_sky(targets)
                            if d2d[0] <= tol:
                                target_name = target_index[int(idx)]
                                target_data = catalog[target_name]
                                target_data['MOA_alert_ID'] = name
                                target_data['MOA_ra'] = moa_events[name]['ra']
                                target_data['MOA_dec'] = moa_events[name]['dec']
                                target_data['MOA_t0'] = moa_events[name]['t0']
                                target_data['MOA_tE'] = moa_events[name]['tE']
                                target_data['MOA_amax'] = moa_events[name]['amax']
                                catalog[target_name] = target_data

    return catalog

def fetch_kmtnet_event_parameters(catalog, targets, target_index):

    BROKER_URL = 'https://kmtnet.kasi.re.kr/~ulens/event/'
    tol = 2.0 / 3600.0 * u.deg

    # Fetch the online KMTNet event catalogs for the years of the ROME/REA survey
    survey_years = [2017, 2018, 2019]
    kmtnet_events = {}
    for year in survey_years:
        par_file_url = path.join(BROKER_URL, str(year), 'listpage.dat')
        response = requests.request('GET', par_file_url)
        if response.status_code == 200:
            for line in response.iter_lines():
                line = str(line)
                if str(year) in line:
                    entries = line.replace("b'",'').replace('\n','').split()
                    name = entries[0]
                    try:
                        if year == 2017:
                            s = SkyCoord(entries[3], entries[4], frame='icrs', unit=(u.hourangle, u.deg))
                            to = float(entries[5]) + 2450000.0
                            te = float(entries[6])
                            uo = float(entries[7])
                            kmtnet_events[name] = {'t0': to, 'tE': te, 'u0': uo, 'ra': s.ra.deg, 'dec': s.dec.deg}
                        else:
                            s = SkyCoord(entries[4], entries[5], frame='icrs', unit=(u.hourangle, u.deg))
                            to = float(entries[6]) + 2450000.0
                            te = float(entries[7])
                            uo = float(entries[8])
                            kmtnet_events[name] = {'t0': to, 'tE': te, 'u0': uo, 'ra': s.ra.deg, 'dec': s.dec.deg}
                    except ValueError:
                        kmtnet_events[name] = {'t0': np.nan, 'tE': np.nan, 'u0': np.nan, 'ra': 'None', 'dec': 'None'}

                    if name in catalog.keys():
                        target_data = catalog[name]
                        target_data['t0'] = kmtnet_events[name]['t0']
                        target_data['tE'] = kmtnet_events[name]['tE']
                        target_data['u0'] = kmtnet_events[name]['u0']
                        target_data['ra'] = kmtnet_events[name]['ra']
                        target_data['dec'] = kmtnet_events[name]['dec']
                        catalog[name] = target_data

                    else:
                        (idx, d2d, d3d) = s.match_to_catalog_sky(targets)
                        if d2d[0] <= tol:
                            target_name = target_index[int(idx)]
                            target_data = catalog[target_name]
                            target_data['KMTNet_alert_ID'] = name
                            target_data['KMTNet_ra'] = kmtnet_events[name]['ra']
                            target_data['KMTNet_dec'] = kmtnet_events[name]['dec']
                            target_data['KMTNet_t0'] = kmtnet_events[name]['t0']
                            target_data['KMTNet_tE'] = kmtnet_events[name]['tE']
                            target_data['KMTNet_u0'] = kmtnet_events[name]['u0']
                            catalog[target_name] = target_data

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
