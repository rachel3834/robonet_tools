import numpy as np
import json
from astropy import units
from astropy.coordinates import Angle, SkyCoord

# Function that prints out the event details (name, rome-field, ra, dec, baseline) for any event above a given baseline mag threshold in every ROME field.
def get_events(events_file='romerea_events.json', 
               field_events_file='events_in_each_field.json', 
               thresh=17.5):
    
    # Read the files
    with open(events_file, 'r') as f:
        romerea_events = json.load(f)
    
    with open(field_events_file, 'r') as f:
        events_in_each_field = json.load(f)
    
    # Loop over all 20 ROME/REA fields
    for field in np.arange(20)+1:
        # Check each event 
        for event in events_in_each_field[str(field)]:
            yr = event.split('-')[1]
            ev_nr = event.split('-')[-1]
            coords_ra = Angle(romerea_events[event][0]/15.0, units.hourangle)
            coords_dec = Angle(romerea_events[event][1], units.deg)
            baseline = romerea_events[event][-1]
            params = {}
            params['name'] = event
            params['rome-field'] = field
            params['ra'] = coords_ra.to_string(unit=units.hourangle, sep=':')
            params['dec'] = coords_dec.to_string(unit=units.deg, sep=':')
            params['baseline'] = baseline
            format_url = 'http://ogle.astrouw.edu.pl/ogle4/ews/YYYY/blg-NNNN.html'.replace('YYYY',yr)
            format_url = format_url.replace('NNNN', ev_nr)
            params['ogle_url'] = format_url
            # Threshold check
            if baseline < thresh:
                print (params)

# Run it
get_events(thresh=17.5)
