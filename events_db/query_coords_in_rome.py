from sys import argv
from os import path
import requests
import json
from astropy.coordinates import SkyCoord
from astropy import units as u

def query_coords_in_rome(ra_str, dec_str, debug=False):
    """Function to query ROME's public API to verify whether a specific
    set of coordinates lie within the survey footprint.
    
    Accepts sexigesimal coordinate strings
    """
    
    c = SkyCoord(ra_str+' '+dec_str, frame='icrs', unit=(u.hourangle, u.deg))
    
    headers = {}
    
    if debug:
        url = 'http://127.0.0.1:8000/db/'
        
    else:
        
        url = 'https://robonet.lco.global/db/'
    
    url = path.join(url,'query_event_in_survey',
                    'ra='+str(c.ra.value)+'&dec='+str(c.dec.value))

    print(url)
    
    response = requests.get(url, headers=headers, timeout=20)
    
    lines = response.text.split('\n')
    for l in lines:
        if len(l) > 0:
            result = l
    
    print(result)
    
if __name__ == '__main__':
    
    if len(argv) >= 3:
        
        ra_str = argv[1]
        dec_str = argv[2]
        
    else:
        
        ra_str = raw_input('Please enter the RA in sexigesimal format: ')
        dec_str = raw_input('Please enter the Dec in sexigesimal format: ')
        
    if len(argv) == 4:
        
        debug = True
    
    else:
    
        debug = False
        
    query_coords_in_rome(ra_str, dec_str, debug)
    