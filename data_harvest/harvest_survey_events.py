# -*- coding: utf-8 -*-
"""
A set of routines to extract the coordinates of known survey events
and compare them with a specified footprint 

@author: ytsapras
"""

import numpy as np
from astropy import units
from astropy.coordinates import Angle, SkyCoord
from urllib import request, parse
from sys import exit
from bs4 import BeautifulSoup
from tqdm import tqdm
import json

#############################################################
##################### User definitions ######################
# A dictionary with the survey names and websites
survey_dict = {'OGLE':'http://ogle.astrouw.edu.pl/ogle4/ews/',
               'MOA':'https://www.massey.ac.nz/~iabond/moa/',
               'KMTNet':'http://kmtnet.kasi.re.kr/ulens/event/'}

# A list of the years to look for events
years = [2017, 2018, 2019]
#############################################################
#############################################################

#############################################################
def construct_url_list(survey_dict, years):
    '''
    Contruct the appropriate url strings for each year.
    Returns a list of urls.
    '''
    urls = []
    for year in years:
        for survey in survey_dict.keys():
            if survey == 'OGLE':
                urlStr = survey_dict[survey]+str(year)+'/ews.html'
            elif survey == 'MOA':
                urlStr = survey_dict[survey]+'alert'+str(year)+'/alert.php'
            elif survey == 'KMTNet':
                urlStr = survey_dict[survey]+str(year)+'/'
            else:
                print(('Undefined survey: %s' % survey))
                exit()
            
            urls.append(urlStr)
    
    return urls

#############################################################
def crawl_table(htmlsoup):
    ''' 
    Function to parse out the rows of an HTML table.
    '''
    return [ [ ''.join(col(text=True)) for col in row.findAll('td') ]\
            for row in htmlsoup.find('table').findAll('tr') ]

#############################################################
def check_exists(event_dict, new_radeg, new_decdeg):
    '''
    Cross-survey identification check.
    Check if an event at these coordinates already exists in the dictionary.
    If the coordinates are already in the dictionary, returns True
    '''
    new_coords = SkyCoord(str(new_radeg)+' '+str(new_decdeg), unit=(units.hourangle, units.deg),
                          frame='icrs')
    
    match_found = False
    
    for ra, dec, ibase in event_dict.values():
        known_coords = SkyCoord(str(ra)+' '+str(dec), unit=(units.hourangle, units.deg), 
                                frame='icrs')
        # Calculate separation in arcsec
        separation = new_coords.separation(known_coords).arcsec
        if separation < 2.5:
            match_found = True
            break
    return match_found
   
#############################################################
def romerea_check(radeg, decdeg):
    '''
    Check if the input coordinates are within the ROME/REA footprint.
    Returns field number (1 to 20) or -1 (if outside footprint). 
    '''
    lhalf = 0.220833333333  # 26.5/(120.)
    field, rate = -1, -1
    fields = [[267.835895375, -30.0608178195, 64.0],
              [269.636745458, -27.9782661111, 49.0],
              [268.000049542, -28.8195573333, 46.0],
              [268.180171708, -29.27851275, 58.0],
              [268.35435, -30.2578356389, 64.0],
              [268.356124833, -29.7729819283, 90.0],
              [268.529571333, -28.6937071111, 72.0],
              [268.709737083, -29.1867251944, 83.0],
              [268.881108542, -29.7704673333, 83.0],
              [269.048498333, -28.6440675, 75.0],
              [269.23883225, -29.2716684211, 70.0],
              [269.39478875, -30.0992361667, 42.0],
              [269.563719375, -28.4422328996, 49.0],
              [269.758843, -29.1796030365, 67.0],
              [269.78359875, -29.63940425, 61.0],
              [270.074981708, -28.5375585833, 61.0],
              [270.81, -28.0978333333, -99.0],
              [270.290886667, -27.9986032778, 52.0],
              [270.312763708, -29.0084241944, 48.0],
              [270.83674125, -28.8431573889, 49.0]]
    
    for idx in range(len(fields)):
        if radeg < fields[idx][0] + lhalf and\
           radeg > fields[idx][0] - lhalf and\
           decdeg < fields[idx][1] + lhalf and\
           decdeg > fields[idx][1] - lhalf:
               return idx+1#, fields[idx][2]
    
    return field#, rate

# Separation between two events (example)
#k1 = known_events['2018-BLG-1659']
#k2 = known_events['2018-BLG-1663']
#c1 = SkyCoord(k1[0],k1[1], frame='icrs', unit=(units.hourangle, units.deg))
#c2 = SkyCoord(k2[0],k2[1], frame='icrs', unit=(units.hourangle, units.deg))
#ang = c1.separation(c2).value # In decimal degrees

#############################################################
def collect_known_events(urls_list):
    ''' 
        Will go over the given survey websites (and years) and 
        collect all unique events. Returns a dictionary 
        of event names, sky coordinates and baseline magnitudes.
    '''
    known_events = {}
    # Loop over all given urls
    print("Harvesting known events ...")
    for urlStr in urls_list:
        # Open the url and extract the event information
        try:
            fileHandle = request.urlopen(urlStr)
            htmlFile = fileHandle.read()
            fileHandle.close()
        except IOError:
            print(('Cannot open URL %s for reading' % urlStr))
        
        soup = BeautifulSoup(htmlFile, features="lxml")
        table_contents = crawl_table(soup)
        
        # Populate the event dictionary with the events
        # OGLE or MOA or KMTNet specific
        if 'ogle.astrouw.edu.pl' in urlStr:
            for line in tqdm(table_contents[1:]):
                #print(("Adding OGLE-"+line[1].split("\n")[0]+" to dictionary."))
                # The event name is the 1st column in the line
                event = 'OGLE-'+line[1].split("\n")[0]
                # The RA and Dec are the 4th and 5th columns in the line
                ra = line[4].split("\n")[0]
                radeg = Angle(ra, unit=units.hourangle).value * 15.0
                dec = line[5].split("\n")[0]
                decdeg = Angle(dec, unit=units.deg).value
                # Ibase is the last column in the line
                ibase = float(line[-1].split("\n")[0])
                # If it is not an already known OGLE event, add it
                # Warning: check_exists is more accurate but very slow
                #if check_exists(known_events, radeg, decdeg) == False:
                known_events[event] = [radeg, decdeg, ibase]
        elif 'massey.ac.nz' in urlStr:
            crossmatch = ''
            # Get the OGLE/MOA crossmatch information
            try:
                matchurl = 'https://www.massey.ac.nz/~iabond/moa/alertXXXX/fetchtxt.php?path=moa/alertXXXX/moa2ogle.txt'.replace('alertXXXX',urlStr.split('/')[-2])
                crossmatch = BeautifulSoup(request.urlopen(matchurl).read(), features="lxml").prettify()
            except:
                print('Failed to find a moa2ogle.txt match file.')
            for line in tqdm(table_contents[1:]):
                # The event name is the 0th column in the line
                event = 'MOA-'+line[0].split("\n")[0]
                # The RA and Dec are the 1st and 2nd columns in the line
                ra = line[1].split("\n")[0]
                radeg = Angle(ra, unit=units.hourangle).value * 15.0
                dec = line[2].split("\n")[0]
                decdeg = Angle(dec, unit=units.deg).value
                # Ibase is the 6th column in the line
                ibase = float(line[6].split("\n")[0])
                # If it is not an already known OGLE event, add it
                # Warning: check_exists is more accurate but very slow
                #if check_exists(known_events, radeg, decdeg) == False:
                if event not in crossmatch:
                    #print(("Adding MOA-"+line[0].split("\n")[0]+" to dictionary."))
                    known_events[event] = [radeg, decdeg, ibase]
        elif '/kmtnet.kasi.re.kr/' in urlStr:
            for line in tqdm(table_contents[1:]):
                # The KMTNet HTML Table lines have special characters.
                # They also have different column formats for each year.
                # Remove special characters from line
                line = [i.strip() for i in line]
                # The event name is the 0th column in the line
                event = line[0].split("\n")[0]
                # The RA, Dec and ibase are the
                # 2017: 8th, 9th and 7th columns in the line
                # 2018: 5th, 6th and 11th columns in the line
                # 2019: 4th, 5th and 10th columns in the line
                if 'KMT-2017-' in event:
                    ra = line[4].split("\n")[0]
                    dec = line[5].split("\n")[0]
                    ibase = float(line[7].split("\n")[0])
                elif 'KMT-2018-' in event:
                    ra = line[5].split("\n")[0]
                    dec = line[6].split("\n")[0]
                    ibase = float(line[11].split("\n")[0])
                elif 'KMT-2019-' in event:
                    ra = line[4].split("\n")[0]
                    dec = line[5].split("\n")[0]
                    ibase = float(line[10].split("\n")[0])
                radeg = Angle(ra, unit=units.hourangle).value * 15.0
                decdeg = Angle(dec, unit=units.deg).value
                # The last KMTNet column in the line is always the 
                # crossmatch information
                xmatch = line[-1].split("\n")[0]
                # If it is not an already known OGLE or MOA event, add it
                # Warning: check_exists is more accurate but very slow
                #if check_exists(known_events, radeg, decdeg) == False:
                if xmatch == '':
                    #print(("Adding "+line[0].split("\n")[0]+" to dictionary."))
                    known_events[event] = [radeg, decdeg, ibase]
        else:
            print(('Undefined survey url string: %s' % urlStr))
            exit()
    
    return known_events

#############################################################    
def select_events_within_footprint(known_events_dictionary, romerea=1):
    '''
        Will go over all events in the input dictionary and only keep those
        that are within the defined survey footprint. Returns a dictionary 
        of event names, sky coordinates and baseline magnitudes. 
        Optionally also returns a dictionary of ROME/REA field numbers with 
        all the event names in each field. [Set romerea=1]
    '''
    import copy
    # Make a deep copy of the known_events_dictionary
    events_in_footprint = copy.deepcopy(known_events_dictionary)
    if romerea==1:
        field_dict = {}
    for event in known_events_dictionary.keys():
        # Check if it is in the survey footprint
        radeg, decdeg = known_events_dictionary[event][0:2]
        val = romerea_check(radeg, decdeg)
        if val == -1:
            events_in_footprint.pop(event, None)
        else:
            #print(event, val)
            if romerea==1:
                if val not in field_dict.keys():
                    field_dict[val]=[event]
                else:
                    field_dict[val].append(event)
    if romerea==1:
        return events_in_footprint, field_dict
    else:
        return events_in_footprint


# Example run
#urls = construct_url_list(survey_dict, years)
#known_events = collect_known_events(urls)
#romerea_events, field_dict = select_events_within_footprint(known_events)

# Write file to output
#with open('romerea_events.json', 'w') as f:
#    json.dump(romerea_events, f)

# How to read a json file
#with open('romerea_events.json', 'r') as f:
#    data = json.load(f)

# How to serve the date in a pandas dataframe
#import pandas as pd
#df = pd.DataFrame(data=romerea_events)
#df = df.fillna(' ').T
#print(df)


#print(field_dict)

