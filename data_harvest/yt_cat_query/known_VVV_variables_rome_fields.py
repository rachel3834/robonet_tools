import numpy as np
from astropy import units
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord
from tqdm import tqdm
import json

#############################################################
##################### User definitions ######################
# The file containing all Spitzer targets
file_of_all_VVV_variables = '/home/ytsapras/ROME-REA/vvv_vizier_votable.tsv'

# Do you also want a plot with the locations of these targets?
plot_fields_VVV_variables = True

#############################################################
def filter_lines_containing_word(file_path, word, output_file_path):
    """
    Reads a file, filters out lines that do not contain the specified word,
    and writes the filtered lines to a new file.

    Args:
        file_path (str): Path to the input file.
        word (str): Word to search for in each line.
        output_file_path (str): Path to the output file.

    Returns:
        None
    """
    filtered_lines = []
    with open(file_path, 'r') as input_file:
        for line in input_file:
            if word in line:
                filtered_lines.append(line)
    
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(filtered_lines)

def romerea_check(radeg, decdeg):
    '''
    Check if the input coordinates are within the ROME/REA footprint.
    Returns field number (1 to 20) or -1 (if outside footprint). 
    '''
    lhalf = 0.232 #0.220833333333  # 26.5/(120.)
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

#############################################################
def collect_known_VVV_variables(file_of_all_VVV_variables):
    ''' 
        Will collect all unique targets. Returns a dictionary 
        of target names, sky coordinates.
    '''
    # create an empty dictionary
    known_VVV_variables_dictionary = {}
    # open the file
    with open(file_of_all_VVV_variables) as f:
        # loop through each line
        counter = 0
        for line in tqdm(f):
            if counter >112 and counter <200852:
                #print (line)
                # split the line into columns
                columns = line.split(';')
                # extract the name and galactic coords
                name = columns[0]
                classification = columns[1]
                period = columns[3]
                ra = float(columns[4])
                dec = float(columns[5])
                gaiaEDR3 = columns[8]
                coords = SkyCoord(ra * units.degree, dec * units.degree, frame='icrs')
                # convert to decimal degrees using astropy.coordinates
                ra_deg = coords.ra.deg
                dec_deg = coords.dec.deg
                # check if the key already exists in the dictionary
                if name not in known_VVV_variables_dictionary:
                    # store the information in the dictionary
                    known_VVV_variables_dictionary[name] = [ra_deg, dec_deg, classification, period, gaiaEDR3]
            counter = counter + 1
    return known_VVV_variables_dictionary

#############################################################    
def select_VVV_variables_within_footprint(known_VVV_variables_dictionary):
    '''
        Will go over all Spitzer in the input dictionary and only keep those
        that are within the defined survey footprint. Returns a dictionary 
        of names, sky coordinates and ROME/REA field number where
        the target is located.
    '''
    VVV_variables_in_footprint = {}
    
    for alertname, alertinfo in known_VVV_variables_dictionary.items():
        # Check if it is in the survey footprint
        radeg, decdeg, classification, period, gaiaEDR3 = alertinfo[0:]
        romefield = romerea_check(radeg, decdeg)
        if romefield != -1:
            VVV_variables_in_footprint[alertname] = [radeg, decdeg, classification, period, gaiaEDR3, romefield]
            
    return VVV_variables_in_footprint

ROME_FIELDS={'ROME-FIELD-01':[ 267.835895375 , -30.0608178195 , '17:51:20.6149','-30:03:38.9442' ],
            'ROME-FIELD-02':[ 269.636745458 , -27.9782661111 , '17:58:32.8189','-27:58:41.758' ],
            'ROME-FIELD-03':[ 268.000049542 , -28.8195573333 , '17:52:00.0119','-28:49:10.4064' ],
            'ROME-FIELD-04':[ 268.180171708 , -29.27851275 , '17:52:43.2412','-29:16:42.6459' ],
            'ROME-FIELD-05':[ 268.35435 , -30.2578356389 , '17:53:25.044','-30:15:28.2083' ],
            'ROME-FIELD-06':[ 268.356124833 , -29.7729819283 , '17:53:25.47','-29:46:22.7349' ],
            'ROME-FIELD-07':[ 268.529571333 , -28.6937071111 , '17:54:07.0971','-28:41:37.3456' ],
            'ROME-FIELD-08':[ 268.709737083 , -29.1867251944 , '17:54:50.3369','-29:11:12.2107' ],
            'ROME-FIELD-09':[ 268.881108542 , -29.7704673333 , '17:55:31.4661','-29:46:13.6824' ],
            'ROME-FIELD-10':[ 269.048498333 , -28.6440675 , '17:56:11.6396','-28:38:38.643' ],
            'ROME-FIELD-11':[ 269.23883225 , -29.2716684211 , '17:56:57.3197','-29:16:18.0063' ],
            'ROME-FIELD-12':[ 269.39478875 , -30.0992361667 , '17:57:34.7493','-30:05:57.2502' ],
            'ROME-FIELD-13':[ 269.563719375 , -28.4422328996 , '17:58:15.2927','-28:26:32.0384' ],
            'ROME-FIELD-14':[ 269.758843 , -29.1796030365 , '17:59:02.1223','-29:10:46.5709' ],
            'ROME-FIELD-15':[ 269.78359875 , -29.63940425 , '17:59:08.0637','-29:38:21.8553' ],
            'ROME-FIELD-16':[ 270.074981708 , -28.5375585833 , '18:00:17.9956','-28:32:15.2109' ],
            'ROME-FIELD-17':[ 270.81 , -28.0978333333 , '18:03:14.4','-28:05:52.2' ],
            'ROME-FIELD-18':[ 270.290886667 , -27.9986032778 , '18:01:09.8128','-27:59:54.9718' ],
            'ROME-FIELD-19':[ 270.312763708 , -29.0084241944 , '18:01:15.0633','-29:00:30.3271' ],
            'ROME-FIELD-20':[ 270.83674125 , -28.8431573889 , '18:03:20.8179','-28:50:35.3666' ]}

#############################################################

all_VVV_variables = collect_known_VVV_variables(file_of_all_VVV_variables)
VVV_variables_within_rome_footprint = select_VVV_variables_within_footprint(all_VVV_variables)
print(str(len(VVV_variables_within_rome_footprint))+' of '+str(len(all_VVV_variables))+' VVV variables are within the ROME/REA fields.')

# Write the identified alerts file to disk
with open('VVV_variables_in_ROMEREA_2.json', 'w') as f:
    json.dump(VVV_variables_within_rome_footprint, f)

# Open existing json file to read?
#with open('VVV_variables_in_ROMEREA.json', 'r') as f:
#    VVV_variables_within_rome_footprint = json.load(f)

print(len(VVV_variables_within_rome_footprint))

# Plot it?
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text

if plot_fields_VVV_variables == True:
    pixel_scale = 0.389
    naxis1 = 4096
    naxis2 = 4096
    
    FIELD_HALF_WIDTH = (( naxis2 * pixel_scale ) / 3600.0) / 2.0 # Deg
    FIELD_HALF_HEIGHT = (( naxis1 * pixel_scale ) / 3600.0) / 2.0 # Deg
    
    fig = plt.figure(2, edgecolor="k", figsize=[12,12])
    ax = fig.add_subplot(1, 1, 1)
    
    for field_id,field_data in ROME_FIELDS.items():
        ar = Rectangle((field_data[0] - FIELD_HALF_WIDTH, field_data[1] - FIELD_HALF_HEIGHT), 
                    width=2*FIELD_HALF_WIDTH, height=2*FIELD_HALF_HEIGHT, angle=0.0, 
                    fill=True, color='red',ec='blue', ls='-', lw=2.0, alpha=0.2)
        txt = Text(field_data[0],field_data[1],field_id,fontsize=8,ha='center',va='center')
        ax.add_patch(ar)
        ax.add_artist(txt)
    
    ras, decs, classification, period, gaiaEDR3, fieldnr = zip(*VVV_variables_within_rome_footprint.values())
    ras = [v[0] for v in VVV_variables_within_rome_footprint.values()]
    decs = [v[1] for v in VVV_variables_within_rome_footprint.values()]
    
    ax.scatter(np.array(ras), np.array(decs), marker='x', color='black', s=2, alpha=0.6)    
    
    ax.grid(True, color='cyan', ls='solid')
    ax.set_aspect('equal')
    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('Dec [deg]')
    ax.set_xlim(271.2,267.5)
    ax.set_ylim(-30.7,-27.6)
    nevs = len(ras)
    plt.text(270.6, -30.52, 'N=%s' % (nevs), bbox=dict(fill=False, edgecolor='blue', linewidth=1))
    plt.show()
