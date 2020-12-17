import numpy as np
from os import path
from sys import argv

def select_stars_by_photometry():

    params = get_args()

    photometry = read_colour_photometry_file(params)

    star_index = select_stars(params, photometry)

    output_image_regions(params,photometry,star_index)

def output_image_regions(params,photometry,star_index, colour='green'):

    f = open(params['region_file'], 'w')

    for j in star_index:
        f.write('circle '+str(photometry['x'][j])+' '+str(photometry['y'][j])+
                ' # color='+colour+'\n')

    f.close()
    
def select_stars(params, photometry):

    xdata = photometry[params['xcol']]
    ydata = photometry[params['ycol']]

    xdx1 = np.where(xdata >= params['xcol_min'])
    xdx2 = np.where(xdata <= params['xcol_max'])
    ydx1 = np.where(ydata >= params['ycol_min'])
    ydx2 = np.where(ydata <= params['ycol_max'])
    idx = set(xdx1).intersection(set(xdx2))
    idx = idx.intersection(set(ydx1))
    star_index = list(idx.intersection(set(ydx2)))

    return star_index

def read_colour_photometry_file(params):

    if not path.isfile(params['phot_file']):
        raise IOError('Cannot find photometry file: '+params['phot_file'])

    file_data = open(params['phot_file']).readlines()
    data = []
    for line in file_data:
        if '#' not in line and len(line.replace('\n',''))>0:
            items = line.replace('\n','').split()
            entry = []
            for item in items:
                entry.append(float(item))
            data.append(entry)
        elif '#' in line:
            headers = line.replace('\n','').replace('#','').split()

    data = np.array(data)

    table_data = []

    for i, name in headers:
        table_data.append(table.Column(name=name, data=data[:,i]))

    photometry = table.Table(data=table_data)

    return photometry

def get_args():
    params = {}
    if len(argv) > 1:
        params['phot_file'] = argv[1]
        params['region_file'] = argv[2]
        params['xcol'] = argv[3]
        params['xcol_min'] = argv[4]
        params['xcol_max'] = argv[5]
        params['ycol'] = argv[6]
        params['ycol_min'] = argv[7]
        params['ycol_max'] = argv[8]
    else:
        params['phot_file'] = input('Please enter the path to the colour photometry file: ')
        params['region_file'] = input('Please enter the path to the output region file: ')
        params['xcol'] = input('Please enter the data column to use for selection in the x-axis: ')
        params['xcol_min'] = input('Please enter the minimum acceptable value for the x-value: ')
        params['xcol_max'] = input('Please enter the maximum acceptable value for the x-value: ')
        params['ycol'] = input('Please enter the data column to use for selection in the y-axis: ')
        params['ycol_min'] = input('Please enter the minimum acceptable value for the y-value: ')
        params['ycol_max'] = input('Please enter the maximum acceptable value for the y-value: ')
    return params
