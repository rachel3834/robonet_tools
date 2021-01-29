from os import path
from sys import argv
import numpy as np
from pyDANDIA import logs
from pyDANDIA import metadata
from pyDANDIA import match_utils
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt

def compare_catalog_xmatch_between_reductions():

    params = get_args()

    log = logs.start_stage_log( params['log_dir'], 'compare_crossmatch' )

    meta1 = metadata.MetaData()
    meta1.load_all_metadata(params['red_dir1'], 'pyDANDIA_metadata.fits')
    log.info('Loaded metadata from '+params['red_dir1'])
    meta2 = metadata.MetaData()
    meta2.load_all_metadata(params['red_dir2'], 'pyDANDIA_metadata.fits')
    log.info('Loaded metadata from '+params['red_dir2'])

    # Crossmatch the catalogs by x, y pixel positions
    (meta1_matches, meta2_matches) = crossmatch_pixel_positions(meta1, meta2,
                                                        log, threshold=1.0)

    matched_data = build_matched_arrays(meta1, meta2, meta1_matches,
                                        meta2_matches, log)

    matched_data = calculate_separations_on_sky(matched_data, log)

    # For each matching star, compare RA, Dec and Gaia ID if available
    compare_coordinates(params, matched_data, log)
    compare_gaia_ids(params, matched_data, log)

    logs.close_log(log)

def compare_coordinates(params, matched_data, log):

    fig = plt.figure(1,(10,10))
    plt.rcParams.update({'font.size': 18})

    plt.hist(matched_data['separations'])
    plt.xlabel('Angular separation [deg]')
    plt.ylabel('Frequency')
    plt.title('Separation between world coordinates of matched stars')

    plt.grid()
    plt.savefig(path.join(params['log_dir'],'angular_separations_matched_stars.png'))
    log.info('Plotted the separation between world coordinates of matched stars')

def compare_gaia_ids(params, matched_data, log):

    delta_id = matched_data['gaia_id1'] - matched_data['gaia_id2']
    idx = np.where(delta_id != 0)

    log.info('Comparing identifications of Gaia source IDs')
    log.info('Found '+str(len(idx))+' stars with different identifications')
    log.info('Found '+str(len(matched_data)-len(idx))+' stars with consistent identifications')

    fig = plt.figure(1,(10,10))
    plt.rcParams.update({'font.size': 18})

    plt.hist(delta_id)
    plt.xlabel('Difference in Gaia ID')
    plt.ylabel('Frequency')
    plt.title('Difference in Gaia ID strings from same data release')

    plt.grid()
    plt.savefig(path.join(params['log_dir'],'gaia_source_ids_matched_stars.png'))

def crossmatch_pixel_positions(meta1, meta2, log, threshold=None):
    """Function to crossmatch cartesian pixel positions.
    Returns two indices of corresponding matched stars: the first is
    the index for the first catalog"""

    matched_stars = match_utils.StarMatchIndex()

    x1 = meta1.star_catalog[1]['x']
    y1 = meta1.star_catalog[1]['y']
    x2 = meta2.star_catalog[1]['x']
    y2 = meta2.star_catalog[1]['y']
    nstars1 = len(meta1.star_catalog[1])
    nstars2 = len(meta2.star_catalog[1])
    index1 = np.arange(0,nstars1,1)
    index2 = np.arange(0,nstars2,1)
    log.info('Cross-matching '+str(nstars1)+' stars from catalogue 1 with '+\
                        str(nstars2)+' stars from catalogue 2')

    idx1 = np.full( (nstars2, nstars1), index1 )
    xv1 = np.full( (nstars2, nstars1), x1)
    yv1 = np.full( (nstars2, nstars1), y1)
    idx2 = np.full( (nstars1, nstars2), index2 )
    print('Got here')
    idx2 = np.transpose(idx2)
    xv2 = np.full( (nstars1, nstars2), x2)
    xv2 = np.transpose(xv2)
    yv2 = np.full( (nstars1, nstars2), y2)
    yv2 = np.transpose(yv2)
    print('Got here 2')

    separations = np.sqrt( (xv1-xv2)**2 + (yv1-yv2)**2 )
    print('Got here 3')
    
    if threshold:
        min_seps = np.where( separations <= threshold )
    else:
        min_seps = np.where( separations == separations.min(axis=0))

    meta1_matches = idx1[min_seps]
    meta2_matches = idx2[min_seps]
    log.info('-> Found '+str(len(meta1_matches))+' matching entries')

    return meta1_matches, meta2_matches

def build_matched_arrays(meta1, meta2, meta1_matches, meta2_matches, log):

    data = [ Column(name='x1', data=meta1.star_catalog[1]['x'][meta1_matches]),
             Column(name='y1', data=meta1.star_catalog[1]['y'][meta1_matches]),
             Column(name='ra1', data=meta1.star_catalog[1]['ra'][meta1_matches]),
             Column(name='dec1', data=meta1.star_catalog[1]['dec'][meta1_matches]),
             Column(name='gaia_id1', data=meta1.star_catalog[1]['gaia_source_id'][meta1_matches]),

             Column(name='x2', data=meta2.star_catalog[1]['x'][meta2_matches]),
              Column(name='y2', data=meta2.star_catalog[1]['y'][meta2_matches]),
              Column(name='ra2', data=meta2.star_catalog[1]['ra'][meta2_matches]),
              Column(name='dec2', data=meta2.star_catalog[1]['dec'][meta2_matches]),
              Column(name='gaia_id2', data=meta2.star_catalog[1]['gaia_source_id'][meta2_matches]),

              Column(name='separation', data=np.zeros(len(meta1_matches))) ]

    log.info('Built table of matched star data')

    return Table(data)

def calculate_separations_on_sky(matched_data, log):

    stars1 = SkyCoord(matched_data['ra1'], matched_data['dec1'],
                        frame='icrs', unit=(u.deg, u.deg))

    stars2 = SkyCoord(matched_data['ra2'], matched_data['dec2'],
                        frame='icrs', unit=(u.deg, u.deg))

    matched_data['separations'] = stars1.separation(stars2)

    log.info('Calculated angular separations between crossmatched star sky positions')

    return matched_data

def get_args():

    params = {}
    if len(argv) < 4:
        params['red_dir1'] = input('Please enter the path to the first reduction directory: ')
        params['red_dir2'] = input('Please enter the path to the second reduction directory: ')
        params['log_dir'] = input('Please enter the path to the log directory: ')
    else:
        params['red_dir1'] = argv[1]
        params['red_dir2'] = argv[2]
        params['log_dir'] = argv[3]

    return params

if __name__ == '__main__':
    compare_catalog_xmatch_between_reductions()
