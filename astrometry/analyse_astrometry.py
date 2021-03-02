from sys import argv
from os import path
from astropy.coordinates import SkyCoord
from astropy import units as u
from pyDANDIA import crossmatch
from pyDANDIA import metadata
import numpy as np
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

def run_analysis():

    (option, input_file, output_dir) = get_args()

    if option == 'C':
        xmatch = read_crossmatch_data(input_file)

        (separations, masked_ra, masked_dec) = calc_astrometric_residuals_xmatch(xmatch)

    elif option == 'M':
        meta = metadata.MetaData()
        meta.load_a_layer_from_file(path.dirname(input_file),
                                    path.basename(input_file),
                                    'star_catalog')

        (separations, masked_ra, masked_dec) = calc_astrometric_residuals_metadata(meta)

    plot_angular_separations(output_dir, separations, masked_ra, masked_dec)

def plot_angular_separations(output_dir, separations, masked_ra, masked_dec):
    """Based on code from Markus Hundertmark"""

    nxbins = 50
    nybins = 50
    print('Plotting stars between RA='+str(masked_ra.min())+' - '+str(masked_ra.max())+\
            ' and Dec='+str(masked_dec.min())+' - '+str(masked_dec.max()))

    binned_stat = binned_statistic_2d(masked_ra, masked_dec,
                                      separations,
                                      statistic='median',
                                      bins = [nxbins,nybins],
                                      range=[[masked_ra.min(), masked_ra.max()],
                                            [masked_dec.min(), masked_dec.max()]])
    fig, ax1 = plt.subplots(1,1)
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9)
    im = ax1.imshow(binned_stat.statistic.T,
                    cmap = 'gist_rainbow', origin='bottom',
                    extent=(0, nxbins, 0, nybins),
                    vmin = 0.0, vmax = separations.max())

    nticks = 5
    bin_incr = float(nxbins)/float(nticks)
    xincr = (masked_ra.max() - masked_ra.min())/nticks
    yincr = (masked_dec.max() - masked_dec.min())/nticks
    xticks = []
    xlabels = []
    yticks = []
    ylabels = []
    for i in range(0,nticks,1):
        xticks.append(i*bin_incr)
        xlabels.append(str(np.round((masked_ra.min()+i*xincr),3)))
        yticks.append(i*bin_incr)
        ylabels.append(str(np.round((masked_dec.min()+i*yincr),3)))

    plt.xticks(xticks, xlabels, rotation=45.0)
    plt.yticks(yticks, ylabels, rotation=45.0)
    #ax1.set_xticks(binned_stat.x_edge,3)
    #ax1.set_xticklabels(np.round(binned_stat.x_edge,3))
    #ax1.set_yticklabels(np.round(binned_stat.y_edge,3))

    cb = fig.colorbar(im, ax = ax1, label = 'Distance WCS - Gaia match [deg]')
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.savefig(path.join(output_dir,'ang_separations.png'))

def calc_astrometric_residuals_xmatch(xmatch):

    mask = np.where(xmatch.stars['gaia_ra'] > 0.0)
    gaia_positions = SkyCoord(xmatch.stars['gaia_ra'][mask], xmatch.stars['gaia_dec'][mask],
                                frame='icrs', unit=(u.deg,u.deg))
    rome_positions = SkyCoord(xmatch.stars['ra'][mask], xmatch.stars['dec'][mask],
                                frame='icrs', unit=(u.deg,u.deg))

    separations = calc_astrometric_residuals_spherical(gaia_positions, rome_positions, mask)

    return separations, xmatch.stars['ra'][mask].data, xmatch.stars['dec'][mask].data

def calc_astrometric_residuals_metadata(meta):

    mask = np.where(meta.star_catalog[1]['gaia_ra'].data > 0.0)
    gaia_positions = SkyCoord(meta.star_catalog[1]['gaia_ra'][mask], meta.star_catalog[1]['gaia_dec'][mask],
                                frame='icrs', unit=(u.deg,u.deg))
    rome_positions = SkyCoord(meta.star_catalog[1]['ra'][mask], meta.star_catalog[1]['dec'][mask],
                                frame='icrs', unit=(u.deg,u.deg))

    separations = calc_astrometric_residuals_spherical(gaia_positions, rome_positions, mask)

    return separations, meta.star_catalog[1]['ra'][mask].data, meta.star_catalog[1]['dec'][mask].data

def calc_astrometric_residuals_spherical(gaia_positions, rome_positions, mask,
                                        small_angles=True):

    if not small_angles:
        a1 = (90.0 - rome_positions.dec.value) * (np.pi/180.0)
        a2 = (90.0 - gaia_positions.dec.value) * (np.pi/180.0)
        a3 = (rome_positions.ra.value - gaia_positions.ra.value) * (np.pi/180.0)
        cos_gamma = np.cos(a1)*np.cos(a2) + np.sin(a1)*np.sin(a2)*np.cos(a3)

        idx = np.where(cos_gamma-1.0 < 1e-15)[0]
        cos_gamma[idx] = 1.0

        separations = np.arccos(cos_gamma) * (180.0/np.pi)

    else:
        delta_ra = rome_positions.ra.value - gaia_positions.ra.value
        delta_dec = rome_positions.dec.value - gaia_positions.dec.value
        separations = np.sqrt( (delta_ra*np.cos(rome_positions.dec.value))**2 +
                                delta_dec*delta_dec )

    return separations

def calc_astrometric_residuals_cartesian(gaia_positions, rome_positions, mask):

    separations = np.sqrt( (xmatch.stars['ra'][mask] - xmatch.stars['gaia_ra'][mask])**2 +
                        (xmatch.stars['dec'][mask] - xmatch.stars['gaia_dec'][mask])**2 )

    return separations

def read_crossmatch_data(input_file):

    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(input_file)

    return xmatch

def get_args():
    if len(argv) == 1:
        option = input('Use [M]etadata or [C]rossmatch file for calculations? ')
        option = str(option).upper()
        if option.upper() == 'C':
            input_file = input('Please enter the path to the crossmatch file: ')
        elif option.upper() == 'M':
            input_file = input('Please enter the path to the metadata file: ')

    else:
        option = argv[1]
        option = str(option).upper()
        input_file = argv[2]

    output_dir = path.dirname(input_file)

    return option, input_file, output_dir

if __name__ == '__main__':
    if len(argv) > 1:
        option = argv[1]
    else:
        option = 'Use metadata '
    run_analysis()
