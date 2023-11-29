from os import path
import argparse
import healpy as hp
from astropy import wcs
from astropy_healpix import HEALPix
from astropy.coordinates import Galactic, TETE, SkyCoord
import numpy as np
from pyDANDIA import plot_all_fields
import matplotlib.pyplot as plt

# Definition of the healpixel resolution to use, which should correspond to a resolution of
# about 1deg
NSIDE=256
RESOLUTION = hp.nside2resol(NSIDE, arcmin=True) / 60.0 # degrees
print('Working with resolution ' + str(RESOLUTION) + 'deg')
BULGE_REGION = {'ra_range': [260.0, 275.0],
                'dec_range': [-35.0, -25.0]}

# TRILEGAL data maps accessed from Rubin archive:
# https://s3df.slac.stanford.edu/data/rubin/sim-data/

def plot_footprint(args):

    # Load the data on the density of stars as a function of sky location based on
    # the TRILEGAL galactic model
    (ahp, hp_log_star_density) = create_star_density_map(args)

    # Extract from the density map the healpixels for the Bulge region as a 2D image array
    bulge_star_density = []
    for ra in np.arange(BULGE_REGION['ra_range'][0], BULGE_REGION['ra_range'][1], RESOLUTION):
        row = []
        for dec in np.arange(BULGE_REGION['dec_range'][0], BULGE_REGION['dec_range'][1], RESOLUTION):
            pixels = calc_hp_healpixels_for_coords(ra, dec, RESOLUTION, NSIDE)
            row.append(np.median(hp_log_star_density[pixels]))
        bulge_star_density.append(row)
    bulge_star_density = np.array(bulge_star_density)

    # Create a WCS header for the center of the Bulge region
    w = create_bulge_wcs(BULGE_REGION, RESOLUTION, bulge_star_density)

    plt.subplot(projection=w)
    plt.imshow(bulge_star_density)
    plt.grid(color='white', ls='solid')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.show()

    #hp.mollview(hp_log_star_density, title="Stellar density",
    #            min=0.0, max=hp_log_star_density.max())
    #hp.graticule()
    #plt.tight_layout()
    #plt.savefig(path.join(args.output_dir, 'stellar_density_map.png'))
    #plt.close(3)

def create_bulge_wcs(BULGE_REGION, RESOLUTION, image_data):

    # Use the center of the image and Bulge region as the reference pixel
    crval = np.zeros(2)
    crval[0] = BULGE_REGION['ra_range'][0] + (BULGE_REGION['ra_range'][1] - BULGE_REGION['ra_range'][0])/2.0
    crval[1] = BULGE_REGION['dec_range'][0] + (BULGE_REGION['dec_range'][1] - BULGE_REGION['dec_range'][0])/2.0

    crpix = np.zeros(2)
    crpix[0] = image_data.shape[1]/2.0
    crpix[1] = image_data.shape[0]/2.0

    # Create the WCS object
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = crpix
    w.wcs.cdelt = np.array([RESOLUTION, RESOLUTION])
    w.wcs.crval = crval
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    return w

def calc_hp_healpixels_for_coords(ra, dec, resolution, nside):
    """For a given set of coordinates, with RA, Dec in decimal degrees, this function calculates the
    healpixels overlapping those coordinates for a map with the given NSIDE"""

    phi = np.deg2rad(ra)
    theta = (np.pi/2.0) - np.deg2rad(dec)
    radius = np.deg2rad(resolution)
    xyz = hp.ang2vec(theta, phi)
    pixels = hp.query_disc(nside, xyz, radius)

    return pixels

def create_star_density_map(args):
    """Function to create a HEALPix object based on a map of Milky Way stellar
    density generated from the Trilegal galactic model.  These data are in
    Galactic coordinates, so this needs to be rotated in order to map it
    to healpix"""

    star_density_map = load_star_density_data(args, limiting_mag=24.7)
    hp_star_density = rotateHealpix(star_density_map)
    hp_log_star_density = np.zeros(len(hp_star_density))
    idx = hp_star_density > 0.0
    hp_log_star_density[idx] = np.log10(hp_star_density[idx])

    ahp = HEALPix(nside=NSIDE, order='ring', frame=TETE())

    return ahp, hp_log_star_density

def load_star_density_data(args, limiting_mag=28.0):

    data_file = path.join(args.star_map_file)
    if path.isfile(data_file):
        with np.load(data_file) as npz_file:
            star_map = npz_file['starDensity']
            mag_bins = npz_file['bins']

            dmag = abs(mag_bins - limiting_mag)
            idx = np.where(dmag == dmag.min())[0]

            star_density_map = np.copy(star_map[:,idx]).flatten()
            star_density_map = hp.reorder(star_density_map, n2r=True)

        return star_density_map

    else:
        raise IOError('Cannot find star density map data file at '+data_file)

    return None

def rotateHealpix(hpmap, transf=['C','G'], phideg=0., thetadeg=0.):
    """Rotates healpix map from one system to the other. Returns reordered healpy map.
    Healpy coord transformations are used, or you can specify your own angles in degrees.
    To specify your own angles, ensure that transf has length != 2.
    Original code by Xiaolong Li
    """

    # For reasons I don't understand, entering in ['C', 'G'] seems to do the
    # transformation FROM galactic TO equatorial. Possibly something buried in
    # the conventions used by healpy.

    # Heavily influenced by stack overflow solution here:
    # https://stackoverflow.com/questions/24636372/apply-rotation-to-healpix-map-in-healpy

    nside = hp.npix2nside(len(hpmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

    # Define a rotator
    if len(transf) == 2:
        r = hp.Rotator(coord=transf)
    else:
        r = hp.Rotator(deg=True, rot=[phideg,thetadeg])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hpmap, trot, prot)

    return rot_map

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('star_map_file', help='Path to star density map file, e.g. /path/TRIstarDensity_r_nside_64.npz')
    parser.add_argument('output_dir', help='Path to output directory')
    args = parser.parse_args()

    return args



if __name__ == '__main__':
    args = get_args()
    plot_footprint(args)