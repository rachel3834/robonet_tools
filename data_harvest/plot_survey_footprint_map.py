from os import path
import argparse
import healpy as hp
from astropy import wcs
from astropy_healpix import HEALPix
from astropy.coordinates import Galactic, TETE, SkyCoord
import numpy as np
from pyDANDIA import plot_all_fields
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
from PIL import Image

# Definition of the healpixel resolution to use, which should correspond to a resolution of
# about 1deg
NSIDE=512
RESOLUTION = hp.nside2resol(NSIDE, arcmin=True) / 60.0 # degrees
print('Working with resolution ' + str(RESOLUTION) + 'deg')
BULGE_REGION = {'ra_range': [260.0, 275.0],
                'dec_range': [-35.0, -25.0]}

# TRILEGAL data maps accessed from Rubin-sim archive:
# https://s3df.slac.stanford.edu/data/rubin/sim-data/
# The highest resolution available has NSIDE=1024, but this map does not have data covering the whole sky and omits
# the Galactic Plane.  So using the NSIDE=512 map as the best available.

# Background maps
BACKGROUND = {'DECam':[ 271.014958, -28.75124167,'18:04:03.59', '-28:45:04.47', 'decam_bulge_image.tiff', 0.263],
              'DSS': [269.335, -29.178, '17:57:20.4', '-29:10:40.8', 'DSS_ROME_footprint_car.png', 6.95],
              'Aladdin': [269.78358, -29.639406, '17:59:08.06', '-29:38:21.86', 'aladdin_rome_footprint.png', 99.999]}



def plot_footprint(args):

    # Load the data on the density of stars as a function of sky location based on
    # the TRILEGAL galactic model
    #(ahp, hp_log_star_density) = create_star_density_map(args)

    # hp.mollview(hp_log_star_density, title="Stellar density",
    #            min=0.0, max=hp_log_star_density.max())
    # hp.graticule()
    # plt.tight_layout()
    # plt.savefig(path.join(args.output_dir, 'stellar_density_map.png'))
    # plt.close(3)

    # Extract from the density map the healpixels for the Bulge region as a 2D image array
    plot_trilegal = False
    if plot_trilegal:
    # Load the data on the density of stars as a function of sky location based on
    # the TRILEGAL galactic model
        (ahp, hp_log_star_density) = create_star_density_map(args)
        bulge_star_density = []
        for ra in np.arange(BULGE_REGION['ra_range'][0], BULGE_REGION['ra_range'][1], RESOLUTION):
            row = []
            for dec in np.arange(BULGE_REGION['dec_range'][0], BULGE_REGION['dec_range'][1], RESOLUTION):
                pixels = calc_hp_healpixels_for_coords(ra, dec, RESOLUTION, NSIDE)
                row.append(np.median(hp_log_star_density[pixels]))
            bulge_star_density.append(row)
        bulge_star_density = np.array(bulge_star_density)

    # Plot a background image
    plot_bkgd = False
    if plot_bkgd:
        bkgd_data = BACKGROUND['Aladdin']
        BKGD_WIDTH = 4.284
        BKGD_HEIGHT = 3.054
        bkgd_image = np.array(Image.open(path.join(args.output_dir, bkgd_data[4])))
        #bkgd_image = np.fliplr(bkgd_image)
        #bkgd_image = np.flipud(bkgd_image)
        dra = (BKGD_WIDTH * 3600.0) / bkgd_image.shape[1]
        x_center = bkgd_image.shape[1]/2.0
        ddec = (BKGD_HEIGHT * 3600.0) / bkgd_image.shape[0]
        y_center = bkgd_image.shape[0]/2.0
        BKGD_HALF_WIDTH = BKGD_WIDTH / 2.0
        BKGD_HALF_HEIGHT = BKGD_HEIGHT / 2.0
        extent = [bkgd_data[0] - BKGD_HALF_WIDTH,
                  bkgd_data[0] + BKGD_HALF_WIDTH,
                  bkgd_data[1] - BKGD_HALF_HEIGHT,
                  bkgd_data[1] + BKGD_HALF_HEIGHT]

    # Create a WCS header for the center of the Bulge region
    #w = create_wcs(bkgd_data[0], bkgd_data[1], dra, ddec, x_center, y_center)

    fig = plt.figure(1, (39, 27))
    ax = plt.subplot()

    if plot_bkgd:
        plt.imshow(bkgd_image, extent=extent)

    plot_rome = True
    if plot_rome:
        for field_id,field_data in plot_all_fields.ROME_FIELDS.items():

            file_name = path.join(args.data_dir, field_id+'_colour.png')

            if path.isfile(file_name):

                image = plt.imread(file_name)
                #image = np.fliplr(image)
                image = np.flipud(image)

                extent = [ field_data[0] - plot_all_fields.FIELD_HALF_WIDTH,
                            field_data[0] + plot_all_fields.FIELD_HALF_WIDTH,
                           field_data[1] - plot_all_fields.FIELD_HALF_HEIGHT,
                            field_data[1] + plot_all_fields.FIELD_HALF_HEIGHT ]

                plt.imshow(image, extent=extent)

                # text = plt.text(field_data[0], field_data[1], field_id,
                #                 fontdict={'fontsize': 14, 'fontweight': 'bold', 'color': 'white'},
                #                           backgroundcolor='none')
                # text.set_path_effects([withStroke(linewidth=3, foreground='black')])

    ax.axis('equal')
    #plot_ranges = plot_all_fields.calc_survey_boundaries()
    #print(plot_ranges)
    #ax.set_xlim([plot_ranges[1], plot_ranges[0]])
    #ax.set_ylim([plot_ranges[2], plot_ranges[3]])
    #plt.grid(color='white', ls='solid')
    (xmin,xmax,ymin,ymax) = plt.axis()
    plt.axis([xmax,xmin,ymin,ymax])
    plt.grid(linestyle='--',color='gray', linewidth=0.5)
    ax.tick_params(axis='x', colors='gray')
    ax.tick_params(axis='y', colors='gray')
    plt.xlabel('RA [deg]', fontsize=60)
    plt.ylabel('Dec [deg]', fontsize=60)
    ax.tick_params(labelsize=60, labelcolor='gray')
    plt.tight_layout()
    plt.savefig(path.join(args.output_dir, 'survey_map.png'))


def create_wcs(ra_center, dec_center, dra, ddec, x_center, y_center):

    # Use the center of the image and Bulge region as the reference pixel
    crval = np.zeros(2)
    crval[0] = ra_center
    crval[1] = dec_center

    crpix = np.zeros(2)
    crpix[0] = y_center
    crpix[1] = x_center

    # Create the WCS object
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = crpix
    w.wcs.cdelt = np.array([dra, ddec])
    w.wcs.crval = crval
    w.wcs.ctype = ["RA---AIT", "DEC--AIT"]

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
    parser.add_argument('data_dir', help='Path to directory with ROME survey field images')
    parser.add_argument('output_dir', help='Path to output directory')
    args = parser.parse_args()

    return args



if __name__ == '__main__':
    args = get_args()
    plot_footprint(args)