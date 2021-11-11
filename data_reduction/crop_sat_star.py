from os import path, mkdir
from sys import argv
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS as aWCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import glob
import shutil

def crop_sat_star_from_images():

    params = get_args()

    params['bkup_dir'] = path.join(params['data_dir'],'fullframe')
    if path.isdir(params['bkup_dir']) == False:
        mkdir(params['bkup_dir'])

    image_list = glob.glob(path.join(params['data_dir'],'*.fits'))

    for image_path in image_list:
        hdul = fits.open(image_path)
        col_min = find_column_bleed(params,hdul)
        row_max = find_sat_row_limit(params,hdul)
        crop_image(params,image_path,hdul,col_min,row_max)

def crop_image(params,image_path,hdul,col_min,row_max):

    data = hdul[0].data[0:row_max, col_min:]

    new_header = update_wcs(hdul[0].header, data.shape[0], data.shape[1])

    hdu_out = fits.PrimaryHDU(data=data, header=new_header)

    bkup_file = path.join(params['bkup_dir'], path.basename(image_path))
    shutil.move(image_path, bkup_file)

    hdu_out.writeto(image_path)
    print('Cropped image '+path.basename(image_path))

def update_wcs(header, new_naxis1, new_naxis2):

    image_wcs = aWCS(header)

    xcentre = int(new_naxis2 / 2.0)
    ycentre = int(new_naxis1 / 2.0)

    centre_world_coords = image_wcs.wcs_pix2world(np.array([[xcentre, ycentre]]), 1)

    pointing = SkyCoord(centre_world_coords[0,0], centre_world_coords[0,1],
                        frame='icrs',unit=(u.deg, u.deg))
    sexigesimal_coords = pointing.to_string(style='hmsdms',sep=':').split()

    header['RA'] = sexigesimal_coords[0]
    header['DEC'] = sexigesimal_coords[1]
    header['CAT-RA'] = sexigesimal_coords[0]
    header['CAT-DEC'] = sexigesimal_coords[1]
    header['CRVAL1'] = float(centre_world_coords[0,0])
    header['CRVAL2'] = float(centre_world_coords[0,1])
    header['CRPIX1'] = xcentre
    header['CRPIX2'] = ycentre

    return header

def find_column_bleed(params,hdul,diagnostics=True):
    image = hdul[0].data

    x_spectrum = image.sum(axis=0)

    if diagnostics:
        plt.plot(np.arange(0,image.shape[0]), x_spectrum, 'r-')
        plt.xlabel('X pixel')
        plt.ylabel('Pixel counts [ADU]')
        plt.savefig(path.join(params['data_dir'],'x_spectrum.png'))
        plt.close()
    xbleed = np.where(x_spectrum == x_spectrum.max())[0][0]
    print('Column bleed at x='+str(xbleed)+' pix')

    # Add a boundary to avoid wider column bleeds:
    xbleed += 5

    return xbleed

def find_sat_row_limit(params,hdul,diagnostics=True):
    image = hdul[0].data

    y_spectrum = image.sum(axis=1)

    # Crop bias ranges off the top and bottom of the frame
    y_spectrum = y_spectrum[30:4000]

    avg = np.median(y_spectrum)
    hilimit = y_spectrum.mean() + 0.5*y_spectrum.std()

    sat_rows = np.where(y_spectrum > hilimit)[0]

    # Centroid on largest clump of saturated rows, then find the lower boundary
    # Weed out non-contiguous rows
    dsat_rows = sat_rows[1:] - sat_rows[0:-1]
    idx = (dsat_rows == 1)
    jdx = np.empty(len(sat_rows), dtype='bool')
    jdx.fill(True)
    jdx[0:len(idx)] = idx
    sat_rows = sat_rows[jdx]
    centroid = int(np.median(sat_rows))

    row_max = centroid
    for r in range(centroid-1,0,-1):
        row = y_spectrum[r]
        if r < centroid and row > hilimit:
            row_min = r
        if r < centroid and row < hilimit:
            break

    #row_min = sat_rows.min()
    print('Lower limit of saturated rows y='+str(row_max)+' pix')

    if diagnostics:
        plt.plot(np.arange(0,len(y_spectrum)), y_spectrum, 'r-')
        plt.plot(np.arange(0,len(y_spectrum)),[avg]*len(y_spectrum),'b-.')
        plt.plot(np.arange(0,len(y_spectrum)),[hilimit]*len(y_spectrum),'b--')
        [xmin, xmax,ymin,ymax] = plt.axis()
        plt.plot([centroid, centroid],[ymin,ymax], 'k-')
        plt.xlabel('X pixel')
        plt.ylabel('Pixel counts [ADU]')
        plt.savefig(path.join(params['data_dir'],'y_spectrum.png'))
        plt.close()

    return row_max

def get_args():
    params = {}
    if len(argv) == 1:
        params['data_dir'] = input('Please enter the path to the data directory: ')
    else:
        params['data_dir'] = argv[1]
    return params

if __name__ == '__main__':
    crop_sat_star_from_images()
