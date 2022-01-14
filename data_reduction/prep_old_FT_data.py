from astropy.io import fits
import numpy as np
import glob
from os import path
from sys import argv
import copy

def update_image_structure(data_dir):

    image_list = glob.glob(path.join(data_dir,'*.fits'))

    for image in image_list:
        hdu = fits.open(image)
        new_image = image.replace('.fits','p.fits')

        if len(hdu) == 1 and hdu[0].name == 'PRIMARY':
            sci = copy.copy(hdu[0])
            sci.name = 'SCI'
            sci.header['PIXSCALE'] = sci.header['CCDSCALE']
            bpm = fits.ImageHDU(np.zeros(sci.data.shape))
            bpm.name = 'BPM'
            new_hdu = hdu = fits.HDUList([sci,bpm])
            new_hdu.writeto(new_image, overwrite=True)
            print('Restructured '+new_image)

if __name__ == '__main__':

    if len(argv) > 1:
        data_dir = argv[1]
    else:
        data_dir = input('Please enter the path to data directory: ')

    update_image_structure(data_dir)
