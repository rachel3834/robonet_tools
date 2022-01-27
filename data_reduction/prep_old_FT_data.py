from astropy.io import fits
import numpy as np
import glob
from os import path, rename
from sys import argv
import copy

def update_image_structure(data_dir):

    image_list = glob.glob(path.join(data_dir,'*.fits'))

    for image in image_list:
        hdu = fits.open(image)
        new_image = image.replace('.fits','p.fits')
        new_image = zero_pad_filenames(new_image)

        if len(hdu) == 1 and hdu[0].name == 'PRIMARY':
            sci = copy.copy(hdu[0])
            sci.name = 'SCI'
            sci.header['PIXSCALE'] = sci.header['CCDSCALE']
            bpm = fits.ImageHDU(np.zeros(sci.data.shape))
            bpm.name = 'BPM'
            new_hdu = hdu = fits.HDUList([sci,bpm])
            new_hdu.writeto(new_image, overwrite=True)
            print('Restructured '+new_image)

def zero_pad_filenames(file_path):

    def padd_id(idstr):
        while len(idstr) < 3:
            idstr = '0' + idstr
        return idstr

    file_name = path.basename(file_path)

    elements = file_name.split('_')

    elements[3] = padd_id(elements[3])
    elements[4] = padd_id(elements[4])

    file_name = '_'.join(elements)
    file_path = path.join(path.dirname(file_path), file_name)

    return file_path

def rename_files(data_dir):
    image_list = glob.glob(path.join(data_dir,'*.fits'))

    for image in image_list:
        new_image = zero_pad_filenames(image)
        rename(image, new_image)

if __name__ == '__main__':

    if len(argv) > 1:
        data_dir = argv[1]
    else:
        data_dir = input('Please enter the path to data directory: ')

    update_image_structure(data_dir)
    #rename_files(data_dir)
