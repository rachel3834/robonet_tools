from os import path, rename
from sys import argv
import numpy as np
from astropy.io import fits
import glob

def invert_images():

    params = get_args()

    image_list = glob.glob(path.join(params['data_dir'],'*.fits'))

    for image in image_list:
        hdul_in = fits.open(image)
        out_hdul = []
        for i,hdu in enumerate(hdul_in):
            if hdu.is_image:
                data = np.flip(hdu.data, axis=0)
                if i == 0:
                    new_hdu = fits.PrimaryHDU(data=data, header=hdu.header)
                else:
                    new_hdu = fits.ImageHDU(data=data, header=hdu.header)
                out_hdul.append(new_hdu)
            else:
                out_hdul.append(hdu)
        out_hdul = fits.HDUList(out_hdul)

        rename(image, image.replace('.fits', '_orig.fits'))
        output_image = image
        out_hdul.writeto(output_image)
        print('Inverted image '+path.basename(output_image))

def get_args():
    params = {}
    if len(argv) == 1:
        params['data_dir'] = input('Please enter the path to the data directory: ')
    else:
        params['data_dir'] = argv[1]
    return params

if __name__ == '__main__':
    invert_images()
