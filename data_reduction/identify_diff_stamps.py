from os import path
from sys import argv
import glob
from astropy.io import fits

def run(params):

    # Get list of difference image directories:
    diff_dirs = glob.glob(path.join(params['diff_dir'],'*.fits'))

    # Add the image name to the fits header of all products from each
    # diff image directory
    for ddir in diff_dirs:
        identify_products(ddir)

def identify_products(diff_image_dir):

    # Get the image name from the directory name:
    imagename = path.basename(diff_image_dir)
    if '.fits' not in imagename:
        raise IOError('Image name '+imagename+' doesnt have the expected fits suffix')

    # Get a list of the reduction products:
    file_list = glob.glob(path.join(diff_image_dir, '*.fits'))

    # Add the image name to the FITS headers:
    for product_path in file_list:
        hdul = fits.open(product_path)
        new_header = hdul[0].header
        new_header['filename'] = imagename

        new_hdu = fits.PrimaryHDU(hdul[0].data, header=new_header)
        new_hdu.writeto(product_path, overwrite=True)

def get_args():
    params = {}
    if len(argv) == 1:
        params['red_dir'] = input('Please enter the path to the reduction directory: ')
    else:
        params['red_dir'] = argv[1]
    params['diff_dir'] = path.join(params['red_dir'],'diffim')
    return params


if __name__ == '__main__':
    params = get_args()
    run(params)
