from os import path
import glob
from astropy.io import fits
import argparse

def reset_imageset_bpms(args):
    """Function to reset the pyDANDIA Bad Pixel Masks of a set of images by
    removing the the relevant image extensions from a set of multi-extension FITS images.
    """

    # Get list of images to process:
    image_list = glob.glob(path.join(args.data_dir, '*.fits'))
    if len(image_list) == 0:
        raise IOError('No images found to process in ' + args.data_dir)

    # For each image, load the HDU list and check for an extension named
    # 'PYDANDIA_PIXEL_MASK'.  If it is present, remove it and overwrite the original
    # file.  If not the output image is a copy of the old one.
    for image_file in image_list:

        with fits.open(image_file) as hdul:
            new_hdul = []
            for extn in hdul:
                if 'PYDANDIA_PIXEL_MASK' not in extn.name:
                    new_hdul.append(extn)
            new_hdul = fits.HDUList(new_hdul)
            new_hdul.writeto(image_file)

        print('Removed pyDANDIA BPM from ' + image_file)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Path to image data directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    reset_imageset_bpms(args)