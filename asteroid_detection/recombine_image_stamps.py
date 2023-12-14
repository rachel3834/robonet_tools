from os import path
from pyDANDIA import metadata
from astropy.io import fits
from astropy.table import Table, Column
import numpy as np
import argparse

def stamps_to_fullframe_image(stamps_table, stamp_dir, filename):
    """
    Function to take a set of image stamps in FITS format, as output by pyDANDIA, and recombine them to
    produce a fullframe image.

    Parameters:
        stamps_table    Table       Dimensions of stamp images, extracted from a pyDANDIA metadata object  stamps table
        stamp_dir       string      Path to directory containing a set image stamps in FITS format
        filename        string      Prefix of the stamp filenames, which are distinguished with numerical indices
                                    e.g. 'diff_stamp', without the .fits extension
        output_image    string      Path to the output fullframe image

    Returns:
        data            array       Fullframe image pixel data
    """

    # Create a holding array for the fullframe image data, sized according to the maximum
    # image pixel boundaries given in the stamps table
    data = np.zeros((stamps_table['Y_MAX'].max(), stamps_table['X_MAX'].max()))

    # Looping over all of the stamps, read in each stamp in turn, combining the pixel data
    for i in meta.stamps[1]['PIXEL_INDEX']:
        stamp_image = path.join(stamp_dir, filename+'_'+str(i)+'.fits')
        stamp_data = fits.getdata(stamp_image)
        ymin = stamps_table['Y_MIN'][i]
        ymax = stamps_table['Y_MAX'][i]
        xmin = stamps_table['X_MIN'][i]
        xmax = stamps_table['X_MAX'][i]
        data[ymin:ymax, xmin:xmax] = stamp_data

    return data

def parse_stamps_table(meta):
    """
    Function to extract the stamps table from a pyDANDIA metadata file in the form of an Astropy table
    with integer values rather than the default strings

    Parameters:
        meta    Metadata object     pyDANDIA metadata file contents, including the stamps table

    Returns:
        stamps  Table object        Extracted table of stamp dimensions
    """

    col_list = [
        Column(name='PIXEL_INDEX', data=meta.stamps[1]['PIXEL_INDEX'], dtype=int),
        Column(name='Y_MIN', data=meta.stamps[1]['Y_MIN'], dtype=int),
        Column(name='Y_MAX', data=meta.stamps[1]['Y_MAX'], dtype=int),
        Column(name='X_MIN', data=meta.stamps[1]['X_MIN'], dtype=int),
        Column(name='X_MAX', data=meta.stamps[1]['X_MAX'], dtype=int)
    ]

    return Table(col_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('metafile', help='Path to pyDANDIA meta data file')
    parser.add_argument('stamp_dir', help='Path to a directory of stamp images')
    parser.add_argument('filename', help='Prefix of stamp images to be combined')
    parser.add_argument('output_image', help='Path of output fullframe images')
    args = parser.parse_args()

    meta = metadata.MetaData()
    meta.load_all_metadata(metadata_directory=path.dirname(args.metafile),
                           metadata_name=path.basename(args.metafile))
    stamps = parse_stamps_table(meta)

    data = stamps_to_fullframe_image(stamps, args.stamp_dir, args.filename)

    fits.writeto(args.output_image, data, overwrite=True)
