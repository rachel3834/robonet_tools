from pyDANDIA import metadata
from sys import argv
from os import path
import glob
import numpy as np

def verify_dataset(red_dir):
    reduction_metadata = metadata.MetaData()
    reduction_metadata.load_all_metadata(metadata_directory=red_dir,
                                         metadata_name='pyDANDIA_metadata.fits')

    data_dir = path.join(red_dir, 'data')
    for image in reduction_metadata.headers_summary[1]['IMAGES']:
        if not path.isfile(path.join(data_dir,image)):
            print('WARNING: Cannot find image from metadata: '+image)
        else:
            print('Found '+image+' from metadata')

    image_list = glob.glob(path.join(data_dir, '*.fits'))
    for image in image_list:
        idx = np.where(reduction_metadata.headers_summary[1]['IMAGES'] == path.basename(image))
        print(idx)
        if len(idx) == 0:
            print('WARNING: Cannot find image '+path.basename(image)+' in metadata')
        else:
            print('Found '+path.basename(image)+' from data directory')
            
if __name__ == '__main__':
    if len(argv) > 1:
        red_dir = argv[1]
    else:
        red_dir = input('Please enter the path to the reduction directory: ')
    verify_dataset(red_dir)
