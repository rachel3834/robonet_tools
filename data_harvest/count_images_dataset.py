from os import path
from sys import argv
import log_utils
import glob

def count_image_data(params):

    log = log_utils.start_day_log(params,'count_images')

    datasets = glob.glob(path.join(params['top_dir'],'*_???-dom?-1m0-??-f???_?p'))

    nimages = 0
    log.info('# Dataset        N images')
    for subdir in datasets:
        image_list = glob.glob(path.join(subdir, 'data', '*.fits'))
        log.info(path.basename(subdir)+' '+str(len(image_list)))
        nimages += len(image_list)
    log.info('Total number of images: '+str(nimages))

    log_utils.close_log(log)

if __name__ == '__main__':

    params = {}
    if len(argv) == 1:
        params['top_dir'] = input('Please enter the path to the top-level data directory: ')
        params['log_dir'] = input('Please enter the log directory path: ')
    else:
        params['top_dir'] = argv[1]
        params['log_dir'] = argv[2]

    count_image_data(params)
