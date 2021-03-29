from os import path, stat, mkdir
from sys import argv
from shutil import move
import glob
import log_utils

def search_invalid_data(params):

    log = log_utils.start_day_log(params,'data_validator')

    red_dir_list = glob.glob(path.join(params['top_level_dir'],'*'))
    exclude_dirs = ['logs', 'config', 'downloads']

    for red_dir in red_dir_list:
        if path.isdir(red_dir) and red_dir not in exclude_dirs:
            review_dataset(red_dir, log)

    log_utils.close_log(log)

def review_dataset(red_dir,log):

    log.info('Reviewing data for '+red_dir)

    min_image_size_bytes = 100000000
    image_list = glob.glob(path.join(red_dir, 'data', '*.fits'))
    unused_dir = path.join(red_dir, 'data', 'unused')

    for image in image_list:
        status = stat(image)
        if status.st_size < min_image_size_bytes:
            if not path.isdir(unused_dir):
                mkdir(unused_dir)
            dest = path.join(unused_dir, path.basename(image))
            move(image, dest)
            log.info('Moving '+path.basename(image)+' from '+\
                    path.basename(red_dir)+' to storage, file size='+\
                    str(status.st_size))

if __name__ == '__main__':
    params = {}
    if len(argv) == 1:
        params['top_level_dir'] = input('Please enter the path to the top-level data reduction area: ')
    else:
        params['top_level_dir'] = argv[1]

    params['log_dir'] = path.join(params['top_level_dir'], 'logs')

    search_invalid_data(params)
