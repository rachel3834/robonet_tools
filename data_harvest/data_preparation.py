from os import path, makedirs, remove, rmdir
from sys import argv, exit
import glob
from shutil import move, rmtree
from pyDANDIA import sort_data
import config_utils
import log_utils
import subprocess

CONFIG_FILE = '/data/omega/configs/data_preparation_config.json'
if path.isfile(CONFIG_FILE) == False:
    CONFIG_FILE = path.join(path.expanduser('~'), 'software', 'robonet_tools',
                        'configs', 'data_preparation_config.json')

def prepare_data_for_reduction():

    config = config_utils.get_config(CONFIG_FILE)

    log = log_utils.start_day_log(config, 'data_preparation')

    compressed_frames = check_for_new_frames(config, log)

    if len(compressed_frames) > 0:

        decompressed_frames = decompress_new_frames(config, log, compressed_frames)

        sort_data.sort_data(config['data_download_dir'],config['separate_instruments'],log=log)

        datasets = get_dataset_list(config, log)

        for dataset_dir in datasets:
            transfer_to_reduction_directory(config, log, dataset_dir)

    log_utils.close_log(log)

def check_for_new_frames(config, log):

    if path.isdir(config['data_download_dir']) == False:
        log.info('ERROR: Cannot find data download directory')
        exit()

    new_frames = glob.glob(path.join(config['data_download_dir'],'*fits*'))

    log.info('Found '+str(len(new_frames))+' new frames to process')

    return new_frames

def decompress_new_frames(config, log, compressed_frames):

    decompressed_frames = []

    unused_dir = path.join(config['data_download_dir'], 'unused')
    if path.isdir(unused_dir) == False:
        makedirs(unused_dir)

    for frame in compressed_frames:

        if frame.split('.')[-1] == 'fz':
            args = ['funpack', frame]
            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            p.wait()

            new_frame = frame.replace('.fz','')

            if path.isfile(new_frame):
                remove(frame)
                log.info('Decompressed '+path.basename(frame))

                decompressed_frames.append(new_frame)

            else:
                log.info('ERROR decompressing '+path.basename(frame))

                move(frame, path.join(unused_dir, path.basename(frame)))

        elif frame.split('.')[-1] in ['fits', 'fts']:
            decompressed_frames.append(frame)

        else:
            log.info('Error: Found unrecognized image compression extention: '+path.basename(frame))
            move(frame, path.join(unused_dir, path.basename(frame)))

    return decompressed_frames

def get_dataset_list(config, log):

    entries = glob.glob(path.join(config['data_download_dir'],'*'))

    datasets = []

    for item in entries:
        if path.isdir(item) and 'unused' not in str(item).lower():
            datasets.append(item)

    log.info('Identified '+str(len(datasets))+' datasets to prepare')

    return datasets

def transfer_to_reduction_directory(config, log, dataset_dir):

    dataset_id = path.basename(dataset_dir)

    red_dir = path.join(config['data_reduction_dir'], dataset_id)

    if path.isdir(red_dir) == False:
        log.info('-> Creating a new reduction directory for '+dataset_id)
        move(dataset_dir, red_dir)

    else:
        frames = glob.glob(path.join(dataset_dir,'data','*fits'))

        log.info('-> Moving '+str(len(frames))+' frames to '+red_dir)

        for f in frames:
            move(f, path.join(red_dir, 'data', path.basename(f)))

        frames_left = glob.glob(path.join(dataset_dir,'data','*fits'))
        if len(frames_left) == 0:
            rmtree(dataset_dir)

if __name__ == '__main__':
    prepare_data_for_reduction()
