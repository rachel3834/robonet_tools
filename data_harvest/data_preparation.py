from os import path, makedirs, remove, rmdir
from sys import argv, exit
import glob
from shutil import move, rmtree
from pyDANDIA import sort_data
import config_utils
import log_utils
import subprocess
from astropy.io import fits

def prepare_data_for_reduction(CONFIG_FILE):

    config = config_utils.get_config(CONFIG_FILE)

    log = log_utils.start_day_log(config, 'data_preparation')

    compressed_frames = check_for_new_frames(config, log)

    if len(compressed_frames) > 0:

        decompressed_frames = decompress_new_frames(config, log, compressed_frames)

        transform_frames(decompressed_frames, log)

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

def transform_frames(decompressed_frames, log):

    for frame in decompressed_frames:
        hdr = fits.getheader(frame)
        if 'WCSERR' not in hdr.keys() or hdr['WCSERR'] != 0:
            hdu = fits.open(frame)
            hdu[0].data = hdu[0].data[::-1,::-1]
            for i in range(1,len(hdu),1):
                if 'BPM' in hdu[i].header['EXTNAME'] or 'ERR' in hdu[i].header['EXTNAME']:
                    hdu[i].data = hdu[i].data[::-1,::-1]
            hdu.writeto(frame, overwrite=True)
            hdu.close()
            log.info('-> Transformed frame '+path.basename(frame))

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

    if len(argv) == 1:
        CONFIG_FILE = input('Please enter the path to the configuration file: ')
    else:
        CONFIG_FILE = argv[1]

    prepare_data_for_reduction(CONFIG_FILE)
