from os import path, makedirs, remove, rmdir
from sys import argv, exit
import glob
from shutil import move, rmtree
from pyDANDIA import sort_data, automatic_pipeline
import config_utils
import log_utils
import subprocess
from astropy.io import fits
from pyDANDIA import automatic_pipeline
import tarfile

def prepare_data_for_reduction(CONFIG_FILE):

    config = config_utils.get_config(CONFIG_FILE)

    log = log_utils.start_day_log(config, 'data_preparation')

    # Identify newly downloaded data, handling the different formats of
    # compressed imaging data (single fits.fz files) and FLOYDS pipeline output
    # (tar.gz)
    compressed_frames = check_for_new_single_frames(config, log)
    compressed_floyds_frames = check_for_new_floyds_tarballs(config, log)

    if len(compressed_frames) > 0:
        log.info('Preparing and sorting imaging data')

        decompressed_images = decompress_new_single_frames(config, log, compressed_frames)

        #transform_frames(decompressed_images, log)

    else:
        log.info('No imaging data found to process')

    if len(compressed_floyds_frames) > 0:
        log.info('\n Preparing and sorting FLOYDS data')

        decompressed_spectra = decompress_floyds_tarballs(config, compressed_floyds_frames, log)

    # Sort the uncompress/extracted image files in the download directory:
    sort_data.sort_data(config['data_download_dir'], config['separate_instruments'], log=log)

    # Sort the uncompressed SOAR spectra in the download directory:
    sort_soar_data(config, log)

    datasets = get_dataset_list(config, log)

    for dataset_dir in datasets:
        transfer_to_reduction_directory(config, log, dataset_dir)

    log_utils.close_log(log)

def sort_soar_data(config, log):
    """
    Function to sort the SOAR/Goodman spectra into reduction directories
    """
    log.info('Sorting SOAR data')

    # Make a file list of all SOAR/Goodman spectra files
    frame_list = glob.glob(path.join(config['data_download_dir'], 'wecfzst_*SOAR*fits'))
    log.info('Identified ' + str(len(frame_list)) + ' SOAR/Goodman frames')

    # Sort the spectra according to different targets
    for frame in frame_list:

        # Identify the dataset for this frame
        try:
            hdr = fits.getheader(frame)

            ds = sort_data.Dataset()
            ds.target = hdr['OBJECT'].replace('/', '').replace(' ', '-')
            ds.site = hdr['SITEID'].replace('/', '')
            ds.enclosure = 'SOAR'
            ds.tel = hdr['TELESCOP'].replace('/', '')
            ds.instrument = hdr['INSTRUME'].replace('/', '')
            ds.filter = hdr['GRATING'] + '-' + hdr['SLIT']

            ds.id_code = str(ds.target)+'_' + \
                    str(ds.site).lower()+'-'+\
                    str(ds.instrument).lower()+'_'+\
                    ds.filter

            # Transfer the data to the appropriate reduction directory
            red_dir = path.join(config['data_download_dir'], ds.id_code)
            dest_dir = path.join(config['data_download_dir'], ds.id_code, 'data')
            unlocked = automatic_pipeline.check_dataset_dir_unlocked(red_dir, log)

            if not path.isdir(dest_dir):
                makedirs(dest_dir)

            if unlocked:
                move(frame, path.join(dest_dir, path.basename(frame)))
                message = path.basename(frame) + ' --> ' + dest_dir
                log.info(message)
            else:
                message = 'Reduction directory ' + red_dir + ' is locked.  Data will remain in incoming directory'
                log.info(message)

        except OSError:
            if log != None:
                log.info('ERROR opening spectrum ' + frame)

        return ds


def check_for_new_single_frames(config, log):

    if path.isdir(config['data_download_dir']) == False:
        log.info('ERROR: Cannot find data download directory')
        exit()

    new_frames = glob.glob(path.join(config['data_download_dir'],'*fits*'))

    log.info('Found '+str(len(new_frames))+' new imaging frames to process')

    return new_frames

def check_for_new_floyds_tarballs(config, log):
    new_frames = glob.glob(path.join(config['data_download_dir'], '*tar.gz*'))

    log.info('Found ' + str(len(new_frames)) + ' new FLOYDS frames to process')

    return new_frames

def decompress_new_single_frames(config, log, compressed_frames):

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

def decompress_floyds_tarballs(config, compressed_floyds_frames, log):

    # Make a storage directory for processed tarballs, if none already exists
    tarball_dir = path.join(config['data_download_dir'], 'tarballs')
    if path.isdir(tarball_dir) == False:
        makedirs(tarball_dir)
    unused_dir = path.join(config['data_download_dir'], 'unused')
    if path.isdir(unused_dir) == False:
        makedirs(unused_dir)

    decompressed_frames = []

    for tarball in compressed_floyds_frames:
        try:
            with tarfile.open(tarball) as tf:
                # Identify the 1D, extracted and wavelength calibrated FITS spectrum
                # in the list of files in the tarball
                files = [name for name in tf.getnames() if '_ex.fits' in name]
                for spectrum in files:
                    tf.extract(spectrum, config['data_download_dir'])
                    log.info('Extracting reduced spectrum ' + spectrum + ' from ' + tarball)
                    decompressed_frames.append(spectrum)
                tf.close()

            # Move processed tarballs to storage
            move(tarball, path.join(tarball_dir, path.basename(tarball)))
            log.info('Moved processed tarball to storage')

        except:
            log.info('ERROR extracting tarball ' + tarball)
            move(tarball, path.join(unused_dir, path.basename(tarball)))

    return decompressed_frames

def transform_frames(decompressed_frames, log):
    in_use = False

    if in_use:
        for frame in decompressed_frames:
            hdr = fits.getheader(frame)
            if 'WCSERR' not in hdr.keys() or hdr['WCSERR'] != 0:
                hdu = fits.open(frame)
                if hdu != None:
                    hdu[0].data = hdu[0].data[::-1,::-1]
                    for i in range(1,len(hdu),1):
                        if 'BPM' in hdu[i].header['EXTNAME'] or 'ERR' in hdu[i].header['EXTNAME']:
                            hdu[i].data = hdu[i].data[::-1,::-1]
                    hdu.writeto(frame.replace('.fits','_t.fits'), overwrite=True)
                    hdu.close()
                    log.info('-> Transformed frame '+path.basename(frame))
            else:
                log.info('No transform necessary for '+path.basename(frame))

def get_dataset_list(config, log):

    entries = glob.glob(path.join(config['data_download_dir'],'*'))

    datasets = []

    for item in entries:
        if path.isdir(item) and str(path.basename(item)).lower() not in ['unused', 'tarballs', 'tmp']:
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
        unlocked = automatic_pipeline.check_dataset_dir_unlocked(red_dir,log)

        if unlocked:
            frames = glob.glob(path.join(dataset_dir,'data','*fits'))

            log.info('-> Moving '+str(len(frames))+' frames to '+red_dir)

            for f in frames:
                move(f, path.join(red_dir, 'data', path.basename(f)))

        # If the directory is locked for an ongoing reduction, frames can be
        # left over.  These should be moved back to the incoming directory
        # where they will be picked up next time the code runs
        frames_left = glob.glob(path.join(dataset_dir,'data','*fits'))
        log.info(str(len(frames_left))+' frames remain that cannot be moved to reduction directories yet (reductions in progress)')
        log.info('Moving to holding area:')
        if len(frames_left) > 0:
            for f in frames_left:
                dest = path.join(dataset_dir, '..')
                move(f, dest)
                log.info('Returned '+path.basename(f)+' to '+dest)

        rmtree(dataset_dir)

if __name__ == '__main__':

    if len(argv) == 1:
        CONFIG_FILE = input('Please enter the path to the configuration file: ')
    else:
        CONFIG_FILE = argv[1]

    prepare_data_for_reduction(CONFIG_FILE)
