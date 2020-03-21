from os import path, getcwd, remove
from sys import argv
import glob
import compress_data_products
import upload_aws
import logging

def process_datasets():

    log = start_log()

    (datasets,local_root,aws_root) = get_args(log)

    for data_set in datasets:
        log.info('\nProcessing dataset '+data_set)

        compress_data_for_all_cameras(data_set,log)

        upload_aws.upload_directory(data_set, local_root, aws_root)

        log.info('Completed upload of '+data_set)
        
    close_log(log)

def get_args(log):

    if len(argv) == 1:
        file_path = input('Please enter the path to the list of dataset reduction paths:')
        local_root = input('Please enter the local root path (will be stripped off the path uploaded to AWS): ')
        aws_root = input('Please enter the AWS root path (will be prefixed to the path uploaded to AWS): ')
    else:
        file_path = argv[1]
        local_root = argv[2]
        aws_root = argv[3]

    if not path.isfile(file_path):
        raise IOError('Cannot find input list of datasets to process')

    datasets = []
    file_lines = open(file_path,'r').readlines()
    log.info('Read list of '+str(len(file_lines))+' datasets to process:')

    for l in file_lines:
        datasets.append(l.replace('\n',''))
        log.info(datasets[-1])

    return datasets, local_root, aws_root

def compress_data_for_all_cameras(red_dir,log):
    """
    Function to run the compression process for all datasets for each different
    camera within a reduction directory.

    DanDIA reductions are typically organized into two subdirectories:
    11/ containing the reductions for Sinistro cameras and
    22/ containing the reductions for SBIG cameras.

    The compression function loops over all camera subdirectories within this tree.
    """

    for camera_class in ['11', '22']:
        camera_dir = path.join(red_dir, camera_class)

        if path.isdir(camera_dir):
            log.info(' -> Found data for camera class '+camera_class)

            camera_dir_list = glob.glob(path.join(camera_dir,'*'))

            for event_dir in camera_dir_list:
                log.info(' -> Processing '+event_dir)

                compress_data_products.compress_dandia_reduced_data_products(event_dir)

            log.info(' -> Completed compression of '+camera_class+' datasets')

    log.info(' -> Completed compression of all data for '+red_dir)

def start_log():

    # Console output not captured, though code remains for testing purposes
    console = False

    log_file = path.join(getcwd(), 'datasets_compressed_uploaded.log')
    if path.isfile(log_file) == True:
        remove(log_file)

    # To capture the logging stream from the whole script, create
    # a log instance together with a console handler.
    # Set formatting as appropriate.
    log = logging.getLogger( 'datasets_compressed_uploaded' )

    if len(log.handlers) == 0:
        log.setLevel( logging.INFO )
        file_handler = logging.FileHandler( log_file )
        file_handler.setLevel( logging.INFO )

        if console == True:
            console_handler = logging.StreamHandler()
            console_handler.setLevel( logging.INFO )

        formatter = logging.Formatter( fmt='%(asctime)s %(message)s', \
                                    datefmt='%Y-%m-%dT%H:%M:%S' )
        file_handler.setFormatter( formatter )

        if console == True:
            console_handler.setFormatter( formatter )

        log.addHandler( file_handler )
        if console == True:
            log.addHandler( console_handler )

    log.info( 'Started run')

    return log

def close_log(log):
    log.info( 'Processing complete\n' )
    logging.shutdown()

if __name__ == '__main__':
    process_datasets()
