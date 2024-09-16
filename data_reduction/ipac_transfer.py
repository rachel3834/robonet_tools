from os import path, walk
import subprocess
import json
import argparse
import logging
from datetime import datetime

def transfer_field_dataset(args):
    """
    Function to upload to IPAC all files associated with the photometry output for a single ROME field,
    including the source catalog and lightcurve files.
    The transfer uses IPAC's FDT protocol.
    """

    # Initialize the program log file to keep a record of the transfer
    log = start_log(args.log_dir)

    # Sanity checks
    if not path.isdir(args.data_dir):
        raise IOError('Cannot find the input data directory at ' + args.data_dir)

    # IPAC's FDT protocol will be used for the transfer, so load the configuration info
    config = load_fdt_config(args)

    # Loop over all files in the data directory and send them to IPAC
    for (root, dirs, files) in walk(args.data_dir):
        for file in files:
            file_path = path.join(root, file)

            upload_file_to_ipac(args, config, file_path, log)

    close_log(log)

def upload_file_to_ipac(args, config, file_path, log):
    """
    Function to upload a single file to IPAC
    """

    # Build the command to execute
    cmd = (['/bin/java', '-jar', args.fdt_exec, '-c', config['url'],
            '-p', str(config['port']), '-d', 'INBOX', file_path])

    # Run this command synchronously
    log.info('Upload cmd: ' + ' '.join(cmd))
    result = subprocess.run(cmd)

    log.info(result.stdout)

def load_fdt_config(args):
    """
    Function to read the FDT configuration information from a JSON file
    """

    if not path.isfile(args.config):
        raise IOError('Cannot find the FDT configuration file at ' + args.config)

    with open(args.config, 'r') as f:
        config = json.load(f)

    return config

def start_log(log_dir):
    """Function to initialize a log file for the pyDANDIA pipeline.

    The naming convention for the file is [log_name]_<date_string>.log.

    The new file will automatically overwrite any previously-existing logfile
    for the given reduction.

    This function also configures the log file to provide timestamps for
    all entries.

    Parameters:
        log_dir   string        Directory path
                                log_root_name  Name of the log file
    Returns:
        log       open logger object
    """
    log_name = 'ipactransfer'

    # Console output not captured, though code remains for testing purposes
    console = False

    ts = datetime.utcnow()

    if path.isdir(log_dir) == False:
        makedirs(log_dir)

    log_file = path.join(log_dir, log_name + '_' + ts.strftime("%Y-%m-%d") + '.log')

    # To capture the logging stream from the whole script, create
    # a log instance together with a console handler.
    # Set formatting as appropriate.
    log = logging.getLogger(log_name)

    if len(log.handlers) == 0:
        log.setLevel(logging.INFO)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)

        if console == True:
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.INFO)

        formatter = logging.Formatter(fmt='%(asctime)s %(message)s', \
                                      datefmt='%Y-%m-%dT%H:%M:%S')
        file_handler.setFormatter(formatter)

        if console == True:
            console_handler.setFormatter(formatter)

        log.addHandler(file_handler)
        if console == True:
            log.addHandler(console_handler)

    log.info('Started run of ' + log_name + '\n')

    return log


def close_log(log):
    """Function to cleanly shutdown logging functions with a final timestamped
    entry.
    Parameters:
        log     logger Object
    Returns:
        None
    """

    log.info('Processing complete\n')
    logging.shutdown()

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Path to the fields top-level data directory')
    parser.add_argument('log_dir', help='Path to the logging directory')
    parser.add_argument('fdt_exec', help='Path to FDT executable')
    parser.add_argument('field_name', help='Name of the field')
    parser.add_argument('config', help='Path to configuration file for FDT transfer')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    transfer_field_dataset(args)