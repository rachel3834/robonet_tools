from os import path
from sys import argv
import subprocess
import glob


def checksum_phot_products(red_dir, log_path):

    if path.isdir(red_dir) == False:
        raise IOError('Cannot find input photometry directory')

    log = open(log_path, 'w')

    data_products = glob.glob(path.join(red_dir, '*'))

    for product in data_products:
        if path.isfile(product):
            calc_checksum(product, log)
        else:
            dataset_products = glob.glob(path.join(product, '*'))

            for dproduct in dataset_products:
                calc_checksum(dproduct, log)

    log.close()

def calc_checksum(file_path, log):

    cl = ['md5sum', file_path]

    p = subprocess.Popen(cl)
    stdoutput = p.communicate()[0]

    if stdoutput:
        log.write(stdoutput+'\n')
        import pdb; pdb.set_trace()
    else:
        print(file_path)

if __name__ == '__main__':

    if len(argv) == 1:
        red_dir = input('Please enter the path to the top-level photometry archive directory: ')
        log_path = input('Please enter the path to the log file: ')
    else:
        red_dir = argv[1]
        log_path = argv[2]

    checksum_phot_products(red_dir, log_path)