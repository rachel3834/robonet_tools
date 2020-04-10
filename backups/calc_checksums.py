from os import path
from sys import argv
import subprocess
import glob


def checksum_phot_products(red_dir):

    if path.isdir(red_dir) == False:
        raise IOError('Cannot find input photometry directory')

    data_products = glob.glob(path.join(red_dir, '*'))
    print(data_products)

    for product in data_products:
        print(product)
        if path.isfile(product):
            calc_checksum(product)
        else:
            dataset_products = glob.glob(path.join(product, '*'))
            print(dataset_products)
            for dproduct in dataset_products:
                calc_checksum(dproduct)

def calc_checksum(file_path):

    cl = ['md5sum', file_path]

    p = subprocess.Popen(cl)
    stdoutput = p.communicate()[0]

    if stdoutput:
        print(stdoutput)

if __name__ == '__main__':

    if len(argv) == 1:
        red_dir = input('Please enter the path to the top-level photometry archive directory: ')
    else:
        red_dir = argv[1]

    checksum_phot_products(red_dir)
    
