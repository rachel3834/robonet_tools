from os import path
from sys import argv
import glob
import hashlib
import calc_checksums

def verify_checksum_phot_products(red_dir, log_path):

    original_checksums = read_checksum_log(log_path)

    for dataset,product in original_checksums.items():

        local_path = path.join(red_dir, dataset, product[0])

        if path.isfile(local_path):
            local_sum = calc_checksums.calc_blake2_checksum(local_path)

            if str(local_sum) != str(product[1]):
                print('WARNING: '+dataset+' '+product[0]+' fails checksum test')
            else:
                print(dataset+' '+product[0]+' OK ')

        else:
            print('WARNING: Cannot find expected data product '+dataset+' '+product[0])

def read_checksum_log(log_path):

    if path.isfile(log_path) == False:
        raise IOError('Cannot find checksum log file '+log_path)

    with open(log_path,'r') as f:
        file_lines = f.readlines()

    original_checksums = {}

    for line in file_lines:
        entries = line.replace('\n','').split()
        original_checksums[entries[0]] = [ entries[1], entries[2] ]

    f.close()

    return original_checksums

if __name__ == '__main__':

    if len(argv) == 1:
        red_dir = input('Please enter the path to the top-level photometry archive directory: ')
        log_path = input('Please enter the path to the checksum log file: ')
    else:
        red_dir = argv[1]
        log_path = argv[2]

    verify_checksum_phot_products(red_dir, log_path)
