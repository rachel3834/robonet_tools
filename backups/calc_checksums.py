from os import path
from sys import argv
import glob
import hashlib

def checksum_phot_products(red_dir, log_path):

    if path.isdir(red_dir) == False:
        raise IOError('Cannot find input photometry directory')

    log = open(log_path, 'w')

    data_products = glob.glob(path.join(red_dir, '*'))

    for product in data_products:
        if not path.isfile(product):
            dataset = path.basename(product)

            dataset_products = glob.glob(path.join(product, '*'))

            for dproduct in dataset_products:
                file_name = path.basename(dproduct)

                sum = calc_blake2_checksum(dproduct)

                log.write(dataset+'  '+file_name+'  '+str(sum)+'\n')
                
    log.close()

def calc_blake2_checksum(file_path,log):
    """Function to calculate the checksum of a file using the BLAKE2b algorithm.
    Note that this algorithm must read data in binary bytes format.
    The read function appears to be able to handle regular file types and HDF5
    this way but not SQLITE3 DB, so this function necessarily skips that file"""

    with open(file_path,'rb') as f:
        file_data = f.read()

    sum = hashlib.blake2b(file_data).hexdigest()

    print(file_path +'  '+str(sum))

    f.close()

    return sum

if __name__ == '__main__':

    if len(argv) == 1:
        red_dir = input('Please enter the path to the top-level photometry archive directory: ')
        log_path = input('Please enter the path to the log file: ')
    else:
        red_dir = argv[1]
        log_path = argv[2]

    checksum_phot_products(red_dir, log_path)
