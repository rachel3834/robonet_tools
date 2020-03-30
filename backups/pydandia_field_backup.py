import os
import sys
import glob
import subprocess
import upload_aws

def backup_field_photometry_products(params):
    """Function to gather all the photometric data products for all datasets
    for a given field and back them up to another disk"""

    if not os.path.isdir(params['output_dir']):
        raise IOError('Cannot find output directory '+params['output_dir'])

    if os.path.isfile(params['db_path']):
        rsync_file(params['db_path'], os.path.join(params['output_dir'], os.path.basename(params['db_path'])))
    else:
        print('WARNING: Could not find field database at '+params['db_path'])

    red_dirs = glob.glob(params['dir_path'],params['field_name']+'*')
    print('Found '+str(len(red_dirs))+' reduction directories')

    for dir in red_dirs:
        dset = os.path.basename(dir)
        phot_source_file = os.path.join(dir,'photometry.hd5')
        phot_dest_file = os.path.join(params['output_dir'], dset+'_photometry.hd5')
        if os.path.isfile(phot_source_file):
            rsync_file(phot_source_file, phot_dest_file)
        else:
            print('WARNING: Could not find photometry file '+phot_source_file)

    print('Backed-up field photometric data products to '+params['output_dir'])

    upload_aws.upload_directory(data_set, local_root, aws_root)

    print('Backed-up field photometric data products to AWS Cloud')

def rsync_file(file_source_path, file_dest_path):

    cl = ['rsync', '-av', file_source_path, file_dest_path]

    p = subprocess.Popen(cl)
    stdoutput = p.communicate()[0]

    if stdoutput:
        print(stdoutput)

def get_args():
    """Function to acquire the necessary commandline arguments"""

    params = {}

    if len(sys.argv) == 1:

        params['dir_path'] = input('Please enter the path to the directory of a reduced dataset:')
        params['field_name'] = input('Please give the name of field: ')
        params['db_path'] = input('Please give the path to the photometry DB: ')
        params['output_dir'] = input('Please enter the directory path for output: ')
        params['local_root'] = input('Please enter the local root path (will be stripped off the path uploaded to AWS): ')
        params['aws_root'] = input('Please enter the AWS root path (will be prefixed to the path uploaded to AWS): ')

    else:

        params['dir_path'] = sys.argv[1]
        params['field_name'] = sys.argv[2]
        params['db_path'] = sys.argv[3]
        params['output_dir'] = sys.argv[4]
        params['local_root'] = sys.argv[5]
        params['aws_root'] = sys.argv[6]

    return params

if __name__ == '__main__':

    params = get_args()
    backup_field_photometry_products(params)
