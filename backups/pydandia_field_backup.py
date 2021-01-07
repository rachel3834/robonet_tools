import os
import sys
import glob
import subprocess
import upload_aws
from pyDANDIA import metadata

def backup_field_photometry_products(params):
    """Function to gather all the photometric data products for all datasets
    for a given field and back them up to another disk"""

    if params['output_dir'] != None and not os.path.isdir(params['output_dir']):
        raise IOError('Cannot find output directory '+params['output_dir'])

    if params['db_path'] and os.path.isfile(params['db_path']):
        rsync_file(params['db_path'], os.path.join(params['output_dir'], os.path.basename(params['db_path'])))
    elif not params['db_path']:
        print('No photometry DB to be backed up')
    else:
        print('WARNING: Could not find field database at '+params['db_path'])

    red_dirs = glob.glob(os.path.join(params['dir_path'],params['field_name']+'*'))
    print('Found '+str(len(red_dirs))+' reduction directories')

    for dir in red_dirs:
        dset = os.path.basename(dir)
        staging_dir = os.path.join(params['output_dir'],dset)
        if not os.path.isdir(staging_dir):
            os.mkdir(staging_dir)

        # Backup the main photometry file
        phot_source_file = os.path.join(dir,'photometry.hdf5')
        phot_dest_file = os.path.join(staging_dir, 'photometry.hdf5')
        if os.path.isfile(phot_source_file):
            rsync_file(phot_source_file, phot_dest_file)
        else:
            print('WARNING: Could not find photometry file '+phot_source_file)

        # Backup the metadata
        meta_source_file = os.path.join(dir,'pyDANDIA_metadata.fits')
        meta_dest_file = os.path.join(staging_dir, 'pyDANDIA_metadata.fits')
        if os.path.isfile(meta_source_file):
            rsync_file(meta_source_file, meta_dest_file)
        else:
            print('WARNING: Could not find metadata file '+meta_source_file)

        # Backup the reference image
        mdata = metadata.MetaData()
        mdata.load_a_layer_from_file(dir, 'pyDANDIA_metadata.fits', 'data_architecture')
        ref_source_file = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],mdata.data_architecture[1]['REF_IMAGE'][0])
        ref_dest_file = os.path.join(staging_dir, 'ref', mdata.data_architecture[1]['REF_IMAGE'][0])

        if not os.path.isdir(os.path.join(staging_dir,'ref')):
            os.mkdir(os.path.join(staging_dir,'ref')
        if os.path.isfile(ref_source_file):
            rsync_file(ref_source_file, ref_dest_file)
        else:
            print('WARNING: Could not find reference image '+ref_source_file)

        # Compress reference image
        if os.path.isfile(ref_dest_file+'.bz2') == False:
            args = ['bzip2', ref_dest_file]
            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            p.wait()
        else:
            print('Skipping compression of '+ref_dest_file+' (compressed product already exists)')

        print('Backed-up data products for '+dset+' to '+staging_dir)

    upload_aws.upload_directory_nochecks(params['output_dir'],
                                params['local_root'],
                                params['aws_root'])

    print('Backed-up field photometric data products to AWS Cloud')

def rsync_file(file_source_path, file_dest_path):

    print('Rsyncing '+file_source_path+' to '+file_dest_path)

    cl = ['rsync', '-av', file_source_path, file_dest_path]

    p = subprocess.Popen(cl)
    stdoutput = p.communicate()[0]

    if stdoutput:
        print(stdoutput)

def get_args():
    """Function to acquire the necessary commandline arguments"""

    def verify_path_slashes(dir_path):

        if not dir_path[0:1] == '/':
            dir_path = '/' + dir_path

        if not dir_path[-1:] == '/':
            dir_path = dir_path + '/'

        return dir_path

    params = {'db_path': None}

    if len(sys.argv) == 1:

        params['dir_path'] = input('Please enter the path to the directory of a reduced dataset:')
        params['field_name'] = input('Please give the name of field: ')
        params['db_path'] = input('Please give the path to the photometry DB or None: ')
        params['output_dir'] = input('Please enter the directory path for output: ')
        params['local_root'] = input('Please enter the local root path (will be stripped off the path uploaded to AWS): ')
        params['aws_root'] = input('Please enter the AWS root path (will be prefixed to the path uploaded to AWS): ')

        if 'none' in str(params['db_path']).lower():
            params['db_path'] = None

    else:

        params['dir_path'] = sys.argv[1]
        params['field_name'] = sys.argv[2]
        params['db_path'] = sys.argv[3]
        params['output_dir'] = sys.argv[4]
        params['local_root'] = sys.argv[5]
        params['aws_root'] = sys.argv[6]

        if 'none' in str(params['db_path']).lower():
            params['db_path'] = None

        for a in sys.argv:
            if 'db_path=' in a:
                params['db_path'] = str(a).split('=')[-1]

    params['local_root'] = verify_path_slashes(params['local_root'])
    params['aws_root'] = verify_path_slashes(params['aws_root'])

    return params

if __name__ == '__main__':

    params = get_args()
    backup_field_photometry_products(params)
