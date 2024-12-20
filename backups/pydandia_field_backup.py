import os
import sys
import glob
import subprocess
import upload_aws
import compress_data_products
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

    red_dirs = glob.glob(os.path.join(params['dir_path'],params['field_name']+'_???-dom?-1m0-??-f???_??'))
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
            comp_phot_dest_file = compress_data_products.bzip2_file(phot_dest_file)
        else:
            print('WARNING: Could not find photometry file '+phot_source_file)

        # Backup the metadata
        meta_source_file = os.path.join(dir,'pyDANDIA_metadata.fits')
        meta_dest_file = os.path.join(staging_dir, 'pyDANDIA_metadata.fits')
        if os.path.isfile(meta_source_file):
            rsync_file(meta_source_file, meta_dest_file)
        else:
            print('WARNING: Could not find metadata file '+meta_source_file)

        # Backup the reference directory image data
        if not os.path.isdir(os.path.join(staging_dir,'ref')):
            os.mkdir(os.path.join(staging_dir,'ref'))

        if os.path.isfile(meta_source_file) and os.path.isdir(os.path.join(dir, 'ref')):
            mdata = metadata.MetaData()
            mdata.load_a_layer_from_file(dir, 'pyDANDIA_metadata.fits', 'data_architecture')
            source_files = []
            dest_files = []
            ref_source_file = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],mdata.data_architecture[1]['REF_IMAGE'][0])
            dref_source_file = ref_source_file.replace('.fits', '_res.fits')
            psfstamp_source = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'final_psf_master_stamp.fits')
            psfstampvar_source = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'final_psf_master_stamp_varience.fits')
            mask_source = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'master_mask.fits')
            maskedref_source = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'masked_ref_image.fits')
            psfmodel_source = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'psf_model.fits')
            psfmodelnorm_source = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'psf_model_normalized.fits')
            psfmodelres_source = os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'psf_model_residuals.fits')
            flist = [ref_source_file, dref_source_file,
                    psfstamp_source, psfstampvar_source,
                    mask_source, maskedref_source,
                    psfmodel_source, psfmodelnorm_source, psfmodelres_source]

            for f in flist:
                source_files.append(f)
                dest_files.append(os.path.join(staging_dir, 'ref', os.path.basename(f)))

            for i,f in enumerate(source_files):
                rsync_file(f, dest_files[i])

            refpngs = glob.glob(os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'ref','*.png'))
            for f_source in refpngs:
                f_dest = os.path.join(staging_dir, 'ref', os.path.basename(f_source))
                rsync_file(f_source, f_dest)

            # Backup the DS9 overlay datafiles
            overlays = glob.glob(os.path.join(mdata.data_architecture[1]['REF_PATH'][0],'*.reg'))
            for f_source in overlays:
                f_dest = os.path.join(staging_dir, 'ref', os.path.basename(f_source))
                rsync_file(f_source, f_dest)

        # Compress reference image
        for image in glob.glob(os.path.join(staging_dir, 'ref', '*.fits')):
            if os.path.isfile(image+'.bz2') == False:
                args = ['bzip2', image]
                p = subprocess.Popen(args, stdout=subprocess.PIPE)
                p.wait()
            else:
                print('Skipping compression of '+image+' (compressed product already exists)')

        # Backup the reduction logs
        log_list = glob.glob(os.path.join(dir,'*.log'))
        for log_source in log_list:
            log_dest = os.path.join(staging_dir,os.path.basename(log_source))
            rsync_file(log_source, log_dest)

        # Backup the output plots
        png_list = glob.glob(os.path.join(dir,'*.png'))
        for png_source in png_list:
            png_dest = os.path.join(staging_dir,os.path.basename(png_source))
            rsync_file(png_source, png_dest)

    # Backup the configuration directory
    config_source = os.path.join(params['dir_path'],'config')
    config_dest = params['output_dir']
    rsync_file(config_source, config_dest)

    # Backup the field crossmatch file and log:
    source_files = [ os.path.join(params['dir_path'], params['field_name']+'_field_crossmatch.fits'),
                     os.path.join(params['dir_path'], 'crossmatch.log'),
                     os.path.join(params['dir_path'], 'crossmatch_gaia.log'),
                     os.path.join(params['dir_path'], 'field_photometry.log')]
    logs_dest = params['output_dir']
    for f in source_files:
        rsync_file(f, logs_dest)
        if '.fits' in f:
            comp_file = compress_data_products.bzip2_file(os.path.join(logs_dest, os.path.basename(f)))

    # Backup the field photometry files:
    source_files = []
    for q in [1,2,3,4]:
        source_files.append( os.path.join(params['dir_path'], params['field_name']+'_quad'+str(q)+'_photometry.hdf5') )
    source_files.append(os.path.join(params['dir_path'], params['field_name']+'_star_dataset_normalizations.hdf5'))
    logs_dest = params['output_dir']
    for f in source_files:
        rsync_file(f, logs_dest)
        comp_file = compress_data_products.bzip2_file(os.path.join(logs_dest, os.path.basename(f)))

    # Backup the logs directory
    logs_source = os.path.join(params['dir_path'],'logs')
    logs_dest = params['output_dir']
    rsync_file(logs_source, logs_dest)

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
