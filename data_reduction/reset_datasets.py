from os import path, remove, rmdir
from sys import argv
import glob
from pyDANDIA import config_utils
from pyDANDIA import logs
from pyDANDIA import automatic_pipeline
from pyDANDIA import metadata
from pyDANDIA import reset_stage_metadata

def reset_pydandia_reductions():

    params = get_args()

    config = config_utils.build_config_from_json(params['config_file'])

    log = logs.start_pipeline_log(config['log_dir'], 'dataset_reset')

    running_processes = automatic_pipeline.read_process_list(config,log)

    datasets = read_dataset_list(params,log)

    if params['reset_code'] == 'stage2_no_ref':
        reset_stage2_no_ref(params, datasets, log)

    elif params['reset_code'] == 'full_reset':
        full_reset(params, datasets, log)

    else:
        print('Unrecognized reset code.  No action taken.')

    log.info('Completed reset process')
    logs.close_log(log)

def get_args():

    params = {}
    if len(argv) == 1:
        print("""Supported reset codes:
                stage2_no_ref: Reset reduction where stage 2 has failed to find a reference image
                full_reset: Remove all previous reduction products
              """)
        params['config_file'] = input('Please enter the path to the auto configuration file: ')
        params['datasets_file'] = input('Please enter the path to the file of datasets to reset: ')
        params['reset_code'] = input('Please enter the reset code: ')
    else:
        params['config_file'] = argv[1]
        params['datasets_file'] = argv[2]
        params['reset_code'] = argv[3]

    return params

def get_dataset_list(params,log):

    entries = glob.glob(path.join(params['red_dir'],'*'))

    exclude_list = ['lightcurves', 'config', 'phot_dbs', 'downloads']

    log.info('Identified the following reduction datasets:')

    datasets = []
    for d in entries:
        if path.isdir(d) and d not in exclude_list:
            datasets.append(d)
            log.info(d)

    return datasets

def read_dataset_list(params,log):

    if path.isfile(params['datasets_file']) == True:

        log.info('Found the list of datasets to be reset: '+params['datasets_file'])

        file_lines = open(params['datasets_file']).readlines()

        datasets = {}

        log.info('Going to reset the following datasets:')

        for line in file_lines:

            if len(line.replace('\n','')) > 0:
                (dataset_code, ref_status) = line.replace('\n','').split()
                datasets[dataset_code] = ref_status

            log.info(dataset_code)

    else:

        raise IOError('Cannot find input list of datasets.  Looking for '+params['datasets_file'])

    return datasets

def reset_stage2_no_ref(params, datasets, log):

    for data_dir in datasets:

        active_reduction = is_reduction_active(data_dir, running_processes)
        later_stage_logs = check_no_later_stage_logs(data_dir, 'reference_astrometry')
        ref_image = has_ref_image_been_selected(data_dir)

        log.info(path.basename(data_dir)+': ')
        log.info(' -> Active reduction? '+repr(active_reduction))
        log.info(' -> Later stage logs present? '+repr(later_stage_logs))
        log.info(' -> Reference image present? '+repr(ref_image))

        if not active_reduction and \
            not later_stage_logs and \
                not ref_image:

            remove(path.join(data_dir,'pyDANDIA_metadata.fits'))
            remove(path.join(data_dir,'stage0.log'))
            remove(path.join(data_dir,'stage1.log'))
            remove(path.join(data_dir,'stage2.log'))
            rmdir(path.join(data_dir, 'logs'))
            remove(path.join(data_dir, 'dataset.lock'))

            log.info(' ==> Reset reduction')

def full_reset(params, datasets, log):

    check = input('WARNING!  You are about to remove ALL data products.  Are you sure?  Y or n: ')

    if check == 'Y':
        for data_dir in datasets:

            active_reduction = is_reduction_active(data_dir, running_processes)

            log.info(path.basename(data_dir)+': ')
            log.info(' -> Active reduction? '+repr(active_reduction))

            if not active_reduction:

                for extn in ['*.dat', '*.png', '*.fits', '*.log', '*.lock', '*.hdf5']:
                    file_list = glob.glob(path.join(data_dir,extn))
                    for entry in file_list:
                        remove(entry)

                for sub_dir in ['logs', 'ref', 'resampled', 'kernel', 'diffim', 'lightcurves']:
                    rmdir(path.join(data_dir,sub_dir))

                log.info(' ==> Reset reduction')
    else:
        print('Reset aborted')


def check_no_later_stage_logs(data_dir, next_stage):
    """Function to check a reduction directory to see if stage logs exist for
    a stage later than the maximum stage indicated, implying that the reduction
    has proceeded beyond that point.

    param data_dir string  Full path to dataset reduction directory
    param next_stage string Name of the stage after the maximum allowed stage
    """

    return path.isfile(path.join(data_dir,next_stage+'.log'))

def is_reduction_active(data_dir, running_processes):

    if path.basename(data_dir) in running_processes.keys():
        return True
    else:
        return False

def has_ref_image_been_selected(data_dir):

    stage2_log = path.join(data_dir, 'stage2.log')
    found_log = path.isfile(stage2_log)

    if found_log:
        found_string = find_entry_in_log(log_file, 'No reference image found')
    else:
        found_string = False

    reduction_metadata = metadata.MetaData()
    reduction_metadata.load_all_metadata(setup.red_dir, 'pyDANDIA_metadata.fits')

    found_meta = 'REF_PATH' in reduction_metadata.data_architecture[1].keys()

    found_dir = path.isdir(data_dir, 'ref')
    if found_dir:
        image_list = glob.glob(data_dir, 'ref', '*fits')
        found_image = len(image_list) >= 1
    else:
        found_image = False

    if found_log and found_string and found_meta and found_dir and found_image:
        return True
    else:
        return False

def find_entry_in_log(log_file, search_string):

    log_data = open(log_file,'r').readlines()

    found_string = False
    for line in log_data:
        if search_string in line:
            found_string = True
            break

    return found_string

if __name__ == '__main__':
    reset_pydandia_reductions()
