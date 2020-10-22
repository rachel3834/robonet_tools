# Import required libraries
from glob import glob
from os import getcwd, path, remove
from sys import argv, exit
from sys import path as systempath
cwd = getcwd()
systempath.append(path.join(cwd, '../'))
from pyDANDIA import pipeline_setup
from pyDANDIA import automatic_pipeline
from pyDANDIA import logs
from pyDANDIA import config_utils
import shutil
import psutil

# Command line calling sequence: pass it the full path to an auto_pipeline_config.json file
# e.g. python clear_locks.py /full/path/to/config/auto_pipeline_config.json

def check_paths(reduction_path):
    ''' Check if the directories exist
    '''
    if path.isdir(reduction_path) == False:
        print ('Directory does not exist', str(reduction_path))
        return 1
    #if path.isdir(path.join(reduction_path, 'ref')) == False:
    #    print ('No ref/ subdirectory found for', str(reduction_path))
    #    return 2
    return 0

def read_last_line(path_to_file):
    ''' Read last line in a text file
    '''
    with open(path_to_file, 'r') as fh:
        last_line = fh.readlines()[-2:-1][0]
    return last_line


def clear_locks(config_file='/home/ytsapras/test_data/data_reduction/config/auto_pipeline_config.json', verbose=False):
    '''This function will perform a recursive sub-directory search and safely remove
       any dataset.lock files if the following conditions are met:
       (Note: The program will avoid all currently active reductions.)
       --> That the highest-number stage executed is stage2, and there are no logs for reference_astrometry
           or stage 3.  
       --> That there is no .fits image in the ref/ sub-directory, which may or may not exist
       --> Check the last entry of stage2.log. It should indicate that stage 2 is complete 
           ("Processing complete").
       --> Check that a dataset.lock file is present. 
       - If these conditions hold, the code should remove the dataset.lock file from the reduction 
           directory. 

    '''
    # Get configuration file
    config = config_utils.build_config_from_json(config_file)
    path_to_top_reduction_directory = config['data_red_dir'] 
    
    # Collect all subdir names, excluding logs/ and config/
    all_sub_directories = glob(path.join(path_to_top_reduction_directory,'[!logs!config]*/'))
    
    # Identify currently active reductions in the list of subdirectories
    log = logs.start_pipeline_log(config['log_dir'], 'clear_lock')
    active_reductions = automatic_pipeline.read_process_list(config, log)
    
    # For each remaining subdirectory, check conditions and
    # remove dataset.lock file if conditions are met
    for subdir in all_sub_directories:
        # Define conditions
        condition1 = False # Is the highest-number stage executed, stage 2?
        condition2 = False # Is there a fits image in ref/ ?
        condition3 = False # Did stage 2 complete successfully?
        condition4 = False # Is a dataset.lock file present? 
        
        # Loop over all subdirs not in active_reductions
        if subdir.split('/')[-2] not in active_reductions.keys():
            
            # If the subdir does not exist or does not have a ref/ subdirectory, skip it
            # If the path exists then this sets ans = 0
            ans = check_paths(subdir)
            
            if (ans == 1 or ans == 2):
                break
            
            if (path.exists(path.join(subdir, 'stage2.log'))):
                condition1 = True
            
            if (path.exists(path.join(subdir, 'stage3.log')) or 
                path.exists(path.join(subdir, 'reference_astrometry.log'))):
                condition1 = False
            
            # Check for fits images in ref/
            if (path.exists(path.join(subdir, 'ref'))):
                ref_img = glob(path.join(subdir, 'ref', '*.fits'))
                if len(ref_img) != 0:
                    condition2 = True
            
            # Check whether stage 2 completed successfully
            last_line_stage2 = read_last_line(path.join(subdir, 'stage2.log'))
            if 'Processing complete' in last_line_stage2:
                condition3 = True
            
            # Check if dataset.lock is present
            if (path.exists(path.join(subdir, 'dataset.lock'))):
                condition4 = True
            
            # Check conditions and apply prescribed remedy
            if (condition1 == True and  condition2 == False and \
                condition3 == True and condition4 == True):
                remove(path.join(subdir, 'dataset.lock'))
                log.info('Removed lock file from %s' % str(subdir))
                
                if verbose == True:
                    print ('Removed lock file from %s' % str(subdir))
                    
            else:
                log.info('Did not remove lock file from %s' % str(subdir))
                
                if verbose == True:
                    print ('Did not remove lock file from%s' % str(subdir))
        
    logs.close_log(log)
    return 1

def reset_data_dirs(config_file='/home/ytsapras/test_data/data_reduction/config/auto_pipeline_config.json', verbose=False):
    '''This function will perform a recursive sub-directory search and reset it for a new reduction from
       scratch. In effect, it will remove everything from the data directory except the directory itself 
       and the data/ folder which contains the images. !!!USE WITH CAUTION!!!

    '''
    # Get configuration file
    config = config_utils.build_config_from_json(config_file)
    path_to_top_reduction_directory = config['data_red_dir'] 
    
    # Collect all subdir names, excluding logs/ and config/
    all_sub_directories = glob(path.join(path_to_top_reduction_directory,'[!logs!config]*/'))
    
    # Keep a log 
    log = logs.start_pipeline_log(config['log_dir'], 'reset_reduction')
    
    # Loop over all subdirs removing everything except the data/ folders
    for subdir in all_sub_directories:
        contents = glob(path.join(subdir, '*'))
        contents.remove(path.join(subdir, 'data'))
        for item in contents:
            if path.isdir(item):
                try:
                    shutil.rmtree(item)
                    log.info("Removed folder: %s" % item)
                    if verbose == True:
                        print ("Removed folder: %s" % item)
                except OSError:
                    log.info("Unable to remove folder: %s" % item)
                    if verbose == True:
                        print ("Unable to remove folder: %s" % item)
                        print ("Check read/write folder permissions.")
            else:
                if path.exists(item):
                    try:
                        remove(item)
                        log.info("Removed file: %s" % item)
                        if verbose == True:
                            print ("Removed file: %s" % item)
                    except OSError:
                        log.info("Unable to remove file: %s" % item)
                        if verbose == True:
                            print("Unable to remove file: %s" % item)
                            print ("Check read/write file permissions.")
    
    logs.close_log(log)
    return 1
    

if __name__ == '__main__':
    config_file='/home/ytsapras/test_data/data_reduction/config/auto_pipeline_config.json'
    if path.exists(config_file):
        clear_locks(config_file=config_file, verbose=False)
    else:
        print('File', config_file, 'was not found.')
