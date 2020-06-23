from os import path, makedirs
import logging
from datetime import datetime

def start_day_log(config,log_name):
    """Function to set up a logger object"""

    if path.isdir(config['log_dir']) == False:
        makedirs(config['log_dir'])

    log_date = datetime.utcnow().strftime("%Y-%m-%d")
    log_file = path.join(config['log_dir'], log_name+'_'+log_date+'.log')

    log = logging.getLogger( log_name )

    if len(log.handlers) == 0:
        log.setLevel( logging.INFO )
        file_handler = logging.FileHandler( log_file )
        file_handler.setLevel( logging.INFO )
        formatter = logging.Formatter( fmt='%(asctime)s %(message)s', \
                                    datefmt='%Y-%m-%dT%H:%M:%S' )
        file_handler.setFormatter( formatter )
        log.addHandler( file_handler )

    log.info('------------------------------------------------------------\n')
    log.info('Started '+log_name+' process\n')

    return log

def close_log(log):
    log.info( 'Processing complete\n' )
    logging.shutdown()
