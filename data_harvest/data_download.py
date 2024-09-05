from os import path, makedirs
from sys import argv
import requests
import json
from datetime import datetime, timedelta
import log_utils
import config_utils
import framelist_utils

def search_archive_for_data(CONFIG_FILE):

    config = config_utils.get_config(CONFIG_FILE)

    log = log_utils.start_day_log(config,'data_download')

    downloaded_frames = read_frame_list(config, log)

    (start_time, end_time) = set_date_range(config, log)

    # Fetch imaging data
    new_frames = fetch_new_imaging_data(config, start_time, end_time, log)

    downloaded_frames = download_new_frames(config,new_frames,downloaded_frames,log)

    # Fetch spectroscopy data
    new_floyds_frames = fetch_new_floyds_data(config, start_time, end_time, log)

    downloaded_frames = download_new_frames(config, new_floyds_frames, downloaded_frames, log)

    framelist_utils.output_frame_list(config, downloaded_frames, log)

    log_utils.close_log(log)

def set_date_range(config, log):
    """Function to set the search date range from the configuration or to
    the last 24hrs

    Expected date-time format is %Y-%m-%d %H:%M"""

    if 'none' not in str(config['start_datetime']).lower() and \
        'none' not in str(config['end_datetime']).lower():
        start_time = datetime.strptime(config['start_datetime'],'%Y-%m-%d %H:%M')
        end_time = datetime.strptime(config['end_datetime'],'%Y-%m-%d %H:%M')
    else:
        end_time = datetime.utcnow()
        start_time = end_time - timedelta(days=1)

    log.info('Searching for data taken between '+start_time.strftime("%Y-%m-%d %H:%M")+\
                ' and '+end_time.strftime("%Y-%m-%d %H:%M"))
    dt = timedelta(hours=24.0)
    start_time -= dt

    return start_time, end_time

def read_frame_list(config, log):
    """Function to read the list of previously-downloaded frames
    Columns are:
    Filename  date-obs   proposal  site  telescope  instrument filter exptime[s] object  reqnum
    """

    downloaded_frames = {}
    if path.isfile(config['frame_list']):
        f = open(config['frame_list'],'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if '#' not in line[0:1]:
                f = framelist_utils.Frame()
                entry = line.replace('\n','').split()
                if len(entry) == 10:
                    f.filename = entry[0]
                    f.dateobs = entry[1]
                    f.proposalid = entry[2]
                    f.site = entry[3]
                    f.telescope = entry[4]
                    f.instrument = entry[5]
                    f.filter = entry[6]
                    f.exptime = entry[7]
                    f.object = entry[8]
                    f.reqnum = entry[9]

                    downloaded_frames[f.filename] = f

                elif entry[2] == 'calibrate':
                    log.info('Skipped calibrate frame '+entry[0])
                else:
                    log.info('ERROR processing frame log entry: '+line)

        log.info('Read list of '+str(len(downloaded_frames))+' frame(s)')
    else:
        log.info('No list of existing frames found')

    return downloaded_frames

def fetch_new_imaging_data(config, start_time, end_time, log):
    """Function to query the archive and retrieve a list of frames,
    described as dictionaries of frame parameters."""

    new_frames = []
    for proposal in config['proposal_ids']:

        ur = { 'PROPID': proposal, 'start': start_time.strftime("%Y-%m-%d %H:%M"),
                                    'end': end_time.strftime("%Y-%m-%d %H:%M"),
                                    'RLEVEL': 91 }

        results = talk_to_lco_archive(config, ur, 'frames', 'GET')

        # Build a list of new frames, excluding calibration data and spectra
        new_frames = framelist_utils.build_imaging_frame_list(config, results, proposal, new_frames, log)

    return new_frames

def fetch_new_floyds_data(config, start_time, end_time, log):
    """
    Function to query the archive specifically for FLOYDS data.
    This is handled separately because the format of the processed data products is quite
    different from imaging data.
    """
    floyds_instruments = ['en06', 'en12']
    new_floyds_frames = []

    for proposal in config['proposal_ids']:
        for inst_id in floyds_instruments:
            ur = { 'PROPID': proposal, 'start': start_time.strftime("%Y-%m-%d %H:%M"),
                                        'end': end_time.strftime("%Y-%m-%d %H:%M"),
                                        'instrument_id': inst_id,
                                        'RLEVEL': 90,
                                        'configuration_type': 'SPECTRUM'}

            results = talk_to_lco_archive(config, ur, 'frames', 'GET')

            # Build a list of new frames, excluding calibration data and spectra
            new_floyds_frames = framelist_utils.build_floyds_frame_list(
                config,
                results,
                proposal,
                new_floyds_frames,
                log)

    return new_floyds_frames

def talk_to_lco_archive(config,ur,end_point,method):
    """Function to communicate with various APIs of the LCO network.
    ur should be a user request while end_point is the URL string which
    should be concatenated to the observe portal path to complete the URL.
    Accepted end_points are:
        "userrequests"
        "userrequests/cadence"
        "userrequests/<id>/cancel"
    Accepted methods are:
        POST GET
    """

    jur = json.dumps(ur)

    headers = {'Authorization': 'Token ' + config['lco_archive_token']}

    if end_point[0:1] == '/':
        end_point = end_point[1:]
    if end_point[-1:] != '/':
        end_point = end_point+'/'
    url = path.join(config['lco_archive_url'],end_point)

    if method == 'POST':
        if ur != None:
            response = requests.post(url, headers=headers, json=ur).json()
        else:
            response = requests.post(url, headers=headers).json()
    elif method == 'GET':
        response = requests.get(url, headers=headers, params=ur).json()

    return response

def download_new_frames(config,new_frames,downloaded_frames,log):
    """Function to download list of frames"""

    headers = {'Authorization': 'Token ' + config['lco_archive_token']}

    for frame in new_frames:

        if frame.filename not in downloaded_frames.keys() and \
            not framelist_utils.is_frame_calibration_data(frame.filename) and \
            frame.proposalid in config['proposal_ids']:
            try:
                response = requests.get(frame.url)

                dframe = open(path.join(config['data_download_dir'],frame.filename),'wb')
                dframe.write(response.content)
                dframe.close()

                if response.status_code == 200:
                    downloaded_frames[frame.filename] = frame
                    log.info('-> Downloaded '+frame.filename)
                else:
                    log.info('->>>> Error downloading '+frame.filename+\
                                ' status code='+str(response.status_code))
            except requests.exceptions.ConnectionError:
                log.info('->>>> ConnectionError downloading '+frame.filename)
                
        else:
            log.info(frame.filename+' already downloaded')

    return downloaded_frames

def start_day_log(config):
    """Function to set up a logger object"""

    if path.isdir(config['log_dir']) == False:
        makedirs(config['log_dir'])

    log_date = datetime.utcnow().strftime("%Y-%m-%d")
    log_file = path.join(config['log_dir'], 'data_download_'+log_date+'.log')

    log = logging.getLogger( 'data_download' )

    if len(log.handlers) == 0:
        log.setLevel( logging.INFO )
        file_handler = logging.FileHandler( log_file )
        file_handler.setLevel( logging.INFO )
        formatter = logging.Formatter( fmt='%(asctime)s %(message)s', \
                                    datefmt='%Y-%m-%dT%H:%M:%S' )
        file_handler.setFormatter( formatter )
        log.addHandler( file_handler )

    log.info('------------------------------------------------------------\n')
    log.info('Started data download process\n')

    return log

def close_log(log):
    log.info( 'Processing complete\n' )
    logging.shutdown()

if __name__ == '__main__':

    if len(argv) == 1:
        CONFIG_FILE = input('Please enter the path to the configuration file: ')
    else:
        CONFIG_FILE = argv[1]

    search_archive_for_data(CONFIG_FILE)
