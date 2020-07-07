from os import path, makedirs
from sys import argv
import requests
import json
from datetime import datetime, timedelta
import log_utils
import config_utils

CONFIG_FILE = '/data/omega/configs/data_download_config.json'
if path.isfile(CONFIG_FILE) == False:
    CONFIG_FILE = path.join(path.expanduser('~'), 'software', 'robonet_tools',
                        'configs', 'data_download_config.json')

class Frame:

    def __init__(self, params = None):
        self.filename = None
        self.url = None
        self.dateobs = None
        self.proposalid = None
        self.site = None
        self.telescope = None
        self.instrument = None
        self.filter = None
        self.exptime = None
        self.object = None
        self.reqnum = None

        if params != None:
            self.set_params(params)

    def set_params(self, params):

        param_mapping = {'url': 'url',
                            'filename': 'filename',
                            'DATE_OBS': 'dateobs',
                            'PROPID': 'proposalid',
                            'INSTRUME': 'instrument',
                            'OBJECT': 'object',
                            'SITEID': 'site',
                            'TELID': 'telescope',
                            'EXPTIME': 'exptime',
                            'FILTER': 'filter',
                            'REQNUM': 'reqnum'}

        for key, attribute in param_mapping.items():
            if key in params.keys():
                setattr(self,attribute,params[key])

    def summary(self):
        return self.filename+' '+self.object+' '+self.dateobs+' '+self.proposalid+' '+\
                    self.site+' '+self.telescope+' '+self.instrument+' '+\
                    self.filter+' '+str(self.exptime)+' '+str(self.reqnum)

def search_archive_for_data():

    config = config_utils.get_config(CONFIG_FILE)

    log = log_utils.start_day_log(config,'data_download')

    downloaded_frames = read_frame_list(config, log)

    (start_time, end_time) = set_date_range(config, log)

    new_frames = fetch_new_frames(config, start_time, end_time, log)

    downloaded_frames = download_new_frames(config,new_frames,downloaded_frames,log)

    output_frame_list(config, downloaded_frames, log)

    log_utils.close_log(log)

def set_date_range(config, log):
    """Function to set the search date range from the configuration or to
    the last 24hrs"""

    if 'none' not in str(config['start_datetime']).lower() and \
        'none' not in str(config['end_datetime']).lower():
        start_time = datetime.strptime(config['start_datetime'],'%Y-%m-%d %H:%M')
        end_time = datetime.strptime(config['end_datetime'],'%Y-%m-%d %H:%M')
    else:
        end_time = datetime.utcnow()
        start_time = end_time - timedelta(days=1)

    log.info('Searching for data taken between '+start_time.strftime("%Y-%m-%d %H:%M")+\
                ' and '+end_time.strftime("%Y-%m-%d %H:%M"))

    return start_time, end_time

def read_frame_list(config, log):
    """Function to read the list of previously-downloaded frames"""

    downloaded_frames = {}
    if path.isfile(config['frame_list']):
        f = open(config['frame_list'],'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if '#' not in line[0:1]:
                f = Frame()
                entry = line.replace('\n','').split()
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

        log.info('Read list of '+str(len(downloaded_frames))+' frame(s)')
    else:
        log.info('No list of existing frames found')

    return downloaded_frames

def fetch_new_frames(config, start_time, end_time, log):
    """Function to query the archive and retrieve a list of frames,
    described as dictionaries of frame parameters."""

    new_frames = []
    for proposal in config['proposal_ids']:

        ur = { 'PROPID': proposal, 'start': start_time.strftime("%Y-%m-%d %H:%M"),
                                    'end': end_time.strftime("%Y-%m-%d %H:%M"),
                                    'RLEVEL': 91 }
        print(ur)
        results = talk_to_lco_archive(config, ur, 'frames', 'GET')
        print(results)

        fcount = 0
        for entry in results['results']:
            f = Frame(params=entry)
            new_frames.append(f)
            fcount += 1

        log.info('Found '+str(fcount)+' frame(s) available from proposal '+proposal)

    return new_frames

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

        if frame.filename not in downloaded_frames.keys():
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

def output_frame_list(config, downloaded_frames, log):

    f = open(config['frame_list'],'w')
    f.write('# Filename  date-obs   proposal  site  telescope  instrument filter exptime[s] object  reqnum\n')
    for fname,frame in downloaded_frames.items():
        f.write(frame.summary()+'\n')
    f.close()

    log.info('Updated list of downloaded data')

if __name__ == '__main__':
    search_archive_for_data()
