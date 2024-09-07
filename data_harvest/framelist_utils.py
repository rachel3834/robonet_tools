from os import path
from sys import argv
import json
import log_utils
import config_utils
from astropy.io import fits
import glob

class Frame:

    def __init__(self, params = None, header = None):
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

        self.param_mapping = {'url': 'url',
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

        self.header_mapping = {'url': 'url',
                            'ORIGNAME': 'filename',
                            'DATE-OBS': 'dateobs',
                            'PROPID': 'proposalid',
                            'INSTRUME': 'instrument',
                            'OBJECT': 'object',
                            'SITEID': 'site',
                            'TELID': 'telescope',
                            'EXPTIME': 'exptime',
                            'FILTER': 'filter',
                            'REQNUM': 'reqnum'}
        if params != None:
            self.set_params(params)

        if header != None:
            self.set_header_params(header)

    def set_params(self, params):

        for key, attribute in self.param_mapping.items():
            if key in params.keys():
                if attribute == 'filename':
                    #value = str(params[key]).replace('.fz','')
                    value = str(params[key])
                else:
                    value = params[key]
                setattr(self,attribute,value)

    def set_header_params(self, header):

        for key, attribute in self.header_mapping.items():
            if key not in ['url']:
                try:
                    if attribute == 'filename':
                        value = str(header[key]).replace('.fz','')
                    else:
                        value = header[key]
                    setattr(self,attribute,value)
                except KeyError:
                    pass

    def summary(self):
        # Filename  date-obs   proposal  site  telescope  instrument filter exptime[s] object  reqnum
        return str(self.filename).replace('.fz','')+' '+self.dateobs+' '+self.proposalid+' '+\
                    self.site+' '+self.telescope+' '+self.instrument+' '+\
                    self.filter+' '+str(self.exptime)+' '+self.object+' '+str(self.reqnum)

def _search_key_in_filename(filename, search_keys):

    status = False
    for key in search_keys:
        if key in filename:
            status = True
            continue

    return status

def is_frame_calibration_data(filename):
    """Function to determine whether or not a given frame is a calibration or
    science frame, based on its filename"""

    search_keys = [ 'bias', 'dark', 'flat', 'skyflat' ]

    status = _search_key_in_filename(filename, search_keys)

    return status

def is_frame_spectrum(filename):
    """
    Function to determine whether the frame contains spectroscopy data, which is handled by
    a separate process.
    """

    # Check for instrument names that indicate FLOYDS data
    search_keys = ['_en']

    status = _search_key_in_filename(filename, search_keys)

    return status

def build_imaging_frame_list(config,query_results, proposal, new_frames, log):
    """Function to add new frames to a list of frames to be downloaded,
    excluding certain data types such as calibration frames.  This function is
    restricted to imaging data only"""

    fcount = 0
    for entry in query_results['results']:
        if entry['PROPID'] in config['proposal_ids'] and \
            not is_frame_calibration_data(entry['filename']) and \
            not is_frame_spectrum(entry['filename']):
            f = Frame(params=entry)
            new_frames.append(f)
            fcount += 1

    log.info('Found '+str(fcount)+' frame(s) available from proposal '+proposal)

    return new_frames

def build_floyds_frame_list(config,query_results, proposal, new_frames, log):
    """
    Function to add the names of new FLOYDS frames to a list of frames to be downloaded.
    This function is designed to handle the output from the FLOYDS data pipeline.
    """

    fcount = 0
    for entry in query_results['results']:
        if entry['PROPID'] in config['proposal_ids']:
            f = Frame(params=entry)
            new_frames.append(f)
            fcount += 1

    log.info('Found '+str(fcount)+' frame(s) available from proposal '+proposal)

    return new_frames

def build_goodman_frame_list(config,query_results, proposal, new_frames, log):
    """
    Function to add the names of new Goodman frames to a list of frames to be downloaded.
    This function is designed to handle the output from the Goodman spectroscopic data pipeline.
    """

    fcount = 0
    for entry in query_results['results']:
        if entry['proposal_id'] in config['soar_proposal_ids'] \
                and 'wecfzst' in entry['basename']:
            f = Frame(params=entry)
            new_frames.append(f)
            fcount += 1

    log.info('Found '+str(fcount)+' frame(s) available from proposal '+proposal)

    return new_frames

def output_frame_list(config, frame_list, log=None):

    f = open(config['frame_list'],'w')
    f.write('# Filename  date-obs   proposal  site  telescope  instrument filter exptime[s] object  reqnum\n')
    for fname,frame in frame_list.items():
        f.write(frame.summary()+'\n')
    f.close()

    if log!=None:
        log.info('Updated list of downloaded data')

def retrieve_image_header(file_path):
    """Images occasionally decompress with two image headers instead of one.
    This code selects the correct one"""

    hdu = fits.open(file_path)

    if len(hdu[0].header) > 100:
        header = hdu[0].header
    else:
        header = hdu[1].header

    return header

def refresh_frame_list():

    if len(argv) == 1:
        configfile = input('Please enter the path to the configuration file: ')
    else:
        configfile = argv[1]

    config = config_utils.get_config(configfile)

    frames_dict = {}

    file_list = glob.glob(path.join(config['data_download_dir'],'*.fits'))
    for file_path in file_list:
        header = retrieve_image_header(file_path)
        f = Frame(header=header)
        frames_dict[path.basename(file_path)] = f

    dir_list = glob.glob(path.join(config['data_reduction_dir'],'*'))
    for dir in dir_list:
        if path.isdir(dir) and path.isdir(path.join(dir,'data')):
            file_list = glob.glob(path.join(dir,'data','*.fits'))
            for file_path in file_list:
                header = retrieve_image_header(file_path)
                f = Frame(header=header)
                f.filename = path.basename(file_path).replace('.fz','')
                if path.basename(file_path) not in frames_dict.keys():
                    frames_dict[path.basename(file_path)] = f

    args = {'frame_list': 'updated_frame_list.txt'}
    output_frame_list(args, frames_dict, log=None)

    print('Updated frame list output to ./updated_frame_list.txt')

if __name__ == '__main__':
    refresh_frame_list()
