import numpy as np
import requests
import os
import time
import sys
import shutil

sys.path.append('/data/robouser/software/robonet_site/scripts/')
import log_utilities


def parse_conf_file(conf_file_adress, conf_file_name):


    configpath = os.path.join(conf_file_adress, conf_file_name)
    
    parameter_dictionnary = {}
    
    with open(configpath,'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                if not line.startswith('\n'):

                    col = line.split()
		    
		    try:

                    	parameter_dictionnary[col[0]] = col[1]

		    except:
	
			parameter_dictionnary[col[0]] = None

    return parameter_dictionnary



def find_images_previously_dowmloaded(already_dowmload_list_adress, already_dowmload_list_name):

	downloaded_images = np.loadtxt(already_dowmload_list_adress + already_dowmload_list_name, dtype=str)
	
	return downloaded_images


def construct_header(archive_adress, api_token,  username, password):


	response = requests.post(archive_adress + api_token,
		   data = {'username': username,'password': password}
		   ).json()

	token = response.get('token')
	header = {'Authorization': 'Token ' + token}

	return header

def find_data_to_process(parameter_dictionnary, header):


	request = parameter_dictionnary['archive']
	request += parameter_dictionnary['api_frames']
	request += '?'


	request += 'limit=50&'

	if parameter_dictionnary['date_start']:
		
		request += 'start=' + parameter_dictionnary['date_start'] + '&'

	else :
		#Last week*2 data
		utc_last_week = time.time()-7*24*3600*2
		day_last_week = time.gmtime(utc_last_week)
		starting_day = time.strftime("%Y-%m-%d", day_last_week)
		
		request += 'start=' + starting_day + '&'

	if parameter_dictionnary['date_end']:
		
		request += 'end=' + parameter_dictionnary['date_end'] + '&'

	if parameter_dictionnary['object']:

		request += 'OBJECT=' + parameter_dictionnary['object'] + '&'
	request += 'RLEVEL=' + parameter_dictionnary['rlevel'] +'&'

	request += 'PROPID=' + parameter_dictionnary['proposal'] +'&' 

	request += 'OBSTYPE=' + parameter_dictionnary['obstype'] 


	

	response = requests.get(request, headers=header).json()

	return response


def download_from_archive(filename):
	#import pdb; pdb.set_trace()
	root_directory = '/archive/engineering/'

	pieces = filename.split('-')
	site = pieces[0][:3]+'/'
	camera = pieces[1]+'/'
	day = pieces[2]
	
	directory = '/archive/engineering/'+site+camera+day+'/processed/'

	
	shutil.copy(directory+filename,'../')

def download_needed_data(conf_file_adress, conf_file_name, already_dowmload_list_adress, already_dowmload_list_name, 
			 output_directory, logger=None):

	try :

		parameter_dictionnary = parse_conf_file(conf_file_adress, conf_file_name)
		
		archive_adress = parameter_dictionnary['archive']
		api_token = parameter_dictionnary['api_token']
		username = parameter_dictionnary['username']
		password = parameter_dictionnary['password']
	
		logger.info('Correctly found and parse the conf file. ')

	except:

		logger.error('No conf file found, aboard!')
		sys.exit(1)

	try :
	
		header = construct_header(archive_adress, api_token,  username, password)
		logger.info('Correctly set the authentification parameters')

	except :
	
		logger.error('Something goes wrong on the authentification, aboard!')
		sys.exit(1)

	try:


		already_downloaded_images = find_images_previously_dowmloaded(already_dowmload_list_adress, already_dowmload_list_name)
		logger.info('Find previous downloaded images')
	except:

		logger.error('Can not find previous dowmloaded images, aboard!')
		sys.exit(1)

	try :
		
		list_of_data = find_data_to_process(parameter_dictionnary, header)
	
		frames = list_of_data['results']
		logger.info('New data to download have been found.')		

	except :

		logger.error('Something goes wrong on new data search, aboard!')
		sys.exit(1)

	if len(frames) !=0:

		while True:
			for frame in frames:
				print 'Start to work on frame :'+frame['filename']
				#if frame['filename'] not in already_downloaded_images :

					
				#	with open(output_directory+frame['filename'], 'wb') as f:

				#		f.write(requests.get(frame['url']).content)

		
				
				#	with open(already_dowmload_list_adress + already_dowmload_list_name, 'a') as f:

				#		f.write(frame['filename']+'\n')

				download_from_archive(frame['filename'])	

			if list_of_data.get('next'):

				list_of_data = requests.get(list_of_data['next'], headers=header).json()
				frames = list_of_data['results']

			else :

				break

def decompress_data(output_directory):

	os.system('funpack '+output_directory+'*.fz')
	os.system('rm '+output_directory+'*.fz')

if __name__ == '__main__':
	
	config = {'log_directory':'/data/romerea/data/logs/2017', 'log_root_name':'data_harvester'}
	logger = log_utilities.start_day_log( config, 'data_collection', console=False )	

	download_needed_data('/data/romerea/data/images/incoming/Scripts/', 'Harverst_Conf.txt', '/data/romerea/data/images/incoming/Scripts/', 'Already_Download_List.txt', '/data/romerea/data/images/incoming/',logger)
	logger.info('Download data success!')
	decompress_data('/data/romerea/data/images/incoming/')
	logger.info('Decompress data success!')

			
