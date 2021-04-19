from os import path
from sys import argv
import config_utils
import boto3

def fetch_realtime_lightcurves_aws():

    config = get_args()

    config = get_credentials(config)

    s3_client = start_s3_client(config)

    file_list = list_available_lightcurves(config, s3_client)

    download_available_lightcurves(config, s3_client, file_list)

def get_args():
    if len(argv) != 2:
        config_file = input('Please enter the path to this scripts configuration file: ')
    else:
        config_file = argv[1]

    config = config_utils.get_config(config_file)

    return config

def get_credentials(config):

    home_dir = path.expanduser("~")

    credentials = open(path.join(home_dir, '.aws', 'credentials')).readlines()
    print(credentials)

    for i,line in enumerate(credentials):
        if config['awsid'] in line:
            config['aws_access_key_id'] = credentials[i+1].split('=')[-1].replace('\n','').lstrip()
            config['aws_secret_access_key'] = credentials[i+2].split('=')[-1].replace('\n','').lstrip()
            break

    if 'aws_access_key_id' not in config.keys():
        raise IOError('No AWS credentials found for user ID '+config['awsid'])

    print(config)

    return config

def start_s3_client(config):

    s3_client = boto3.client(
                's3',
                aws_access_key_id=config['aws_access_key_id'],
                aws_secret_access_key=config['aws_secret_access_key']
            )

    return s3_client

def list_available_lightcurves(config, s3_client):

    aws_path = path.join(config['aws_bucket'], 'OMEGA', 'realtime_lightcurves/')
    aws_path = path.join(config['aws_bucket'])
    log_path = path.join(config['local_data_dir'], 'available_lightcurves.txt')

    response = s3_client.list_objects(Bucket=config['aws_bucket'], Prefix='OMEGA/realtime_lightcurves')

    file_list = []
    for entry in response['Contents']:
        file_list.append(entry['Key'])

    return file_list

def download_available_lightcurves(config, s3_client, file_list):

    for file_object in file_list:
        file_name = path.join(config['local_data_dir'], path.basename(file_object))

        with open(file_name, 'wb') as f:
            s3_client.download_fileobj(config['aws_bucket'], file_object, f)


if __name__ == '__main__':
    fetch_realtime_lightcurves_aws()
