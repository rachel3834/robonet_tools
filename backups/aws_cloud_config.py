# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 17:18:37 2018

@author: rstreet
"""
import os
import sys
import socket

class AWSConfig():
    """Configuration class for AWS Cloud"""

    def __init__(self):
        self.profile = None
        self.bucket = None

    def load_config(self,config_file_path):

        if os.path.isfile(config_file_path):

            file_lines = open(config_file_path,'r').readlines()

            for line in file_lines:

                data = line.replace('<',':::').replace('>',':::').split(':::')

                if len(data) > 3:

                    data = data[2]

                    for key in ['profile','bucket']:

                        if key in line:
                            print(key, data)
                            setattr(self,key,data)

        else:

            print('ERROR: Cannot find configuration file '+config_file_path)

            sys.exit()

def get_aws_config():
    """Function to obtain the necessary credentials to upload to
    a pre-existing AWS Bucket
    """

    config = AWSConfig()

    config_file_path = get_config_path()

    config.load_config(config_file_path)

    #print(config.profile)
    #print(config.bucket)

    return config

def get_config_path():
    """Function to resolve the path to the correct configuration file,
    depending on system environment"""

    host_name = socket.gethostname()

    if 'rachel' in str(host_name).lower():

        config_file_path = '/Users/rstreet1/ROMEREA/configs/aws.xml'

    elif 'einstein' in str(host_name).lower():

        config_file_path = '/data/romerea/configs/aws.xml'

    elif 'einstore' in str(host_name).lower():

        config_file_path = '/data00/robouser/romerea/configs/aws.xml'

    return config_file_path
