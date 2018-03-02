# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:57:58 2018

@author: rstreet
"""

import os
import sys
import socket
from commands import getstatusoutput

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
                            print key, data
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
    print(config.profile)
    print(config.bucket)
    return config

def get_config_path():
    """Function to resolve the path to the correct configuration file, 
    depending on system environment"""
    
    host_name = socket.gethostname()
    
    if 'rachel' in str(host_name).lower():
        
        config_file_path = '/Users/rstreet/software/robonet_tools/configs/aws.xml'
    
    elif 'einstein' in str(host_name).lower():
        
        config_file_path = '/data/romerea/configs/aws.xml'
    
    return config_file_path

def aws_cp(aws_config, local_file_path, local_root, aws_root):
    """Function to compose an awscli commandline for a given file"""
    
    aws_path = os.path.join(aws_config.bucket, local_file_path.replace(local_root,aws_root))
    
    cl = 'aws --profile='+aws_config.profile+' s3 cp '+local_file_path+' '+aws_path

    (iexec,coutput) = getstatusoutput(cl)
    
    print(coutput)
    
def upload_directory(dir_path, local_root, aws_root):
    """Function to upload the contents of a data directory to AWS, 
    including all files and sub-directories.  
    """
    
    aws_config = get_aws_config()
    
    for root, dirs, files in os.walk(dir_path):
    
        for name in files:
            
            local_file_path = os.path.join(root,name)
            
            if not os.path.islink(local_file_path):
                
                aws_cp(aws_config,local_file_path,local_root,aws_root)


if __name__ == '__main__':
    
    if len(sys.argv) < 4:
        
        dir_path = raw_input('Please enter the path to the directory to upload: ')
    
        local_root = raw_input('Please enter the local root path (will be stripped off the path uploaded to AWS): ')
        
        aws_root = raw_input('Please enter the AWS root path (will be prefixed to the path uploaded to AWS): ')
        
    else:
        
        dir_path = sys.argv[1]
        
        local_root = sys.argv[2]
        
        aws_root = sys.argv[3]
        
    upload_directory(dir_path, local_root, aws_root)
    