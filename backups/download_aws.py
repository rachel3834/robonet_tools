# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 17:22:24 2018

@author: rstreet
"""

import os
import sys
from commands import getstatusoutput
import aws_cloud_config

def download_directory(aws_dir_path, local_dir_path):
    """Function to download a data directory from the AWS Cloud"""
    
    aws_config = aws_cloud_config.get_aws_config()
    
    check_sanity(local_dir_path)
    
    file_list = fetch_directory_file_list(aws_config, aws_dir_path)

    for entry in file_list:
        print(entry)
        
def fetch_directory_file_list(aws_config, aws_dir_path):
    """Function to query the AWS Cloud to return a contents listing of 
    the directory path given"""
    
    dir_listing = fetch_directory_listing(aws_config, aws_dir_path)
    
    cont = True
    
    while cont:
        
        cont = False
        
        sub_listing = []
        
        for entry in dir_listing:
            
            if '/' in entry[-1]:
                
                sub_aws_path = os.path.join(aws_dir_path,entry)
                
                sub_dir_list = fetch_directory_listing(aws_config, sub_aws_path)
                
                sub_listing += sub_dir_list
                
                for item in sub_dir_list:
                    print item
                    if '/' in item[-1]:

                        cont = True
        
        dir_listing += sub_listing
    
    return dir_listing
    
def fetch_directory_listing(aws_config, aws_dir_path):
    """Function to query the AWS Cloud to return a contents listing of 
    the directory path given"""
    
    if aws_dir_path[0:1] == '/':
        aws_dir_path = aws_dir_path[1:]
    
    aws_url = os.path.join(aws_config.bucket, aws_dir_path)
    
    cl = 'aws --profile='+aws_config.profile+' s3 ls '+aws_url
        
    (iexec,coutput) = getstatusoutput(cl)
    
    dir_listing = []
    
    for line in coutput.split('\n'):
        entry = line.replace('PRE','').lstrip()
        
        dir_listing.append(entry)

    return dir_listing
    
def check_sanity(local_dir_path):
    
    if os.path.isdir(local_dir_path) == False:
        raise IOError('Cannot find local directory '+local_dir_path)
        exit()
    
if __name__ == '__main__':
    
    if len(sys.argv) < 3:
        
        aws_dir_path = raw_input('Please enter the path to the AWS directory to download (without the AWS bucket prefix): ')
    
        local_dir_path = raw_input('Please enter the local directory to which data will be downloaded: ')
        
    else:
        
        aws_dir_path = sys.argv[1]
        
        local_dir_path = sys.argv[2]
        
    download_directory(aws_dir_path, local_dir_path)
    