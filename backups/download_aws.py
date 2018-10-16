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
    
    dir_listing = fetch_directory_file_list(aws_config, aws_dir_path)

    download_file_list(aws_config, aws_dir_path, local_dir_path, dir_listing)
        
def fetch_directory_file_list(aws_config, aws_dir_path):
    """Function to query the AWS Cloud to return a contents listing of 
    the directory path given"""
    
    print('Querying Cloud directory list the contents to be downloaded (please be patient)...')
    
    dir_listing = fetch_directory_listing(aws_config, aws_dir_path)
    
    sub_dir_list = check_for_sub_directories(dir_listing)
    
    while len(sub_dir_list) > 0:
        
        for entry in sub_dir_list:

            if len(entry) > 1 and '/' in entry[-1]:
                
                sub_aws_path = os.path.join(aws_dir_path,entry)
                
                sub_dir_contents = fetch_directory_listing(aws_config, sub_aws_path)
                
                paths_list = [os.path.join(entry,sub_dir_contents[i]) for i in range(0,len(sub_dir_contents),1)]
                
                dir_listing += paths_list
                
                i = dir_listing.index(entry)
                tmp = dir_listing.pop(i)
                
        sub_dir_list = check_for_sub_directories(dir_listing)
    
    #dir_listing = [os.path.join(aws_dir_path,dir_listing[i]) for i in range(0,len(dir_listing),1)]
    
    print(' -> Found '+str(len(dir_listing))+' files to be downloaded')
    
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
        
        if len(entry) > 0:
            
            if '/' not in entry:
                
                entry_list = entry.split(' ')
                
                dir_listing.append(entry_list[-1])
                
            else:

                dir_listing.append(entry)
            
    return dir_listing

def check_for_sub_directories(dir_listing):
    
    sub_dir_list = []
    
    for item in dir_listing:
        
        if len(item) > 1 and '/' in item[-1]:
            
            sub_dir_list.append(item)
    
    return sub_dir_list
    
def check_sanity(local_dir_path):
    
    if os.path.isdir(local_dir_path) == False:
        raise IOError('Cannot find local directory '+local_dir_path)
        exit()

def download_file_list(aws_config, aws_dir_path, local_dir_path, dir_listing):
    """Function to retrieve a list of files from the AWS Cloud"""
    
    print('Downloading data from the Cloud...')
    
    if '/' in aws_dir_path[-1]:
        base_dir = os.path.basename(aws_dir_path[:-1])
    else:
        base_dir = os.path.basename(aws_dir_path)
    
    if aws_dir_path[0:1] == '/':
        aws_dir_path = aws_dir_path[1:]
    
    for f in dir_listing:
        
        aws_path = os.path.join(aws_config.bucket, aws_dir_path, f)
        local_path = os.path.join(local_dir_path, base_dir, f)
        
        cl = 'aws --profile='+aws_config.profile+' s3 cp '+aws_path+' '+local_path
    
        (iexec,coutput) = getstatusoutput(cl)
    
        print(coutput)

if __name__ == '__main__':
    
    if len(sys.argv) < 3:
        
        aws_dir_path = raw_input('Please enter the path to the AWS directory to download (without the AWS bucket prefix): ')
    
        local_dir_path = raw_input('Please enter the local directory to which data will be downloaded: ')
        
    else:
        
        aws_dir_path = sys.argv[1]
        
        local_dir_path = sys.argv[2]
        
    download_directory(aws_dir_path, local_dir_path)
    