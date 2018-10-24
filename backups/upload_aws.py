# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:57:58 2018

@author: rstreet
"""

import os
import sys
import socket
from commands import getstatusoutput
import aws_cloud_config

def aws_cp(aws_config, local_file_path, local_root, aws_root):
    """Function to compose an awscli commandline for a given file"""

    if aws_root[0:1] == '/':
        aws_root = aws_root[1:]

    aws_path = os.path.join(aws_config.bucket, local_file_path.replace(local_root,aws_root))
        
    cl = 'aws --profile='+aws_config.profile+' s3 cp '+local_file_path+' '+aws_path
    
    print (cl)
    
    (iexec,coutput) = getstatusoutput(cl)
    
    print(coutput)
    
def upload_directory(dir_path, local_root, aws_root):
    """Function to upload the contents of a data directory to AWS, 
    including all files and sub-directories.  
    """
    
    aws_config = aws_cloud_config.get_aws_config()
    
    for dir_paths, subdirs, files in os.walk(dir_path):
    
        for d in dir_paths:
            
            if 'imred' in d or 'rawlc' in d:
                
                (has_tar,tarball_list) = check_tarball_exists(d)
                
                if has_tar and len(tarball_list) > 0:
                    
                    for ball in tarball_list:
                        
                        local_file_path = os.path.join(d,ball)
                        
                        aws_cp(aws_config,local_file_path,local_root,aws_root)
                else:
                    
                    upload_all_files_in_subdir(aws_config,aws_root,local_root,d,files)
                    
            else:
                
                upload_all_files_in_subdir(aws_config,aws_root,local_root,d,files)

def upload_all_files_in_subdir(aws_config,aws_root,local_root,subdir,files):
    """Function to upload all files in a given subdirectory"""
    
    for name in files:
                    
        local_file_path = os.path.join(subdir,name)
        
        if not os.path.islink(local_file_path) and \
            is_original_image(local_file_path) == False:
            
            aws_cp(aws_config,local_file_path,local_root,aws_root)

def is_original_image(file_path):
    """Function to check whether a given file is a raw or BANZAI-processed
    original image, as opposed to one that has been further processed by 
    the DIA pipeline.  Since all LCO data are archived in the Cloud anyway, 
    the original image data should be excluded from the upload."""
    
    result = False
    
    if 'fits' in file_path:
        
        (dir_path,filename) = os.path.split(file_path)
        subdir = os.path.basename(dir_path)
    
        if '20' in subdir[0:2]:
            
            if '-e90' in filename or '-e91' in filename:
                
                if '.red.' not in filename and '.geo.' not in filename and \
                    '.dif.' not in filename and '.bpm.' not in filename:
                        
                    result = True
    return result

def check_tarball_exists(dir_path):
    """Function to check for the existance of a tarball of data from a directory"""
    
    if 'rawlc' in dir_path and \
        os.path.isfile(os.path.join(dir_path,'lightcurves.tar')):
        
        return True, [os.path.isfile(os.path.join(dir_path,'lightcurves.tar'))]
        
    if 'imred' in dir_path:
        flist = glob.glob(os.path.join(dir_path,'.tar'))
        
        return True, flist
    
    return False, []
    
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
    