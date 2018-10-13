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
    