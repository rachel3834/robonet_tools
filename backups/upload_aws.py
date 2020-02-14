# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:57:58 2018

@author: rstreet
"""

import os
import sys
import socket
import subprocess
import aws_cloud_config

def aws_cp(aws_config, local_file_path, local_root, aws_root):
    """Function to compose an awscli commandline for a given file"""

    if aws_root[0:1] == '/':
        aws_root = aws_root[1:]

    aws_path = os.path.join(aws_config.bucket, local_file_path.replace(local_root,aws_root))

    cl = ['aws', '--profile='+aws_config.profile,'s3','cp',local_file_path,aws_path]

    print(local_file_path+' -->  '+aws_root)
    # WARNING: Cannot use stdout-> pipe as the cl output of awscli for larger
    # files can deadline the child process.
    p = subprocess.Popen(cl)
    #p.wait()
    stdoutput = p.communicate()[0]
    if stdoutput:
        print(stdoutput)

def upload_directory(dir_path, local_root, aws_root):
    """Function to upload the contents of a data directory to AWS,
    including all files and sub-directories.
    """

    aws_config = aws_cloud_config.get_aws_config()

    tarred_data = check_for_tarballs(dir_path)

    for root, dirs, files in os.walk(dir_path):

        for name in files:

            local_file_path = os.path.join(root,name)

            if 'fits' in local_file_path and '20160803' in local_file_path:
                print(os.path.islink(local_file_path))
                print(is_original_image(local_file_path))
                print(file_in_tarball(local_file_path, tarred_data))
                print(local_file_path)

            if os.path.islink(local_file_path) == False and \
                is_original_image(local_file_path) == False and \
                file_in_tarball(local_file_path, tarred_data) == False:

                aws_cp(aws_config,local_file_path,local_root,aws_root)

                if 'fits' in local_file_path:
                    print('-> Uploaded')

def check_for_tarballs(dir_path):
    """Function to check for the existance of a tarball of data from a directory"""

    tarred_data = {}

    tarred_data['/lc/'] = os.path.isfile(os.path.join(dir_path,'lightcurves.tar'))
    tarred_data['/imred/'] = os.path.isfile(os.path.join(dir_path,'imred_data.tar'))
    tarred_data['/gimred/'] = os.path.isfile(os.path.join(dir_path,'gimred_data.tar'))
    tarred_data['/dimred/'] = os.path.isfile(os.path.join(dir_path,'dimred_data.tar'))

    return tarred_data

def file_in_tarball(file_path, tarred_data):
    """Function to check whether a given file is contained in an available
    tarball or needs to be uploaded separately"""

    for d in tarred_data.keys():

        if d in file_path:

            return True

    return False

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

if __name__ == '__main__':

    if len(sys.argv) < 4:

        dir_path = input('Please enter the path to the directory to upload: ')

        local_root = input('Please enter the local root path (will be stripped off the path uploaded to AWS): ')

        aws_root = input('Please enter the AWS root path (will be prefixed to the path uploaded to AWS): ')

    else:

        dir_path = sys.argv[1]

        local_root = sys.argv[2]

        aws_root = sys.argv[3]

    upload_directory(dir_path, local_root, aws_root)
