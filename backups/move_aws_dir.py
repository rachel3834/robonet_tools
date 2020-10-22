import os
import sys
import socket
import subprocess
import aws_cloud_config

def archive_old_results_aws():

    (src_dir, dest_dir) = get_args()

    
def get_args():

    if len(sys.argv) == 1:
        src_dir = input('Please enter the AWS path to the directory to move: ')
        dest_dir = input('Please enter the new AWS directory path: ')
    else:
        src_dir = sys.argv[1]
        dest_dir = sys.argv[2]

    return src_dir, dest_dir
