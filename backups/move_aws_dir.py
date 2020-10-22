import os
import sys
import socket
import subprocess
import aws_cloud_config
import aws_utils

def change_dir_path():

    (src_dir, dest_dir) = get_args()

    aws_config = aws_cloud_config.get_aws_config()

    contents = aws_utils.list_all_subdirs(aws_config,src_dir)

    print(contents)

    for src_path in contents['FILES']:
        dest_path = os.path.join(dest_dir, os.path.basename(src_path))
        dest_path = src_path.replace(src_dir, dest_dir)
        aws_utils.aws_move(aws_config, src_path, dest_path)

def get_args():

    if len(sys.argv) == 1:
        src_dir = input('Please enter the AWS path to the directory to move: ')
        dest_dir = input('Please enter the new AWS directory path: ')
    else:
        src_dir = sys.argv[1]
        dest_dir = sys.argv[2]

    return src_dir, dest_dir


if __name__ == '__main__':
    change_dir_path()
