import os
import sys
import aws_cloud_config
import aws_utils
import argparse

def remove_aws_dir():

    parser = argparse.ArgumentParser()
    parser.add_argument('dir_path', help='AWS path to the directory to remove')
    args = parser.parse_args()

    aws_config = aws_cloud_config.get_aws_config()

    print(args.dir_path)
    opt = input('WARNING: This script permanently deletes data from AWS!  Continue Y or n? ')
    if opt == 'Y':
        contents = aws_utils.list_all_subdirs(aws_config, args.dir_path)

        for src_path in contents['FILES']:
            aws_utils.aws_rm(aws_config, src_path)

if __name__ == '__main__':
    remove_aws_dir()