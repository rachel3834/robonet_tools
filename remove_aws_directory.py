
import os
import sys
import subprocess
import aws_cloud_config

def remove_directory(aws_dir_path):
    """Function to delete a 'directory' and all it's contents in AWS.
    Note that since AWS entries are technically not directories in the traditional sense,
    the only way to achieve this is is to delete all entries within the directory.
    """

    aws_dir_path = verify_path_slashes(aws_dir_path)

    aws_config = aws_cloud_config.get_aws_config()



def fetch_aws_entry_list(aws_config, aws_dir_path):
    """Function to fetch a list of the contents of the AWS 'directory'"""

    aws_path = os.path.join(aws_config.bucket, aws_dir_path)

    cl_args_list = ['aws', '--profile='+aws_config.profile,'s3','ls',aws_path]

    aws_output = run_aws_command(aws_config, cl_args_list)
    
def run_aws_command(aws_config, cl_args_list):
    """Function to run an AWS command as a subprocess, when provided with the
    necessary arguments and configuration."""

    # WARNING: Cannot use stdout-> pipe as the cl output of awscli for larger
    # files can deadline the child process.
    p = subprocess.Popen(cl)
    #p.wait()
    stdoutput = p.communicate()[0]
    if stdoutput:
        print(stdoutput)

    return stdoutput

def verify_path_slashes(dir_path):

    if not dir_path[0:1] == '/':
        dir_path = '/' + dir_path

    if not dir_path[-1:] == '/':
        dir_path = dir_path + '/'

    return dir_path

if __name__ == '__main__':
    if len(sys.argv) == 1:
        aws_dir_path = input('Please enter the path to the "directory" in AWS: ')
    else:
        aws_dir_path = sys.argv[1]

    remove_directory(aws_dir_path)
