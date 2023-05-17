from os import path
import subprocess

class AWSConfig():
    """Configuration class for AWS Cloud"""

    def __init__(self, config_file_path):
        self.profile = None
        self.bucket = None
        self.output_directory = None

        if path.isfile(config_file_path):

            file_lines = open(config_file_path,'r').readlines()

            for line in file_lines:

                data = line.replace('<',':::').replace('>',':::').split(':::')

                if len(data) > 3:

                    data = data[2]

                    for key in ['profile','bucket','output_directory']:

                        if key in line:
                            setattr(self,key,data)



def upload_file(aws_config, local_file_path, local_root, aws_root):
    """Function to upload a single file to an AWS S3 bucket"""

    if aws_root[0:1] == '/':
        aws_root = aws_root[1:]

    aws_path = path.join(aws_config.bucket, local_file_path.replace(local_root,aws_root))

    cl = ['aws', '--profile='+aws_config.profile,'s3','cp',local_file_path,aws_path]

    # WARNING: Cannot use stdout-> pipe as the cl output of awscli for larger
    # files can deadline the child process.
    p = subprocess.Popen(cl)
    #p.wait()
    stdoutput = p.communicate()[0]
    if stdoutput:
        print(stdoutput)
