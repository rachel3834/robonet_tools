import os
import sys
import aws_cloud_config

def list_aws_dir(aws_dir_path):

    aws_path = os.path.join(aws_config.bucket, aws_dir_path)

    cl = ['aws', '--profile='+aws_config.profile,'s3','cp',local_file_path,aws_path]

    p = subprocess.Popen(cl)
    #p.wait()
    stdoutput = p.communicate()[0]
    if stdoutput:
        print(stdoutput)


if __name__ == '__main__':

    aws_dir_path = 'ROMEREA/reduced_data/fields/'
    list_aws_dir(aws_dir_path)
