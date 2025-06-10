import os
import sys
import subprocess
import aws_cloud_config

def list_aws_dir(aws_config,aws_dir_path):

    if aws_dir_path[-1:] != '/':
        aws_dir_path = aws_dir_path+'/'
    if aws_dir_path[0:1] == '/':
        aws_dir_path = aws_dir_path[1:]

    aws_path = os.path.join(aws_config.bucket, aws_dir_path)
    print('PATH: ', aws_path, aws_config.bucket, aws_dir_path, type(aws_config.bucket))

    cl = ['aws', '--profile='+aws_config.profile,'s3','ls',aws_path]
    print(cl)

    p = subprocess.Popen(cl, stdout=subprocess.PIPE)
    #p.wait()
    stdoutput = p.communicate()[0].decode("utf-8")

    contents = {'FILES': [], 'DIRECTORIES': []}
    for line in stdoutput.split('\n'):
        entry = line.lstrip()
        if len(entry) > 0 and 'PRE ' in line and entry != 'PRE /':
            contents['DIRECTORIES'].append(os.path.join(aws_dir_path,entry.replace('PRE ','')))
        elif len(entry) > 0:
            contents['FILES'].append(os.path.join(aws_dir_path,entry.split()[-1]))

    return contents

def list_all_subdirs(aws_config, aws_top_dir):

    contents = list_aws_dir(aws_config,aws_top_dir)
    print('Top level contents: ',contents)
    print('Identifing files in all sub-directories of this tree...')

    n_loop = 0
    while len(contents['DIRECTORIES']) > 0:
        initial_dir_list = contents['DIRECTORIES']
        for sub_dir in initial_dir_list:
            sub_contents = list_aws_dir(aws_config,sub_dir)
            contents = combine_contents(contents, sub_contents)
            idx = contents['DIRECTORIES'].index(sub_dir)
            _ = contents['DIRECTORIES'].pop(idx)
        n_loop += 1

    return contents

def combine_contents(contents, sub_contents):

    for cat in ['FILES', 'DIRECTORIES']:
        for entry in sub_contents[cat]:
            contents[cat].append(entry)

    return contents

def aws_move(aws_config, src_path, dest_path):

    aws_src_path = os.path.join(aws_config.bucket, src_path)
    aws_dest_path = aws_src_path.replace(src_path, dest_path)

    cl = ['aws', '--profile='+aws_config.profile,'s3','mv',aws_src_path, aws_dest_path]

    p = subprocess.Popen(cl, stdout=subprocess.PIPE)
    #p.wait()
    stdoutput = p.communicate()[0].decode("utf-8")
    print(stdoutput)

if __name__ == '__main__':

    aws_dir_path = 'ROMEREA/reduced_data/fields/ROME-FIELD-01/'
    #aws_dir_path = 'ROMEREA/reduced_data/fields/ROME-FIELD-01/ROME-FIELD-01_coj-doma-1m0-11-fa12_gp/'
    aws_config = aws_cloud_config.get_aws_config()
    contents = list_all_subdirs(aws_config,aws_dir_path)

    print('Complete file list: ')
    print(contents)
