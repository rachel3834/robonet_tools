# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:52:49 2018

@author: rstreet
"""
import os
import sys
import glob
import subprocess

def compress_pydandia_reduced_dataset(dataset_dir):
    """Function to prepare a dataset process by pyDANDIA for archive storage.

    This function will compress all image data and data products using
    fpack -q 64 compression to avoid losses.

    Input:
        dataset_dir  string  full path to the directory of reduced data
    """

    if not os.path.isdir(dataset_dir):

        print('ERROR: Cannot find the directory of data products '+event_dir)
        sys.exit()

    image_dirs = [ 'data', 'diffim', 'kernel', 'ref', 'resampled' ]

    for d in image_dirs:

        data_path = os.path.join(dataset_dir,d)

        if os.path.isdir(data_path):
            bzip2_image_dir(d)

    if os.path.isfile(os.path.join(dataset_dir,'vphas_catalog.fits')):
        os.remove(os.path.join(dataset_dir,'vphas_catalog.fits'))

def compress_dandia_reduced_data_products(event_dir):
    """Function to prepare a dataset processed by DanDIA for archive storage.

    This function will compress all image data and data products using
    fpack -q 64 compression to avoid losses.

    Input:
        event_dir  string  full path to the directory of reduced data
    """

    if not os.path.isdir(event_dir):

        print('ERROR: Cannot find the directory of data products '+event_dir)
        sys.exit()

    init_dir = os.getcwd()
    os.chdir(event_dir)

    red_filter = identify_reduction_filter(event_dir)

    image_dirs = glob.glob( os.path.join('20??????') )
    image_dirs += glob.glob( os.path.join('imred',red_filter) )
    image_dirs += glob.glob( os.path.join('gimred',red_filter) )
    image_dirs += glob.glob( os.path.join('dimred',red_filter) )

    for d in image_dirs:

        fpack_image_dir(d)

    dir_list = [ os.path.join('imred',red_filter),
                 os.path.join('gimred',red_filter),
                 os.path.join('dimred',red_filter) ]

    for d in dir_list:

        tar_directory(d,event_dir)

    tar_lightcurves(event_dir, red_filter)

    os.chdir(init_dir)
    print('Completed data compression')

def identify_reduction_filter(dir_path):
    """Function to identify which filter was used for the reduction directory"""

    par_files = glob.glob(os.path.join(dir_path, 's2*par'))
    if len(par_files) > 0:
        red_filter = par_files[0].split('.')[1]
    else:
        raise IOError('Cannot find the filter for this reduction directory')
        sys.exit()

    return red_filter

def fpack_image_dir(dir_path):
    """Function to fpack compress all FITS images within a specified directory"""

    file_list = glob.glob( os.path.join(dir_path,'*.fits') )

    for f in file_list:

        if os.path.isfile(f+'.fz') == False:

            args = ['fpack', '-q', '64', f]

            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            p.wait()

            if 'error' in p.stout and not os.path.isfile(f+'.fz'):
                args = ['bzip2', f]

                p = subprocess.Popen(args, stdout=subprocess.PIPE)
                p.wait()

            if os.path.isfile(f+'.fz') or os.path.isfile(f+'.bz2') :
                os.remove(f)
            else:
                print('WARNING: Cannot find compressed data product '+f+'.fz or .bz2, so skipping delete of original')

        else:
            print('Skipping compression of '+f+' (compressed product already exists)')

def bzip2_image_dir(dir_path):
    """Function to bzip2 compress all FITS images within a specified directory"""

    file_list = glob.glob( os.path.join(dir_path,'*.fits') )
    file_list += glob.glob( os.path.join(dir_path,'*.fts') )
    file_list += glob.glob( os.path.join(dir_path,'*.fit') )

    print(' -> Found '+str(len(file_list))+' fits files in '+dir_path)

    if len(file_list) > 0:
        for f in file_list:

            if os.path.isfile(f+'.bz2') == False:

                args = ['bzip2', f]

                p = subprocess.Popen(args, stdout=subprocess.PIPE)
                p.wait()

            else:
                print('Skipping compression of '+f+' (compressed product already exists)')

        print(' -> Completed compression of fits files in '+dir_path)

def bunzip2_dir(dir_path):
    """Function to uncompress all files compressed with bzip2 within a
    specified directory"""

    file_list = glob.glob( os.path.join(dir_path,'*.bz2') )

    print(' -> Found '+str(len(file_list))+' compressed files in '+dir_path)

    if len(file_list) > 0:
        for f in file_list:

            args = ['bunzip2', f]

            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            p.wait()

        print(' -> Completed decompression of files in '+dir_path)

def tar_directory(dir_path,event_dir):
    """Function to build a tarball of all files within a given directory"""

    subdirs = dir_path.split('/')
    subdir = subdirs[-2]+'_'+subdirs[-1]

    tarball = os.path.join(event_dir,subdir+'_data.tar')

    if os.path.isfile(tarball) == False:
        flist = glob.glob( os.path.join(dir_path,'*') )

        args = ['tar', '-cvf',tarball] + flist
        print('Building tarball of '+dir_path)

        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        p.wait()

        print(tarball)

def tar_lightcurves(event_dir, red_filter):
    """Function to build a tarball of the lightcurves within a given
    directory.  This requires special handling because the number of lightcurves
    is often too long for the tar command line."""

    tarball = os.path.join(event_dir,'lightcurves.tar')
    dir_path = os.path.join('lc',red_filter,'rawlc')
    lc_list = os.path.join(dir_path,'..','lc_list.txt')

    if os.path.isfile(tarball) == False:
        flist = glob.glob( os.path.join(dir_path,'*') )

        f = open(lc_list,'w')
        for file in flist:
            f.write(file+'\n')
        f.close()

        cl = ['tar', '-cvf',tarball, '-T', lc_list]

        print(' '.join(cl))
        p = subprocess.Popen(cl)
        stdoutput = p.communicate()[0]

        if stdoutput:
            print(stdoutput)

def prepare_pydandia_reduced_dataset_for_archive(dir_path,output_dir):
    """Function to build a tarball of data products from a pyDANDIA reduction
    of a single dataset ready for archiving.  This will not include the raw
    data files."""

    build_tarball = False

    sub_dir_list = [ 'diffim', 'kernel', 'ref', 'resampled' ]

    data_list = ["*.dat", "*.log", "pyDANDIA_metadata.fits*", "*.png",
                 "diffim/*", "kernel/*", "ref/*", "resampled/*"]

    start_dir = os.getcwd()

    if os.path.isdir(dir_path):

        os.chdir(os.path.join(dir_path,'..'))

        if os.path.isdir(output_dir):

            for d in sub_dir_list:

                subd = os.path.join(os.path.basename(dir_path),d)

                if os.path.isdir(subd):
                    bzip2_image_dir(subd)

            if build_tarball:
                tar_file = os.path.join(output_dir,os.path.basename(dir_path)+'.tar')

                args = ['tar', '-cvf', tar_file]

                for entry in data_list:

                    files = glob.glob(os.path.join(os.path.basename(dir_path),entry))

                    args += files

                p = subprocess.Popen(args, stdout=subprocess.PIPE)
                p.wait()

        else:
            raise IOError('Cannot find output directory '+output_dir)

        os.chdir(start_dir)

    else:
        raise IOError('Cannot find input directory '+dir_path)

def uncompress_pydandia_archived_dataset(dir_path):
    """Function to decompress a pyDANDIA-reduced dataset that has been archived"""

    for sub_dir in os.walk(dir_path):

        print('Decompressing '+sub_dir[0]+'...')

        bunzip2_dir(sub_dir[0])

def get_args():
    """Function to acquire the necessary commandline arguments"""

    if len(sys.argv) == 1:

        dir_path = input('Please enter the path to the directory of a reduced dataset:')
        output_dir = input('Please enter the directory path for output: ')

    else:

        dir_path = sys.argv[1]
        output_dir = sys.argv[2]

    return dir_path, output_dir

if __name__ == '__main__':

    (dir_path,output_dir) = get_args()

    compress_dandia_reduced_data_products(dir_path)
    #prepare_pydandia_reduced_dataset_for_archive(dir_path,output_dir)
    #uncompress_pydandia_archived_dataset(dir_path)
