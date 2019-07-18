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
    
    image_dirs = glob.glob( os.path.join(event_dir,'20??????') )
    image_dirs += glob.glob( os.path.join(event_dir,'imred','*') )
    image_dirs += glob.glob( os.path.join(event_dir,'gimred','*') )
    image_dirs += glob.glob( os.path.join(event_dir,'dimred','*') )
    
    for d in image_dirs:
        
        fpack_image_dir(d)
    
    
    dir_list = [ os.path.join(event_dir,'imred'),
                 os.path.join(event_dir,'gimred'),
                 os.path.join(event_dir,'dimred'),
                 os.path.join(event_dir,'lc/*/rawlc/') ]
    
    for d in dir_list:
        
        tar_directory(d,event_dir)
        
    print('Completed data compression')
    
def fpack_image_dir(dir_path):
    """Function to fpack compress all FITS images within a specified directory"""
    
    file_list = glob.glob( os.path.join(dir_path,'*.fits') )
    
    for f in file_list:
        
        if os.path.isfile(f+'.fz') == False:
                
            args = ['fpack', '-q', '64', f]
            
            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            p.wait()
            
            if os.path.isfile(f+'.fz'):
                os.remove(f)
            else:
                print('WARNING: Cannot find compressed data product '+f+'.fz, so skipping delete of original')
                
        else:
            print('Skipping compression of '+f+' (compressed product already exists)')

def bzip2_image_dir(dir_path):
    """Function to bzip2 compress all FITS images within a specified directory"""
    
    file_list = glob.glob( os.path.join(dir_path,'*.fits') )
    
    for f in file_list:
        
        if os.path.isfile(f+'.fz') == False:
                
            args = ['bzip2', f]
            
            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            p.wait()
            
        else:
            print('Skipping compression of '+f+' (compressed product already exists)')
    
    print(' -> Completed compression of fits files in '+dir_path)
    
def tar_directory(dir_path,event_dir):
    """Function to build a tarball of all files within a given directory"""
    
    subdir = os.path.basename(dir_path)

    if 'lc/' in dir_path:
        tarball = os.path.join(event_dir,'lightcurves.tar')
    else:
        tarball = os.path.join(event_dir,subdir+'_data.tar')

    if os.path.isfile(tarball) == False:
        flist = glob.glob( os.path.join(dir_path,'*') )

        cl = 'tar -cvf '+tarball+' '+dir_path

        print('Building tarball of '+dir_path)

        (iexec,coutput) = getstatusoutput(cl)

        print(tarball)
        print(coutput)
    
def get_args():
    """Function to acquire the necessary commandline arguments"""
    
    if len(sys.argv) == 1:
        
        event_dir = raw_input('Please enter the path to the directory of reduced data [DanDIA-format]:')
        
    else:
        
        event_dir = sys.argv[1]
    
    return event_dir
    
if __name__ == '__main__':
    
    event_dir = get_args()    
    
    compress_dandia_reduced_data_products(event_dir)