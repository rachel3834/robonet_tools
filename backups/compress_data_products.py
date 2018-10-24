# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:52:49 2018

@author: rstreet
"""
import os
import sys
import glob
from commands import getstatusoutput

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
    
    
    dir_list += glob.glob( os.path.join(event_dir,'imred','*') )
    dir_list += glob.glob( os.path.join(event_dir,'gimred','*') )
    dir_list += glob.glob( os.path.join(event_dir,'dimred','*') )
    dir_list = glob.glob( os.path.join(event_dir,'lc/*/rawlc/') )
    
    for d in dir_list:
        
        tar_directory(d)
        
    print('Completed data compression')
    
def fpack_image_dir(dir_path):
    """Function to fpack compress all FITS images within a specified directory"""
    
    file_list = glob.glob( os.path.join(dir_path,'*.fits') )
    
    for f in file_list:
        
        if os.path.isfile(f+'.fz') == False:
            
            (iexec,output) = getstatusoutput('fpack -q 64 '+f)
    
            if len(output.replace(' ','').replace('\n','')) > 0:
                
                print(output)
            
            else:
                
                if os.path.isfile(f+'.fz'):
                
                    os.remove(f)
            
        else:
            print('Skipping compression of '+f+' (compressed product already exists)')

def tar_directory(dir_path):
    """Function to build a tarball of all files within a given directory"""
    
    (p,subd) = os.path.split(d)
    subdir = os.path.basename(p)
    
    if 'lc/' in d:
        tarball = os.path.join(d,'lightcurves.tar')
    else:
        tarball = os.path.join(d,subdir+'_data.tar')
    
    flist = glob.glob( os.path.join(d,'*') )
    
    cl = 'tar -cvf '+tarball+' '+flist
    
    (iexec,output) = getstatusoutput(cl)
    
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