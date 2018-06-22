# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 12:45:16 2018

@author: rstreet
"""

from os import path, walk
from sys import argv
from commands import getstatusoutput

def spyder():
    """Function to walk through a directory tree looking for FITS images.
    If found, all FITS images are compressed using fpack -q 64"""
    
    (top_dir, uncompress) = get_args()
    
    sub_sub_dirs = []
    file_list = []
    
    for d, subd, files in walk(top_dir):
        
        sub_sub_dirs += subd
        for f in files:
            file_list.append(path.join(d,f))
                            
    compress_fits_images(file_list, top_dir, uncompress=uncompress)

            
def compress_fits_images(file_list, dir_path, uncompress=False):
    """Function to compress any FITS images in a list of files of various types"""
    
    fcount = 0
    
    for f in file_list:
        
        if uncompress == False:
            
            if f.endswith('.fits') or f.endswith('.fts'):
                
                fcount += 1
                
                #(iexec,output) = getstatusoutput('fpack -q 64 '+f)
                output = ''
                 
                if len(output.replace(' ','').replace('\n','')) > 0:
                    
                    print(output)
                
                else:
                    
                    if path.isfile(f+'.fz') and path.isfile(f):
                    
                        os.remove(f)
            
        else:
            
            if f.endswith('.fz'):
                
                fcount += 1
                
                #(iexec,output) = getstatusoutput('funpack '+f)
                
                output = ''
                
                if len(output.replace(' ','').replace('\n','')) > 0:
                    
                    print(output)
                
    if uncompress:
        
        print('Uncompressed '+str(fcount)+' FITS files in '+dir_path)
        
    else:
        
        print('Compressed '+str(fcount)+' FITS files in '+dir_path)

def get_args():
    """Function to harvest commandline arguments or ask the user"""
    
    uncompress = False
    
    if len(argv) < 3:
        
        top_dir = raw_input('Please enter the top level path to the directory tree for image compression: ')
        uncompress = raw_input('Compress (Y) or uncompress (N) the data? ')
        
        if 'Y' in str(uncompress).upper():

            uncompress = True
            
    else:
        
        top_dir = argv[1]
        
        if 'uncompress' in argv[2]:
            uncompress = True
        
    return top_dir, uncompress
    
if __name__ == '__main__':
    
    spyder()
    