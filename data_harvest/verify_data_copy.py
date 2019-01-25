# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 14:07:24 2019

@author: rstreet
"""

from sys import argv
from os import path, walk
import filecmp

def verify_data_copy(orig_dir, copy_dir, log_path):
    
    log = open(log_path,'w')
    
    for dpath, subdirs, files in walk(orig_dir):
        for name in files:
            
            orig = path.join(dpath, name)
            dest = orig.replace(orig_dir,copy_dir)
            dif = filecmp.cmp(orig,dest)
            
            if dif == False:
                log.write(orig +' differs from '+dest+': '+repr(dif))
                
    log.close()
    

if __name__ == '__main__':
    
    if len(argv) == 1:
        orig_dir = raw_input('Please enter the original data top-level directory path: ')
        copy_dir = raw_input('Please enter the data copy top-level directory path: ')
        log_path = raw_input('Please enter the path to the log file: ')
        
    else:
        orig_dir = argv[1]
        copy_dir = argv[2]
        log_path = argv[3]

    verify_data_copy(orig_dir, copy_dir, log_path)