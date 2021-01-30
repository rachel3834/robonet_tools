from astropy.io import fits
import os
import shutil  
import subprocess


data_path = '/data/omega/data/data_reduction/'




def decompress_corrupted_files(data_dir):



    list_of_datasets = os.listdir(data_dir)
    
    
    for dataset in list_of_datasets:
    
        path = data_dir+dataset+'/'
        list_of_file = os.listdir(path)
        
        for filename in list_of_file:
        
            decompress_corrupt(path,filename)


def decompress_corrupt(path,filename):

    if '.fits' in filename:
    
        open_fits = fits.open(path+filename)
        
        if type(open_fits[1]).__name__ == 'CompImageHDU':
        
            sourcefile = path+filename
            destfile = path+filename+'.fz'
            shutil.move(sourcefile,destfile)
            
            args = ['funpack', frame]
            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            
            p.wait()
            
            os.remove(destfile)
