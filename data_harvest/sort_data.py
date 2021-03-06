# -*- coding: utf-8 -*-
"""
Created on Thu May 10 09:41:02 2018

@author: rstreet
"""

from sys import argv, exit
from os import path, mkdir
from os import stat as filestatus
from astropy.io import fits
from sys import argv
import glob
from shutil import move

def sort_data(data_dir, option):
    """Function to sort a directory of FITS frames into per-target, per-filter
    sub-directories"""

    image_list = make_image_list(data_dir)

    for image in image_list:

        ds = get_image_dataset(image,option)

        sort_image_to_dataset(image,ds,data_dir)

class Dataset():
    """Class describing a unique dataset from the LCO Network"""

    def __init__(self):
        self.network = None
        self.site = None
        self.enclosure = None
        self.tel = None
        self.instrument = None
        self.filter = None
        self.target = None
        self.id_code = None
        self.facilities = {'Faulkes Telescope South': ['LCOGT','COJ','CLMA','FTS'],
                           'Faulkes Telescope North': ['LCOGT','OGG','CLMA', 'FTN'],
                    	   'Liverpool Telescope': ['LJMU','LAP','LT','LT'] }

    def parse_telescope(self):

        if self.tel in self.facilities.keys():
            (self.network,self.site,self.enclosure,self.tel) = self.facilities[self.tel]

        if self.network == None:
            self.network = 'EXTN'

    def get_dataset_code(self, separate_instruments=False):
        print(separate_instruments)
        if separate_instruments:
            self.id_code = str(self.target)+'_' + \
                    str(self.site).lower()+'-'+\
                    str(self.enclosure).lower()+'-'+\
                    str(self.tel).lower()+'-'+\
                    str(self.instrument).lower()+'_'+\
                    self.filter
        else:
            self.id_code = str(self.target)+ '_' + self.filter

def get_image_dataset(image,option):
    """Function to identify what dataset an image belongs to, based on the
    target, instrument and filter information from its FITS header.

    Inputs:
        :param str image: Full path to image FITS file
    """

    telescopes = {'Faulkes Telescope South': ['LCOGT','COJ','CLMA','FTS'],
                  'Faulkes Telescope North': ['LCOGT','OGG','CLMA', 'FTN'],
        		'Liverpool Telescope': ['LJMU','LAP','LT','LT'] }

    hdr = fits.getheader(image)

    ds = Dataset()
    ds.target = hdr['OBJECT'].replace('/','')
    ds.site = hdr['SITEID'].replace('/','')
    ds.enclosure = hdr['ENCID'].replace('/','')
    ds.tel = hdr['TELESCOP'].replace('/','')
    ds.instrument = hdr['INSTRUME'].replace('/','')
    ds.filter = hdr['FILTER']

    ds.parse_telescope()

    ds.get_dataset_code(separate_instruments=option)

    return ds


def make_image_list(data_dir):
    """Function to check to see if the data directory contains FITS files"""

    if path.isdir(data_dir):

        frame_list = glob.glob(path.join(data_dir,'*.fits'))

        if len(frame_list) == 0:

            print('Cannot find any FITS images in the input directory')
            exit()

    else:

        print('Cannot find input directory: '+data_dir)
        exit()

    return frame_list

def sort_image_to_dataset(image,ds,data_dir):
    """Function to move the given image to a sub-directory determined by
    its dataset ID code"""

    dest_dir = path.join(data_dir,ds.id_code)

    if not path.isdir(dest_dir):
        mkdir(dest_dir)

    move(image,path.join(dest_dir,path.basename(image)))

    print(path.basename(image)+' --> '+path.basename(dest_dir))

if __name__ == '__main__':

    if len(argv) == 1:
        data_dir = input('Please enter data directory path: ')
        option = input('Combine all data for a given target from multiple instruments?  T or F: ')
    else:
        data_dir = argv[1]
        option = argv[2]

    if 't' in str(option).lower():
        option = True
    else:
        option = False

    sort_data(data_dir, option)

    print('Completed data sorting')
