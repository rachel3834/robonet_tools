# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 14:31:41 2018

@author: rstreet
"""

from os import getcwd, path
from sys import exit
from sys import path as systempath
cwd = getcwd()
systempath.append(path.join(cwd,'../'))
import make_colour_image
import numpy as np
from astropy.table import Table

def test_intialise_colour_image():
    """Function to verify the instantiation of a ColourImage object
    and check that the input parameters are correctly set"""
    
    cmd_file = 'make_colour_image.cmd'
    
    test_image = make_colour_image.ColourImage()
    
    cimage = make_colour_image.intialise_colour_image(cmd_file)
    
    assert type(cimage) == type(test_image)
    
    keys = ['b_image_path', 'g_image_path', 'r_image_path', \
            'b_cat_path', 'g_cat_path', 'r_cat_path']
    
    for key in keys:
        assert getattr(cimage,key) != None


def test_read_image_data():
    """Function to test the class method to read in the image data"""
    
    cmd_file = 'make_colour_image.cmd'

    cimage = make_colour_image.intialise_colour_image(cmd_file)
    
    cimage.read_image_data()
    
    for par in ['b_image', 'g_image', 'r_image']:
        
        assert type(getattr(cimage, par)) != None
        assert type(getattr(cimage, par)) == type(np.zeros(1))

def test_read_star_catalogs():
    """Function to test the class method to read in the star catalog data
    from pyDANDIA meta datafiles"""
    
    cmd_file = 'make_colour_image.cmd'

    cimage = make_colour_image.intialise_colour_image(cmd_file)
    
    cimage.read_star_catalogs()
    
    test_table = Table()
    
    for par in ['b_cat', 'g_cat', 'r_cat']:
        
        cat = getattr(cimage, par)
        
        assert type(cat) == type(test_table)
        assert 'star_index' in cat.colnames
        assert 'x' in cat.colnames
        assert 'y' in cat.colnames
        assert 'RA' in cat.colnames
        assert 'DEC' in cat.colnames
        assert len(cat) > 0

def test_calc_offsets_from_red():
    """Function to test the class method which calculates the X,Y shift
    offset of the blue and green catalogues relative to the red. """
    
    cmd_file = 'make_colour_image.cmd'

    cimage = make_colour_image.intialise_colour_image(cmd_file)
    
    cimage.read_star_catalogs()

    cimage.calc_offsets_from_red()
    
    assert cimage.b_offset != None
    assert cimage.g_offset != None

def test_stack_colour_layers():
    """Function to test the class method which combines the 3 images into
    a single, aligned array"""

    cimage = make_colour_image.ColourImage()
    
    cimage.out_image_path = 'test_colour_image.png'
    
    cimage.r_image = np.zeros([100,100])
    cimage.r_image[40:60,40:60] = 1.0
    
    cimage.g_image = np.zeros([100,100])
    cimage.g_image[35:55,35:55] = 1.0
    cimage.g_offset = (-5,-5)
    
    cimage.b_image = np.zeros([100,100])
    cimage.b_image[45:65,45:65] = 1.0
    cimage.b_offset = (-5,-5)
    
    cimage.stack_colour_layers()
    
    assert type(cimage.colour_image) == type(np.zeros(1))

    cimage.output_colour_image()
    
def test_output_colour_image():
    """Function to test the class method to output a colour image in PNG format"""
    
    cimage = make_colour_image.ColourImage()
    
    cimage.out_image_path = 'test_colour_image.png'
    
    cimage.colour_image = np.zeros([100,100,3])
    
    cimage.colour_image[45:55,45:55,0:2] = 1.0
    
    cimage.output_colour_image()
    
    
if __name__ == '__main__':
    
    #test_intialise_colour_image()
    #test_read_image_data()
    #test_read_star_catalogs()
    #test_calc_offsets_from_red()
    test_stack_colour_layers()
    #test_output_colour_image()
    