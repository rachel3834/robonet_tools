# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 14:56:29 2018

@author: rstreet
"""
from sys import argv, exit
from os import path
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import matching
from astropy.table import Table
import astropy.units as u
from pyDANDIA import metadata
import numpy as np
import png
from matplotlib import pyplot as plt


class ColourImage:
    """Class describing the separate images combined to make a colour image"""
    
    def __init__(self):
        
        self.b_image_path = None
        self.b_cat_path = None
        self.g_image_path = None
        self.g_cat_path = None
        self.r_image_path = None
        self.r_cat_path = None
        self.tolerance = 1.0 * u.arcsec
        self.out_image_path = None
    
    def read_image_data(self):
        
        for par in ['b_image', 'g_image', 'r_image']:
            
            image_path = getattr(self, par+'_path')
            
            if image_path != None and path.isfile(image_path):
                
                hdu = fits.open(image_path)
                
                setattr(self, par, hdu[0].data)
                
            else:
                
                print('Error reading in '+image_path)
                
                exit()
    
    def read_star_catalogs(self):
        
        for par in ['b_cat', 'g_cat', 'r_cat']:
            
            meta_file = getattr(self, par+'_path')
            
            if meta_file != None and path.isfile(meta_file):
    
                m = metadata.MetaData()
                m.load_a_layer_from_file( path.dirname(meta_file), 
                                          path.basename(meta_file), 
                                          'star_catalog' )

                star_catalog = Table()
                star_catalog['star_index'] = m.star_catalog[1]['star_index']
                star_catalog['x'] = m.star_catalog[1]['x_pixel']
                star_catalog['y'] = m.star_catalog[1]['y_pixel']
                star_catalog['RA'] = m.star_catalog[1]['RA_J2000']
                star_catalog['DEC'] = m.star_catalog[1]['DEC_J2000']

                setattr(self, par+'_meta', m)
                setattr(self, par, star_catalog)
                
    def calc_offsets_from_red(self):
        """Method to calculate the offset of the blue and green from the red, 
        by matching stars between the star catalogs"""
        
        red_stars = SkyCoord(self.r_cat['RA'], self.r_cat['DEC'], unit="deg")
        
        for key,par in {'b_offset':'b_cat', 'g_offset':'g_cat'}.items():
            
            cat = getattr(self, par)
            
            cat_stars = SkyCoord(cat['RA'], cat['DEC'], unit="deg")
            
            match_data = matching.search_around_sky(red_stars, cat_stars, 
                                                seplimit=self.tolerance)    
            
            
            dx = self.r_cat['x'][match_data[0]] - cat['x'][match_data[1]]
            dy = self.r_cat['y'][match_data[0]] - cat['y'][match_data[1]]
            
            setattr(self, key, (dx.mean(), dy.mean()) )
            print dx.mean(), dy.mean(), np.median(dx), np.median(dy)
    
    def stack_colour_layers(self):
        """Method to apply the calculated offsets and combine the 3 images
        into a single array"""
        
        image = np.zeros([self.r_image.shape[0], self.r_image.shape[1], 3])
        
        image[:,:,0] = self.r_image
        
        x = np.arange(0, len(self.g_image[0]),1) + self.g_offset[0]
        idx = np.where(x > 0)
        x = x[idx]
        idx = np.where(x < len(self.g_image[0]))
        x = x[idx]
        
        print x
        print self.g_offset
        
        y = np.arange(0, len(self.g_image[1]),1) + self.g_offset[1]
        jdx = np.where(y > 0)
        y = y[jdx]
        jdx = np.where(y < len(self.g_image[1]))
        y = y[jdx]
        
        (xx, yy) = np.meshgrid(x, y)
        (xv, yv) = np.meshgrid(idx[0], jdx[0])
        
        image[xx,yy,1] = self.g_image[xv,yv]
        
        x2 = np.arange(0, len(self.b_image[0]),1) + self.b_offset[0]
        idx = np.where(x2 > 0)
        x2 = x2[idx]
        idx = np.where(x2 < len(self.b_image[0]))
        x2 = x2[idx]
        
        print x2
        print self.b_offset
        
        y2 = np.arange(0, len(self.b_image[1]),1) + self.b_offset[1]
        jdx = np.where(y2 > 0)
        y2 = y2[jdx]
        jdx = np.where(y2 < len(self.b_image[1]))
        y2 = y2[jdx]
        
        (xx, yy) = np.meshgrid(x2, y2)
        (xv, yv) = np.meshgrid(idx[0], jdx[0])
        
        image[xx,yy,2] = self.b_image[xv,yv]
        
        self.colour_image = image
        
    def output_colour_image(self):
        """Method to output the colour image in PNG format"""
                
        fig = plt.figure(1,(10,10))
        
        ax=plt.subplot(1,1,1)
        
        plt.imshow(self.colour_image, aspect='equal')

        plt.axis('off')
        
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        
        plt.savefig(self.out_image_path, bbox_inches=extent)
        
def make_colour_movie():
    """Function to combine near-simultaneous images acquired in 3 passbands
    to make a series of colour images of the same field and combine them
    into a movie"""

    # Read in trendlogs for 3 passband datasets of the same field
    # pointing. 

    # Select near-simultaneous trios of images that meet the 
    # quality criteria. 

    # Select a trio to act as the reference frame. 
    # Read the star catalogues for this trio.

    # Calculate offsets between reference trio

    # Transform reference trio images into their common reference frame

    # For all other image trios:
    
        # Read the star catalogues for the trio components. 
    
        # Calculate the offset of all components of the trio from the reference 
        # frame
    
        # Transform all images to the common reference frame and combine
    
        # Store as frame in final movie.
    
    print('Function not yet implemented')

    
def make_colour_image():
    """Function to combine 3 images from different passband into a 
    colourised image"""
    
    cmd_file = get_args()
    
    cimage = intialise_colour_image(cmd_file)
    
    cimage.read_image_data()
    
    cimage.read_star_catalogs()
    
    cimage.calc_offsets_from_red()
    
    # Transform images into a common reference frame
    
    # Colourise and combine into the final image

def get_args():
    """Function to harvest the required arguments from the user, either
    via the commandline arguments or by request.  
    """
    
    if len(argv) == 1:
        
        cmd_file = input('Please enter the path to the input command file: ')
    
    else:
        
        cmd_file = argv[1]
    
    return cmd_file

def intialise_colour_image(cmd_file):
    """Function to extract the required input parameters from a command file
    in ASCII format containing the following parameters in two-column, space 
    or tab-separated format:
    bimage  path to bluest image
    bcat    path to bluest catalogue
    gimage  path to mid-wavelength image
    gcat    path to mid-wavelength catalogue
    rimage  path to redest image
    rcat    path to redest catalogue
    outimage path to output image [PNG]
    The parameter name will be converted to lower case. 
    """
    
    if not path.isfile(cmd_file):
        
        print('ERROR: Cannot find input command file '+cmd_file)
        
        exit()
    
    file_lines = open(cmd_file,'r').readlines()
    
    cimage = ColourImage()
    
    for line in file_lines:
        
        try:
            (param, entry) = line.replace('\n', '').split()
            
            setattr(cimage, str(param).lower(), entry)
        except ValueError:
            
            print('ERROR: Mal-formed command file entry at line:')
            print(line)
            
            exit()
    
    return cimage
    
if __name__ == '__main__':
    
    make_colour_image()
    