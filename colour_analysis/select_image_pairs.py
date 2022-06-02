# -*- coding: utf-8 -*-
"""
@author: rstreet
"""

from os import path
from sys import argv
from pyDANDIA import metadata
from astropy.table import Table
import astropy.units as u
import numpy as np

class BiColourDataset:
    """Class describing the attributes and methods of a combined dataset
    consisting of pairs of images taken in 2 filters"""

    def __init__(self):

        self.f1dir = None
        self.f1 = None
        self.f1_images = []
        self.f2dir = None
        self.f2 = None
        self.f2_images = []
        self.outdir = None
        self.image_pairs = []
        self.fwhm_thresh = 2.0
        self.sky_thresh = 5000.0
        self.use_fwhm_threshold = True

    def read_meta_data(self):
        """Function to read in the meta data for each dataset"""

        for d in [ 'f1dir', 'f2dir' ]:

            m = metadata.MetaData()

            m.load_a_layer_from_file( getattr(self,d),
                                      'pyDANDIA_metadata.fits',
                                      'images_stats' )

            t = Table()
            t['im_name'] = m.images_stats[1]['IM_NAME']
            t['fwhm'] = m.images_stats[1]['FWHM']
            t['sky'] = m.images_stats[1]['SKY']

            setattr(self, d+'_stats', t)

    def make_image_table(self):
        """Function to combine the image statistics from all datasets into
        a single table, sorted into alphanumerical order"""

        names = list(self.f1dir_stats['im_name']) + \
                        list(self.f2dir_stats['im_name'])

        f = [self.f1]*len(list(self.f1dir_stats['im_name'])) + \
            [self.f2]*len(list(self.f2dir_stats['im_name']))

        fwhm = list(self.f1dir_stats['fwhm'])+ \
                        list(self.f2dir_stats['fwhm'])

        sky = list(self.f1dir_stats['sky']) + \
                        list(self.f2dir_stats['sky'])

        qc = ['0']*len(names)

        image_table = np.array(list(zip(names,f,fwhm,sky,qc)))

        idx = np.argsort(image_table[:,0])

        for i in range(0,4,1):

            image_table[:,i] = image_table[idx,i]

        self.image_table = image_table

    def quality_control(self):
        """Function to apply quality control selection to the images"""

        fwhm = np.zeros(len(self.image_table))
        sky = np.zeros(len(self.image_table))

        for i in range(0,len(self.image_table),1):

            if np.isnan(float(self.image_table[i,2])) == False and \
                np.isnan(float(self.image_table[i,3])) == False:

                fwhm[i] = float(self.image_table[i,2])
                sky[i] = float(self.image_table[i,3])

            else:

                fwhm[i] = 99.9
                sky[i] = 9999.9

        jdx = np.where(fwhm <= self.fwhm_thresh)
        ldx = np.where(sky <= self.sky_thresh)

        if self.use_fwhm_threshold:
            idx = list((set(jdx[0]).intersection(set(ldx[0]))))
        else:
            idx = list(jdx[0])

        self.image_table[idx,4] = '1'

    def append_image_list(self,name,f):

        if f == self.f1:
            self.f1_images.append(name)
        elif f == self.f2:
            self.f2_images.append(name)

    def identify_image_pairs(self,txt_output=False):
        """Function to review the sorted image table to identify pairs of
        images in different passbands that were taken sequentially in the
        same observation subrequest."""

        if txt_output:
            output = open(path.join(self.outdir, 'image_pairs.txt'),'w')
            output.write('# Image name   Filter  FWHM   SKY   QC\n')

        t = Table()

        for i in range(0,len(self.image_table[:,0])-2,1):

            name1 = self.image_table[i,0]
            name2 = self.image_table[i+1,0]
            date1 = int(str(self.image_table[i,0]).split('-')[2])
            date2 = int(str(self.image_table[i+1,0]).split('-')[2])
            num1 = int(str(self.image_table[i,0]).split('-')[3])
            num2 = int(str(self.image_table[i+1,0]).split('-')[3])
            f1 = self.image_table[i,1]
            f2 = self.image_table[i+1,1]
            qc1 = self.image_table[i,4]
            qc2 = self.image_table[i+1,4]

            if date1 == date2 and \
                num2 == (num1 + 1) and \
                f1 != f2 and \
                qc1 == '1' and qc2 == '1' :

                self.image_pairs.append( (i,i+1) )

                self.append_image_list(name1,f1)
                self.append_image_list(name2,f2,)

                if txt_output:
                    for j in range(i,i+2,1):

                        text = ''

                        for item in self.image_table[j,:]:

                            text += ' '+str(item)

                        output.write(text+'\n')

                    output.write('\n')

        if txt_output:
            output.close()

        t['f1_images'] = self.f1_images
        t['f2_images'] = self.f2_images

        self.image_pairs_table = t

def select_image_pairs():
    """Function to select pairs of images in different passbands that were
    taken as a sequence"""

    dataset = get_args()

    dataset.read_meta_data()

    dataset.make_image_table()

    dataset.quality_control()
    for entry in dataset.image_table:
        print(entry)

    dataset.identify_image_pairs(txt_output=True)

def get_args():
    """Function to obtain the necessary parameters"""

    dataset = BiColourDataset()

    if len(argv) < 6:

        dataset.f1dir = input('Please enter the path to the reduction directory for the first filters dataset: ')
        dataset.f1 = input('Please enter the name of the first filter bandpass [e.g. g, r or i]: ')
        dataset.f2dir = input('Please enter the path to the reduction directory for the second filters dataset: ')
        dataset.f2 = input('Please enter the name of the second filter bandpass [e.g. g, r or i]: ')
        dataset.outdir = input('Please enter the path to the output directory: ')

    else:

        dataset.f1dir = argv[1]
        dataset.f1 = argv[2]
        dataset.f2dir = argv[3]
        dataset.f2 = argv[4]
        dataset.outdir = argv[5]

    if len(argv) >= 6:
        for arg in argv[6:]:
            if 'no-fwhm-cut' in arg:
                dataset.use_fwhm_threshold = False

    return dataset

if __name__ == '__main__':

    select_image_pairs()
