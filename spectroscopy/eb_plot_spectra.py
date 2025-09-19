import numpy as np
import matplotlib.pyplot as plt
import os

os.environ['PYSYN_CDBS'] =  '/home/bachelet/Work/cdbs'


import scipy.optimize as so

import speclite.filters
import speclite

from astropy.io import fits
import astropy.units as u
import astropy.constants as constantes
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from Spyctres import Spyctres


#Event basics

coord = SkyCoord(ra=266.04333*u.degree, dec= -39.11322*u.degree, frame='icrs')
momo = Spyctres.star_spectrum_new(5000,0,4.5,catalog='k93models')
wave = momo._model.points[0]
mask = (wave<10000) & (wave>4000)
wave_ref = wave[mask]


spec = np.loadtxt('./SOAR/wecfsto_03670368_OGLE-2024-BLG-0034_SOAR_spec_1_02-07-2024_target_1_ws_1.csv',skiprows=0)
SOAR1 = np.c_[spec[:,0],spec[:,1],spec[:,1]*0.01]

spec = np.loadtxt('./SOAR/wecfzst_03160317_OGLE-2024-BLG-0034_SOAR_2_23-04-2025_target_1_ws_1.csv',skiprows=0)
SOAR2 = np.c_[spec[:,0],spec[:,1],spec[:,1]*0.01]


#Find telluric lines
telluric_lines,telluric_mask = Spyctres.load_telluric_lines(0.90)

# Plot
plt.yscale('log')
plt.errorbar(SOAR1[:,0], SOAR1[:,1], SOAR1[:,2],fmt='.',label='SOAR1')
plt.errorbar(SOAR2[:,0], SOAR2[:,1], SOAR2[:,2],fmt='.',label='SOAR2')
plt.fill_between(wave_ref,0,10000,where=telluric_mask(wave_ref),color='grey',alpha=0.25)
plt.legend()
plt.show()


breakpoint()
