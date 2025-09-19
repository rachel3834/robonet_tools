import numpy as np
import matplotlib.pyplot as plt
import os

os.environ['PYSYN_CDBS'] =  '/home/bachelet/Work/cdbs/'


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

def objective_function(params,spectrums,mask,catalog,bounds):
    try:
        for ind,bound in enumerate(bounds):

            if  (params[ind]<bound[0]) | (params[ind]>bound[1]):
                return -np.inf

        chichi = Spyctres.fit_spectra_chichi(params,spectrums,mask,catalog)

        return -chichi

    except:
        return -np.inf




#Initializing
Spyctres.define_2MASS_filters()
Spyctres.define_GAIA_filters()

twomass_filters = speclite.filters.load_filters('MASS-J', 'MASS-H','MASS-K')
sdss_filters = speclite.filters.load_filters('sdss2010-*')
bessel_filters = speclite.filters.load_filters('bessell-*')

gaia_filters = speclite.filters.load_filters('gaiadr2-*')



#Event basics

coord = SkyCoord(ra=277.56025*u.degree, dec=-8.22021*u.degree, frame='icrs')
momo = Spyctres.star_spectrum_new(5000,0,4.5,catalog='k93models')
wave = momo._model.points[0]
mask = (wave<24500) & (wave>4000)
wave_ref = wave[mask]

seed = [ 1.04065605e+00,  8.31050919e+00, -3.10593017e+02,  3.63938949e+00,
        3.74264615e-05,  3.15651712e-01,  2.83040782e-01,  1.94752823e-01,
        1.85565475e+00,  1.87596383e+00]

#Adding new spectra



spectra = {}
SEDS=[]
ABS = []
magnifications = []


spec = './Spectra/Gaia18ajz_NIR_20180325_IDP.fits'
NIR = np.c_[[fits.open(spec)[1].data['WAVE'].astype(float)[0],fits.open(spec)[1].data['FLUX'][0].astype(float),fits.open(spec)[1].data['ERR'][0].astype(float)]].T

spec = './Spectra/Gaia18ajz_VIS_20180325_IDP.fits'
VIS =np.c_[[fits.open(spec)[1].data['WAVE'].astype(float)[0],fits.open(spec)[1].data['FLUX'][0].astype(float),fits.open(spec)[1].data['ERR'][0].astype(float)]].T

spec = './Spectra/Gaia18ajz_UVB_20180325_IDP.fits'
UVB = np.c_[[fits.open(spec)[1].data['WAVE'].astype(float)[0],fits.open(spec)[1].data['FLUX'][0].astype(float),fits.open(spec)[1].data['ERR'][0].astype(float)]].T



spectrum_XShooter = np.r_[UVB,VIS,NIR]
spectrum_XShooter[:,0] *= 10#nm to Angstrom


bin_spec,cov = Spyctres.bin_spectrum(spectrum_XShooter,wave_ref)
SNR = bin_spec[:,1]/bin_spec[:,2]
#offset1 = 10**(offset1(bin_spec[:,0])/2.5)
spectrum_XShooter  = np.c_[bin_spec[:,0],bin_spec[:,1],bin_spec[:,1]/SNR]
#breakpoint()

spectrum_XShooter = spectrum_XShooter[spectrum_XShooter[:,0].argsort(),]
spectrum_XShooter = spectrum_XShooter[np.unique(spectrum_XShooter[:,0],return_index=True)[1]]

mask = (spectrum_XShooter[:,1]>0) & (spectrum_XShooter[:,2]>0) #& (spectrum_XShooter[:,1]/spectrum_XShooter[:,2]>5) & (spectrum_XShooter[:,1]<10**-14)
spectrum_XShooter = spectrum_XShooter[mask]


spectra['XShooter_25_03_2018'] = {}
spectra['XShooter_25_03_2018']['spectrum'] = spectrum_XShooter
spectra['XShooter_25_03_2018']['JD'] = Time(2458202.86309,format='jd')
spectra['XShooter_25_03_2018']['magnification'] = 4.93
spectra['XShooter_25_03_2018']['barycentric_velocity']  = Spyctres.Barycentric_velocity(spectra['XShooter_25_03_2018']['JD'],coord)
spectra['XShooter_25_03_2018']['SED'] = [[gaia_filters[1],17.69,0.1],[bessel_filters[4],16.257+0.45,0.1]]
offset1 = Spyctres.SED_offset(np.array(spectra['XShooter_25_03_2018']['SED']),spectrum_XShooter)


#Find telluric lines
telluric_lines,telluric_mask = Spyctres.load_telluric_lines(0.90)

# Plot
plt.yscale('log')
plt.errorbar(spectrum_XShooter[:,0], spectrum_XShooter[:,1], spectrum_XShooter[:,2],fmt='.',label='SALT')
plt.fill_between(wave_ref,0,1,where=telluric_mask(wave_ref),color='grey',alpha=0.25)
plt.show()


#Fit
bound = [[-2,2],[4,10],[-2000,2000],[3,4],[-0.9,-0.3],[0.0,3.0],[np.log10(0.5*offset1[0]),np.log10(1.5*offset1[0])]]#,[np.log10(0.5*offset2[0]),np.log10(1.5*offset2[0])],[-2,2],[-2,2]]


#res = so.differential_evolution(Spyctres.fit_spectra_chichi, bound, args=(spectra,telluric_mask,'k93models'),disp=True,popsize=2,workers=4,polish=False)


#momo = Spyctres.model_spectra(res['x'],spectra)
momo = Spyctres.model_spectra(seed,spectra)

plt.yscale('log')
plt.errorbar(spectrum_XShooter[:,0], spectrum_XShooter[:,1], spectrum_XShooter[:,2],fmt='.',label='SALT')
plt.fill_between(wave_ref,0,1,where=telluric_mask(wave_ref),color='grey',alpha=0.25)
plt.plot(momo[0][:,0],momo[0][:,1])
plt.xlabel(r'$\lambda [\AA]$')
plt.ylabel(r'$F_\lambda [erg/s/cm^2/\AA]$')
plt.show()


breakpoint()

import emcee
import multiprocessing as mul
nwalkers = 2*len(seed)
ndim = len(seed)
nchains = 10000

objective_function(seed,spectra,telluric_mask,'k93models',bound)

with mul.Pool(processes=8) as pool:


    pos = seed +  len(seed)*[1] * np.random.randn(nwalkers, len(seed))*10**-3
            
    nwalkers, ndim = pos.shape    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, objective_function, args=(spectra,telluric_mask,'k93models',bound),pool = pool)
    final_positions, final_probabilities, state = sampler.run_mcmc(pos, nchains, progress=True)

    breakpoint()


