# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 14:56:22 2019

@author: rstreet
"""

from os import getcwd, path
from sys import path as systempath
cwd = getcwd()
systempath.append(path.join(cwd,'../'))
import numpy as np
import lens_properties

def test_lens():
    
    lens = lens_properties.Lens()
    
    lens.D = 2.8        # Kpc
    lens.sig_D = 0.4
    
    lens.ML = 0.17
    lens.sig_ML = 0.03
    
    thetaS = 9.143
    sig_thetaS = 0.792

    lens.rho = 0.0165
    lens.sig_rho = 0.0001
    
    lens.tE = 44.3
    lens.sig_tE = 0.1
    
    pi_E_N = 0.37
    sig_pi_E_N = 0.05
    pi_E_E = 0.01
    sig_pi_E_E = 0.01
    lens.pi_E = np.array( [pi_E_N, pi_E_E] )
    lens.sig_pi_E = np.array( [sig_pi_E_N, sig_pi_E_E] )
    
    lens.s = 0.7750
    lens.sig_s = 0.0007
    
    lens.dsdt = 0.49 / 365.24           # d^-1
    lens.sig_dsdt = 0.02 / 365.24 
    lens.dalphadt = -0.37 / 365.24       # d^-1
    lens.sig_dalphadt = 0.08 / 365.24 
    
    return lens

def test_calc_angular_einstein_radius():
    
    lens = test_lens()
    
    thetaS = 9.143  # micro-as
    sig_thetaS = 0.792
    
    lens.calc_angular_einstein_radius(thetaS, sig_thetaS)
    
    print('Angular Einstein radius = '+str(lens.thetaE)+'+/-'+\
                                        str(lens.sig_thetaE)+' microarcsec')
    print('Lens-source relative parallax: ' + str(lens.pi_rel)+\
                                '+/-'+str(lens.sig_pi_rel)+' microarcsec')
                                
    assert round(lens.thetaE,2) == 554.12

def test_calc_distance():
    
    lens = test_lens()
    
    DS = 7.48 # Kpc
    sig_DS = 0.0 # Kpc
    thetaS = 9.143
    sig_thetaS = 0.792
    
    lens.calc_angular_einstein_radius(thetaS, sig_thetaS)
    lens.calc_distance(DS,sig_DS,log=None)
    
    print('Distance to the lens: '+str(lens.D)+' +/- '+\
                                         str(lens.sig_D)+' kpc')
    
    assert round(lens.D,2) == 2.95
    assert round(lens.sig_D,2) == 0.29
    
def test_calc_projected_separation():
    
    lens = test_lens()
    
    DS = 7.48 # Kpc
    sig_DS = 0.0 # Kpc
    thetaS = 9.143
    sig_thetaS = 0.792
    
    lens.calc_angular_einstein_radius(thetaS, sig_thetaS)
    lens.calc_distance(DS,sig_DS,log=None)
    lens.calc_projected_separation(log=None)
    
    print('Projected separation of lens masses = '+\
                    str(lens.a_proj)+' +/- '+str(lens.sig_a_proj)+' AU')
    
    assert round(lens.a_proj,2) == 1.27
    assert round(lens.sig_a_proj,2) == 0.17

def test_calc_einstein_radius():
    
    lens = test_lens()
    DS = 7.48 # Kpc
    sig_DS = 0.0 # Kpc
    thetaS = 9.143
    sig_thetaS = 0.792
    
    lens.calc_angular_einstein_radius(thetaS, sig_thetaS)
    lens.calc_einstein_radius(log=None)
    
    print('Einstein radius = '+str(lens.RE)+' +/- '+str(lens.sig_RE)+' AU')
    
    assert round(lens.RE,2) == 1.55
    assert round(lens.sig_RE,2) == 0.26
    
def test_calc_orbital_energies():
    
    lens = test_lens()
    
    DS = 7.48 # Kpc
    sig_DS = 0.0 # Kpc
    thetaS = 9.143
    sig_thetaS = 0.792
    
    lens.calc_angular_einstein_radius(thetaS, sig_thetaS)
    lens.calc_distance(DS,sig_DS,log=None)
    lens.calc_einstein_radius(log=None)
    lens.calc_orbital_energies(log=None)
    
    print('Binary lens ratio of KE/PE = '+\
                    str(lens.kepe)+' +/- '+str(lens.sig_kepe))
    
    assert round(lens.kepe,2) == 0.08
    assert round(lens.sig_kepe,2) == 0.05
    
if __name__ == '__main__':
    
    test_calc_angular_einstein_radius()
    print('\n')
    test_calc_distance()
    print('\n')
    test_calc_projected_separation()
    print('\n')
    test_calc_einstein_radius()
    print('\n')
    test_calc_orbital_energies()
    