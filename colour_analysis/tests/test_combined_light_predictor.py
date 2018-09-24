# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 16:19:37 2018

@author: rstreet
"""


from os import getcwd, path
from sys import exit
from sys import path as systempath
cwd = getcwd()
systempath.append(path.join(cwd,'../'))
import combined_light_predictor

def setup_test_star():
    
    star = combined_light_predictor.Star()
    
    star.MB = 11.0
    star.MV = 10.0
    star.Mg = 10.0
    star.Mr = 9.5
    star.Mi = 9.2
    star.Mz = 9.1
    star.MJ = 8.2
    star.MH = 8.3
    star.MKs = 8.4
    star.W149 = 8.0
    star.Z087 = 9.1
    
    star.distance_modulus = 11.0
    
    return star

def setup_test_binary_star():
    
    star1 = setup_test_star()
    star2 = setup_test_star()
    
    binary = combined_light_predictor.BinaryStar()
    
    binary.star1 = star1
    binary.star2 = star2
    
    binary.distance_modulus = star1.distance_modulus
    
    return binary
    
def test_star_calculate_colours():
    
    star = setup_test_star()
    
    star.calculate_colours()
    
    assert star.gr == (star.Mg - star.Mr)
    assert star.ri == (star.Mr - star.Mi)
    assert star.gi == (star.Mg - star.Mi)
    assert star.JH == (star.MJ - star.MH)
    assert star.HKs == (star.MH - star.MKs)
    assert star.JKs == (star.MJ - star.MKs)
    assert star.BV == (star.MB - star.MV)
    
def test_star_calculate_apparent_magnitudes():
    
    star = setup_test_star()
    
    star.calculate_apparent_magnitudes()
    
    for f in ['B','V','g','r','i','z','J','H','Ks']:
        
        assert getattr(star,'m'+f) == (getattr(star,'M'+f) + \
                                        getattr(star,'distance_modulus'))

def test_star_apply_extinction():
    
    star = setup_test_star()
    
    star.calculate_colours()
    
    star.calculate_apparent_magnitudes()
        
    extinction = { 'AV': 1.0, 'Ag': 1.0, 'Ar': 1.0, 'Ai': 1.0, 
                  'EBV': 0.1,
                  'Rg': 1.0, 'Rr': 1.0, 'Ri': 1.0 }
    
    star.apply_extinction(extinction)
        
    for f in ['V','g','r','i']:
        
        dmag = getattr(star,'m'+f+'_corr') - getattr(star,'m'+f)
        
        assert dmag == extinction['A'+f]
    
    dcol = round(star.BV_corr - star.BV,1)
    
    assert dcol == extinction['EBV']
    
    c1 = ['g', 'r']
    c2 = ['r', 'i']
        
    for i in range(0,2,1):
        
        mcorr1 = getattr(star,'m'+c1[i]) - extinction['A'+c1[i]]
        mcorr2 = getattr(star,'m'+c2[i]) - extinction['A'+c2[i]]
        
        m1 = getattr(star,'m'+c1[i]+'_corr')
        m2 = getattr(star,'m'+c2[i]+'_corr')
        
        assert (m1-m2) == (mcorr1-mcorr2)

def test_mag_to_flux():
    """Function to test the conversion from magnitudes into flux units"""
    
    mag = 10.0
    ZP = 25.0
    
    flux = combined_light_predictor.mag_to_flux(mag,ZP)
    
    assert flux == 1e6

def test_flux_to_mag():
    
    flux = 1e6
    ZP = 25.0
    
    mag = combined_light_predictor.flux_to_mag(flux,ZP)
    
    assert mag == 10.0
    
def test_binary_calculate_combined_light():
    
    binary = setup_test_binary_star()
    
    binary.calculate_combined_light()
    
    for f in ['B','V','g','r','i','z','J','H','Ks']:
        
        m1 = getattr(binary.star1,'M'+f)
        m2 = getattr(binary.star2,'M'+f)
        
        f1 = combined_light_predictor.mag_to_flux(m1,binary.star1.ZP)
        f2 = combined_light_predictor.mag_to_flux(m2,binary.star2.ZP)
    
        cflux = f1 + f2
    
        cmag = combined_light_predictor.flux_to_mag(cflux,binary.star1.ZP)
    
        assert cmag == getattr(binary,'M'+f+'_combined')
        
def test_binary_calculate_apparent_magnitudes():
    
    binary = setup_test_binary_star()
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    for f in ['B','V','g','r','i','z','J','H','Ks']:
        
        dmag = getattr(binary,'M'+f+'_combined') + \
                            getattr(binary,'distance_modulus')
        
        assert getattr(binary,'m'+f+'_combined') == dmag


def test_binary_calculate_combined_colours():
    
    binary = setup_test_binary_star()
    
    binary.calculate_combined_light()
    
    binary.calculate_combined_colours()
    
    c1 = [ 'B', 'g', 'r', 'J', 'H', 'J' ]
    c2 = [ 'V', 'r', 'i', 'H', 'Ks', 'Ks' ]
        
    for i in range(0,len(c1),1):
        
        m1 = getattr(binary,'M'+c1[i]+'_combined')
        m2 = getattr(binary,'M'+c2[i]+'_combined')
        
        assert getattr(binary,c1[i]+c2[i],(m1-m2))

def test_binary_apply_extinction():
    
    binary = setup_test_binary_star()
    
    binary.calculate_combined_light()
    
    binary.calculate_combined_colours()
    
    binary.calculate_apparent_magnitudes()
    
    extinction = { 'AV': 1.0, 'Ag': 1.0, 'Ar': 1.0, 'Ai': 1.0, 
                  'EBV': 0.1,
                  'Rg': 1.0, 'Rr': 1.0, 'Ri': 1.0 }
    
    binary.apply_extinction(extinction)
    
    for f in ['V','g','r','i']:
        
        dmag = getattr(binary,'m'+f+'_combined_corr') - getattr(binary,'m'+f+'_combined')
        
        assert dmag == extinction['A'+f]
    
    dcol = round(binary.BV_combined_corr - binary.BV_combined,1)
    
    assert dcol == extinction['EBV']
    
    c1 = ['g', 'r']
    c2 = ['r', 'i']
        
    for i in range(0,2,1):
        
        mcorr1 = getattr(binary,'m'+c1[i]+'_combined') - extinction['A'+c1[i]]
        mcorr2 = getattr(binary,'m'+c2[i]+'_combined') - extinction['A'+c2[i]]
        
        m1 = getattr(binary,'m'+c1[i]+'_combined_corr')
        m2 = getattr(binary,'m'+c2[i]+'_combined_corr')
        
        assert (m1-m2) == (mcorr1-mcorr2)

    
if __name__ == '__main__':
    
    test_binary_calculate_apparent_magnitudes()