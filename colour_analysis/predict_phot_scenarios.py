# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 15:27:52 2018

@author: rstreet
"""

from os import path
from sys import argv
import combined_light_predictor
import jester_phot_transforms

def predict_photometry(output_path,dm,extinction):
    """Driver code to set up scenarios and generate output"""
    
    #kstar = Kdwarf(name='K star', dm=dm)
    #kstar.calculate_apparent_magnitudes()
    #kstar.apply_extinction(extinction)
    #kstar.transform_extinc_corr_Johnson_SDSS()
    
    m3star = M3dwarf(name='M star', dm=dm)
    m3star.calculate_apparent_magnitudes()
    m3star.apply_extinction(extinction)
    #mstar.transform_extinc_corr_Johnson_SDSS()
    
    m7star = M7dwarf(name='M star', dm=dm)
    m7star.calculate_apparent_magnitudes()
    m7star.apply_extinction(extinction)
    
    #wd_hot = WhiteDwarf(100000.0,dm=dm)
    #wd_hot.calculate_apparent_magnitudes()
    #wd_hot.apply_extinction(extinction)
    #wd_hot.transform_extinc_corr_Johnson_SDSS()
    
    #wd_mid = WhiteDwarf(20000.0,dm=dm)
    #wd_mid.calculate_apparent_magnitudes()
    #wd_mid.apply_extinction(extinction)
    #wd_mid.transform_extinc_corr_Johnson_SDSS()
    
    #wd_cool = WhiteDwarf(4000.0,dm=dm)
    #wd_cool.calculate_apparent_magnitudes()
    #wd_cool.apply_extinction(extinction)
    #wd_cool.transform_extinc_corr_Johnson_SDSS()
    
    ms_binary = scenario1(dm,extinction)
    
    #k_wd_binary = scenario2(dm,extinction)
    
    #m_wd_binary = scenario3(dm,extinction)
    
    output_latex_table(m3star,m7star, ms_binary, output_path)
                       

def Kdwarf(name=None,dm=None):
    """Definition of a single, main sequence K-dwarf star of mass 0.88 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 5.749
    star.Mr = 5.199
    star.Mi = 5.058
    star.MB = 6.368
    star.MV = 5.584
    star.MJ = 4.329
    star.MH = 3.942
    star.MKs = 3.897
    
    star.distance_modulus = dm
    
    star.calculate_colours()
        
    return star

def M4dwarf(name=None,dm=None):
    """Definition of a single, main sequence M-dwarf star of mass 0.25 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 13.650
    star.Mr = 12.081
    star.Mi = 10.714
    star.MB = 14.350
    star.MV = 12.676
    star.MJ = 8.398
    star.MH = 7.910
    star.MKs = 7.676
    
    star.distance_modulus = dm
    
    star.calculate_colours()
    
    return star

def M3dwarf(name=None,dm=None):
    """Definition of a single, main sequence M-dwarf star of mass 0.35 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 11.933
    star.Mr = 10.409
    star.Mi = 9.475
    star.MB = 13.175
    star.MV = 11.574
    star.MJ = 7.566
    star.MH = 7.014
    star.MKs = 6.779
    
    star.distance_modulus = dm
    
    star.calculate_colours()
    
    return star


def M7dwarf(name=None,dm=None):
    """Definition of a single, main sequence M-dwarf star of mass 0.09 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 19.149
    star.Mr = 17.181
    star.Mi = 14.700
    star.MB = 21.124
    star.MV = 18.674
    star.MJ = 11.033
    star.MH = 10.458
    star.MKs = 10.178
    
    star.distance_modulus = dm
    
    star.calculate_colours()
    
    return star

def WhiteDwarf(teff, dm=None):
    """Definition of a single White Dwarf star, with parameters determined
    by the t_eff provided
    teff one of { 4000.0, 20000.0, 100000.0 }
    Photometry from Bergeron 1995 PASP 107 1047.
    """

    star = combined_light_predictor.Star()
    
    star.distance_modulus = dm
    
    if teff == 4000.0:
        
        star.name = 'WD (4000K)'
        star.MV = 16.105
        star.MB = 17.081
        star.BV = 0.976
        star.VI = 1.249
        star.RI = 0.620
        star.VKs = 1.309
        star.JH = -0.152
        star.HKs = -0.126
    
    elif teff == 20000.0:
        
        star.name = 'WD (20000K)'
        star.MV = 10.732
        star.MB = 10.688
        star.BV = -0.044
        star.VI = -0.216
        star.RI = -0.119
        star.VKs = -0.705
        star.JH = -0.063
        star.HKs = -0.106
    
    elif teff == 100000.0:
        
        star.name = 'WD (100000K)'
        star.MV = 8.408
        star.MB = 8.079
        star.BV = -0.329
        star.VI = -0.358
        star.RI = -0.210
        star.VKs = -1.093
        star.JH = -0.131
        star.HKs = -0.134
    
    else:
        
        print('ERROR: No photometry available for a white dwarf of teff='+str(teff)+'K')
        
        return star
        
    star.MI = star.MV - star.VI
    star.MKs = star.MV - star.VKs
    star.MH = star.MKs + star.HKs
    star.MJ = star.MH + star.JH
    star.JKs = star.MJ - star.MKs
    
    phot = jester_phot_transforms.transform_JohnsonCousins_to_SDSS(BV=star.BV, sigBV=0.0, 
                                  RI=star.RI, sigRI=0.0,
                                  V=star.MV, sigV=0.0)
    
    star.gr = phot['g-r']
    star.ri = phot['r-i']
    star.Mg = phot['g']
    star.Mr = phot['r']
    star.Mi = star.Mr - star.ri
    
    return star

    
def scenario1(dm,extinction):
    """Define a main sequence binary consisting of a K-dwarf and an M-dwarf"""

    m3star = M3dwarf(name='K star', dm=dm)
    m7star = M7dwarf(name='M star', dm=dm)
    
    binary = combined_light_predictor.BinaryStar()
    binary.distance_modulus = dm
    
    binary.star1 = m3star
    binary.star2 = m7star
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_combined_colours()
    
    binary.apply_extinction(extinction)
    
    return binary

def scenario2(dm,extinction):
    """Define a binary consisting of a mid-range temperature white dwarf 
    plus a main sequence K-star"""
    
    kstar = Kdwarf(name='K star', dm=dm)
    wd = WhiteDwarf(20000.0,dm=dm)
    
    binary = combined_light_predictor.BinaryStar()
    binary.distance_modulus = dm
    
    binary.star1 = kstar
    binary.star2 = wd
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_combined_colours()

    binary.apply_extinction(extinction)
    
    return binary
   

def scenario3(dm,extinction):
    """Define a binary consisting of a mid-range temperature white dwarf 
    plus a main sequence M-star"""
    
    mstar = Mdwarf(name='M star', dm=dm)
    wd = WhiteDwarf(20000.0,dm=dm)
    
    binary = combined_light_predictor.BinaryStar()
    binary.distance_modulus = dm
    
    binary.star1 = wd
    binary.star2 = mstar
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_combined_colours()

    binary.apply_extinction(extinction)
    
    return binary

def roundfloat(value,ndp):
    """Function to convert a floating point value to a string, rounding
    the the given number of decimal places in the process. 
    The resulting string will be end zero-padded to ensure a consistent
    number of characters are returned.
    """
    
    output = str( round(value,ndp) )
    
    while len(output.split('.')[1]) < ndp:
        output = output + '0'
    
    return output
    
def output_latex_table(m3star,m7star, ms_binary, output_path):
    """Function to output a LaTeX-format table of the stellar photometry"""
    
    t = open(output_path, 'w')
    
    t.write('\\begin{table}[h!]\n')
    t.write('\\centering\n')
    t.write('\\caption{Predicted photometric properties of the lens system.  Apparent magnitudes are calculated for the measured lens distance without extinction or reddening, except for the bottom section.  } \\label{tab:binaryphot}\n')
    t.write('\\begin{tabular}{lccc}\n')
    t.write('\\hline\n')
    t.write('\\hline\n')
    t.write('Quantity	   & M3-dwarf      & M7-dwarf   & MS-binary \\\\\n')
    t.write('\\[mag\\]	   &		       & 	        & M3+M7     \\\\\n')
    t.write('$M_{B}$       & '+roundfloat(m3star.MB,3)+' & '+roundfloat(m7star.MB,3)+' & '+roundfloat(ms_binary.MB_combined,3)+'\\\\\n')
    t.write('$M_{V}$       & '+roundfloat(m3star.MV,3)+' & '+roundfloat(m7star.MV,3)+' & '+roundfloat(ms_binary.MV_combined,3)+'\\\\\n')
    t.write('$M_{g}$       & '+roundfloat(m3star.Mg,3)+' & '+roundfloat(m7star.Mg,3)+' & '+roundfloat(ms_binary.Mg_combined,3)+'\\\\\n')
    t.write('$M_{r}$       & '+roundfloat(m3star.Mr,3)+' & '+roundfloat(m7star.Mr,3)+' & '+roundfloat(ms_binary.Mr_combined,3)+'\\\\\n')
    t.write('$M_{i}$       & '+roundfloat(m3star.Mi,3)+' & '+roundfloat(m7star.Mi,3)+' & '+roundfloat(ms_binary.Mi_combined,3)+'\\\\\n')
    t.write('$M_{J}$       & '+roundfloat(m3star.MJ,3)+' & '+roundfloat(m7star.MJ,3)+' & '+roundfloat(ms_binary.MJ_combined,3)+'\\\\\n')
    t.write('$M_{H}$       & '+roundfloat(m3star.MH,3)+' & '+roundfloat(m7star.MH,3)+' & '+roundfloat(ms_binary.MH_combined,3)+'\\\\\n')
    t.write('$M_{Ks}$      & '+roundfloat(m3star.MKs,3)+' & '+roundfloat(m7star.MKs,3)+' & '+roundfloat(ms_binary.MKs_combined,3)+'\\\\\n')
    t.write('$(B-V)$       & '+roundfloat(m3star.BV,3)+' & '+roundfloat(m7star.BV,3)+' & '+roundfloat(ms_binary.BV_combined,3)+'\\\\\n')
    t.write('$(g-r)$       & '+roundfloat(m3star.gr,3)+' & '+roundfloat(m7star.gr,3)+' & '+roundfloat(ms_binary.gr_combined,3)+'\\\\\n')
    t.write('$(r-i)$       & '+roundfloat(m3star.ri,3)+' & '+roundfloat(m7star.ri,3)+' & '+roundfloat(ms_binary.ri_combined,3)+'\\\\\n')
    t.write('$(J-H)$       & '+roundfloat(m3star.JH,3)+' & '+roundfloat(m7star.JH,3)+' & '+roundfloat(ms_binary.JH_combined,3)+'\\\\\n')
    t.write('$(H-Ks)$      & '+roundfloat(m3star.HKs,3)+' & '+roundfloat(m7star.HKs,3)+' & '+roundfloat(ms_binary.HKs_combined,3)+'\\\\\n')
    t.write('$(J-Ks)$      & '+roundfloat(m3star.JKs,3)+' & '+roundfloat(m7star.JKs,3)+' & '+roundfloat(ms_binary.JKs_combined,3)+'\\\\\n')
    t.write('$m_{B}$       & '+roundfloat(m3star.mB,3)+' & '+roundfloat(m7star.mB,3)+' & '+roundfloat(ms_binary.mB_combined,3)+'\\\\\n')
    t.write('$m_{V}$       & '+roundfloat(m3star.mV,3)+' & '+roundfloat(m7star.mV,3)+' & '+roundfloat(ms_binary.mV_combined,3)+'\\\\\n')
    t.write('$m_{g}$       & '+roundfloat(m3star.mg,3)+' & '+roundfloat(m7star.mg,3)+' & '+roundfloat(ms_binary.mg_combined,3)+'\\\\\n')
    t.write('$m_{r}$       & '+roundfloat(m3star.mr,3)+' & '+roundfloat(m7star.mr,3)+' & '+roundfloat(ms_binary.mr_combined,3)+'\\\\\n')
    t.write('$m_{i}$       & '+roundfloat(m3star.mi,3)+' & '+roundfloat(m7star.mi,3)+' & '+roundfloat(ms_binary.mi_combined,3)+'\\\\\n')
    t.write('$m_{J}$       & '+roundfloat(m3star.mJ,3)+' & '+roundfloat(m7star.mJ,3)+' & '+roundfloat(ms_binary.mJ_combined,3)+'\\\\\n')
    t.write('$m_{H}$       & '+roundfloat(m3star.mH,3)+' & '+roundfloat(m7star.mH,3)+' & '+roundfloat(ms_binary.mH_combined,3)+'\\\\\n')
    t.write('$m_{Ks}$      & '+roundfloat(m3star.mKs,3)+' & '+roundfloat(m7star.mKs,3)+' & '+roundfloat(ms_binary.mKs_combined,3)+'\\\\\n')
    t.write('\hline\n')
    t.write('$m_{V,corr}$  & '+roundfloat(m3star.mV_corr,3)+' & '+roundfloat(m7star.mV_corr,3)+' & '+roundfloat(ms_binary.mV_combined_corr,3)+'\\\\\n')
    t.write('$(B-V)_{corr}$& '+roundfloat(m3star.BV_corr,3)+' & '+roundfloat(m7star.BV_corr,3)+' & '+roundfloat(ms_binary.BV_combined_corr,3)+'\\\\\n')
    t.write('$m_{g,corr}$  & '+roundfloat(m3star.mg_corr,3)+' & '+roundfloat(m7star.mg_corr,3)+' & '+roundfloat(ms_binary.mg_combined_corr,3)+'\\\\\n')
    t.write('$m_{r,corr}$  & '+roundfloat(m3star.mr_corr,3)+' & '+roundfloat(m7star.mr_corr,3)+' & '+roundfloat(ms_binary.mr_combined_corr,3)+'\\\\\n')
    t.write('$m_{i,corr}$  & '+roundfloat(m3star.mi_corr,3)+' & '+roundfloat(m7star.mi_corr,3)+' & '+roundfloat(ms_binary.mi_combined_corr,3)+'\\\\\n')
    t.write('$(g-r)_{corr}$& '+roundfloat(m3star.gr_corr,3)+' & '+roundfloat(m7star.gr_corr,3)+' & '+roundfloat(ms_binary.gr_combined_corr,3)+'\\\\\n')
    t.write('$(r-i)_{corr}$& '+roundfloat(m3star.ri_corr,3)+' & '+roundfloat(m7star.ri_corr,3)+' & '+roundfloat(ms_binary.ri_combined_corr,3)+'\\\\\n')
    t.write('\hline\n')
    t.write('\end{tabular}\n')
    t.write('\end{table}\n')
    
    t.close()

def option_or_none(option):
    
    try:
        value = float(option)
        
    except ValueError:
        value = None
    
    return value
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        
        output_path = raw_input('Please enter the path to the output file: ')
        dm = float(raw_input('Please enter the distance modulus: '))
        AV = option_or_none(raw_input('Please enter the extinction Av [or None]: '))
        EBV = option_or_none(raw_input('Please enter the reddening, E(B-V) [or None]: '))
        Ag = option_or_none(raw_input('Please enter the extinction A(SDSS-g) [or None]: '))
        Ar = option_or_none(raw_input('Please enter the extinction A(SDSS-r) [or None]: '))
        Ai = option_or_none(raw_input('Please enter the extinction A(SDSS-i) [or None]: '))
        
    else:
        output_path = argv[1]
        dm = float(argv[2])
        AV = option_or_none(argv[3])
        EBV = option_or_none(argv[4])
        Ag = option_or_none(argv[5])
        Ar = option_or_none(argv[6])
        Ai = option_or_none(argv[7])
    
    extinction = { 'AV': AV, 'Ag': Ag, 'Ar': Ar, 'Ai': Ai, 'EBV': EBV }
    
    predict_photometry(output_path,dm,extinction)