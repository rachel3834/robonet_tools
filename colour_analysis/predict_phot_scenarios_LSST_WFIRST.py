# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 15:27:52 2018

@author: rstreet
"""

from os import path
from sys import argv
import combined_light_predictor
import jester_phot_transforms

def predict_photometry(dm):
    """Driver code to set up scenarios and generate output"""
    
    gstar = Gdwarf(name='G star', dm=dm)
    gstar.estimate_wfirst_abs_magnitudes()
    gstar.calculate_apparent_magnitudes()
    print('G star at '+str(dm)+': '+gstar.summary_app())
    
    kstar = Kdwarf(name='K star', dm=dm)
    kstar.estimate_wfirst_abs_magnitudes()
    kstar.calculate_apparent_magnitudes()
    print('K star at '+str(dm)+': '+kstar.summary_app())
    
    mstar = Mdwarf(name='M star', dm=dm)
    mstar.estimate_wfirst_abs_magnitudes()
    mstar.calculate_apparent_magnitudes()
    print('M star at '+str(dm)+': '+mstar.summary_app())
    
    ltdwarf = LTdwarf(name='Brown dwarf', dm=dm)
    ltdwarf.estimate_wfirst_abs_magnitudes()
    ltdwarf.calculate_apparent_magnitudes()
    print('LT BD at '+str(dm)+': '+ltdwarf.summary_app())
    
def Gdwarf(name=None,dm=None):
    """Definition of a single, main sequence G-dwarf star of mass 1.0 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 5.011
    star.Mr = 4.574
    star.Mi = 4.475
    star.Mz = 4.465
    star.MJ = 3.908
    star.MH = 3.586
    star.MKs = 3.555
    
    star.distance_modulus = dm
    
    star.calculate_colours()
        
    return star
                   

def Kdwarf(name=None,dm=None):
    """Definition of a single, main sequence K-dwarf star of mass 0.6 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 9.479
    star.Mr = 8.189
    star.Mi = 7.539
    star.Mz = 7.150
    star.MJ = 5.936
    star.MH = 5.333
    star.MKs = 5.196
    
    star.distance_modulus = dm
    
    star.calculate_colours()
        
    return star

def Mdwarf(name=None,dm=None):
    """Definition of a single, main sequence M-dwarf star of mass 0.25 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 13.650
    star.Mr = 12.081
    star.Mi = 10.714
    star.Mz = 9.927
    star.MJ = 8.398
    star.MH = 7.910
    star.MKs = 7.676
    
    star.distance_modulus = dm
    
    star.calculate_colours()
    
    return star

def LTdwarf(name=None,dm=None):
    """Definition of a single L/T brown dwarf star of type T0 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 99.999
    star.Mr = 99.999
    star.Mi = 22.94
    star.Mz = 19.11
    star.MJ = 15.89
    star.MH = 15.09
    star.MKs = 14.46
    
    star.distance_modulus = dm
    
    star.calculate_colours()
    
    return star

    
def scenario1(dm):
    """Define a main sequence binary consisting of a K-dwarf and an M-dwarf"""

    kstar = Kdwarf(name='K star', dm=dm)
    mstar = Mdwarf(name='M star', dm=dm)
    
    binary = combined_light_predictor.BinaryStar()
    binary.distance_modulus = dm
    
    binary.star1 = kstar
    binary.star2 = mstar
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_combined_colours()

    return binary

def scenario2(dm):
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

    return binary
   

def scenario3(dm=None):
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

    return binary

def output_latex_table(kstar,mstar,wd_hot,wd_mid,wd_cool, 
                       ms_binary, k_wd_binary, m_wd_binary, output_path):
    """Function to output a LaTeX-format table of the stellar photometry"""
    
    t = open(output_path, 'w')
    
    t.write('\\begin{table}[h!]\n')
    t.write('\\centering\n')
    t.write('\\caption{Predicted photometric properties of the lens system, if both stars are on the main sequence and have solar age and metallicity.  Apparent magnitudes are calculated for the measured lens distance without extinction or reddening.} \\label{tab:binaryphot}\n')
    t.write('\\begin{tabular}{lccccccc}\n')
    t.write('\\hline\n')
    t.write('\\hline\n')
    t.write('Quantity	   & K-dwarf & M-dwarf & White Dwarf         & White Dwarf          & White Dwarf            &Main sequence & White dwarf +   & White dwarf +  \\\n')
    t.write('		   &		 & 	     & $T_{eff}$=4000\\,K  & $T_{eff}$=20,000\\,K & $T_{eff}$=100,000\\,K  & K+M binary   & K-dwarf binary  & M-dwarf binary \\\n')
    t.write('$M_{g}$ [mag] & '+str(round(kstar.Mg,3))+' & '+str(round(mstar.Mg,3))+' & '+str(round(wd_cool.Mg,3))+' & '+str(round(wd_mid.Mg,3))+' & '+str(round(wd_hot.Mg,3))+' & '+str(round(ms_binary.Mg_combined,3))+' & '+str(round(k_wd_binary.Mg_combined,3))+' & '+str(round(m_wd_binary.Mg_combined,3))+'\n')
    t.write('$M_{r}$ [mag] & '+str(round(kstar.Mr,3))+' & '+str(round(mstar.Mr,3))+' & '+str(round(wd_cool.Mr,3))+' & '+str(round(wd_mid.Mr,3))+' & '+str(round(wd_hot.Mr,3))+' & '+str(round(ms_binary.Mr_combined,3))+' & '+str(round(k_wd_binary.Mr_combined,3))+' & '+str(round(m_wd_binary.Mr_combined,3))+'\n')
    t.write('$M_{i}$ [mag] & '+str(round(kstar.Mi,3))+' & '+str(round(mstar.Mi,3))+' & '+str(round(wd_cool.Mi,3))+' & '+str(round(wd_mid.Mi,3))+' & '+str(round(wd_hot.Mi,3))+' & '+str(round(ms_binary.Mi_combined,3))+' & '+str(round(k_wd_binary.Mi_combined,3))+' & '+str(round(m_wd_binary.Mi_combined,3))+'\n')
    t.write('$M_{J}$ [mag] & '+str(round(kstar.MJ,3))+' & '+str(round(mstar.MJ,3))+' & '+str(round(wd_cool.MJ,3))+' & '+str(round(wd_mid.MJ,3))+' & '+str(round(wd_hot.MJ,3))+' & '+str(round(ms_binary.MJ_combined,3))+' & '+str(round(k_wd_binary.MJ_combined,3))+' & '+str(round(m_wd_binary.MJ_combined,3))+'\n')
    t.write('$M_{H}$ [mag] & '+str(round(kstar.MH,3))+' & '+str(round(mstar.MH,3))+' & '+str(round(wd_cool.MH,3))+' & '+str(round(wd_mid.MH,3))+' & '+str(round(wd_hot.MH,3))+' & '+str(round(ms_binary.MH_combined,3))+' & '+str(round(k_wd_binary.MH_combined,3))+' & '+str(round(m_wd_binary.MH_combined,3))+'\n')
    t.write('$M_{Ks}$ [mag]& '+str(round(kstar.MKs,3))+' & '+str(round(mstar.MKs,3))+' & '+str(round(wd_cool.MKs,3))+' & '+str(round(wd_mid.MKs,3))+' & '+str(round(wd_hot.MKs,3))+' & '+str(round(ms_binary.MKs_combined,3))+' & '+str(round(k_wd_binary.MKs_combined,3))+' & '+str(round(m_wd_binary.MKs_combined,3))+'\n')
    t.write('$(g-r)$ [mag] & '+str(round(kstar.gr,3))+' & '+str(round(mstar.gr,3))+' & '+str(round(wd_cool.gr,3))+' & '+str(round(wd_mid.gr,3))+' & '+str(round(wd_hot.gr,3))+' & '+str(round(ms_binary.gr_combined,3))+' & '+str(round(k_wd_binary.gr_combined,3))+' & '+str(round(m_wd_binary.gr_combined,3))+'\n')
    t.write('$(r-i)$ [mag] & '+str(round(kstar.ri,3))+' & '+str(round(mstar.ri,3))+' & '+str(round(wd_cool.ri,3))+' & '+str(round(wd_mid.ri,3))+' & '+str(round(wd_hot.ri,3))+' & '+str(round(ms_binary.ri_combined,3))+' & '+str(round(k_wd_binary.ri_combined,3))+' & '+str(round(m_wd_binary.ri_combined,3))+'\n')
    t.write('$(J-H)$ [mag] & '+str(round(kstar.JH,3))+' & '+str(round(mstar.JH,3))+' & '+str(round(wd_cool.JH,3))+' & '+str(round(wd_mid.JH,3))+' & '+str(round(wd_hot.JH,3))+' & '+str(round(ms_binary.JH_combined,3))+' & '+str(round(k_wd_binary.JH_combined,3))+' & '+str(round(m_wd_binary.JH_combined,3))+'\n')
    t.write('$(H-Ks)$ [mag] & '+str(round(kstar.HKs,3))+' & '+str(round(mstar.HKs,3))+' & '+str(round(wd_cool.HKs,3))+' & '+str(round(wd_mid.HKs,3))+' & '+str(round(wd_hot.HKs,3))+' & '+str(round(ms_binary.HKs_combined,3))+' & '+str(round(k_wd_binary.HKs_combined,3))+' & '+str(round(m_wd_binary.HKs_combined,3))+'\n')
    t.write('$(J-Ks)$ [mag] & '+str(round(kstar.JKs,3))+' & '+str(round(mstar.JKs,3))+' & '+str(round(wd_cool.JKs,3))+' & '+str(round(wd_mid.JKs,3))+' & '+str(round(wd_hot.JKs,3))+' & '+str(round(ms_binary.JKs_combined,3))+' & '+str(round(k_wd_binary.JKs_combined,3))+' & '+str(round(m_wd_binary.JKs_combined,3))+'\n')
    t.write('$m_{g}$ [mag] & '+str(round(kstar.mg,3))+' & '+str(round(mstar.mg,3))+' & '+str(round(wd_cool.mg,3))+' & '+str(round(wd_mid.mg,3))+' & '+str(round(wd_hot.mg,3))+' & '+str(round(ms_binary.mg_combined,3))+' & '+str(round(k_wd_binary.mg_combined,3))+' & '+str(round(m_wd_binary.mg_combined,3))+'\n')
    t.write('$m_{r}$ [mag] & '+str(round(kstar.mr,3))+' & '+str(round(mstar.mr,3))+' & '+str(round(wd_cool.mr,3))+' & '+str(round(wd_mid.mr,3))+' & '+str(round(wd_hot.mr,3))+' & '+str(round(ms_binary.mr_combined,3))+' & '+str(round(k_wd_binary.mr_combined,3))+' & '+str(round(m_wd_binary.mr_combined,3))+'\n')
    t.write('$m_{i}$ [mag] & '+str(round(kstar.mi,3))+' & '+str(round(mstar.mi,3))+' & '+str(round(wd_cool.mi,3))+' & '+str(round(wd_mid.mi,3))+' & '+str(round(wd_hot.mi,3))+' & '+str(round(ms_binary.mi_combined,3))+' & '+str(round(k_wd_binary.mi_combined,3))+' & '+str(round(m_wd_binary.mi_combined,3))+'\n')
    t.write('$m_{J}$ [mag] & '+str(round(kstar.mJ,3))+' & '+str(round(mstar.mJ,3))+' & '+str(round(wd_cool.mJ,3))+' & '+str(round(wd_mid.mJ,3))+' & '+str(round(wd_hot.mJ,3))+' & '+str(round(ms_binary.mJ_combined,3))+' & '+str(round(k_wd_binary.mJ_combined,3))+' & '+str(round(m_wd_binary.mJ_combined,3))+'\n')
    t.write('$m_{H}$ [mag] & '+str(round(kstar.mH,3))+' & '+str(round(mstar.mH,3))+' & '+str(round(wd_cool.mH,3))+' & '+str(round(wd_mid.mH,3))+' & '+str(round(wd_hot.mH,3))+' & '+str(round(ms_binary.mH_combined,3))+' & '+str(round(k_wd_binary.mH_combined,3))+' & '+str(round(m_wd_binary.mH_combined,3))+'\n')
    t.write('$m_{Ks}$ [mag] & '+str(round(kstar.mKs,3))+' & '+str(round(mstar.mKs,3))+' & '+str(round(wd_cool.mKs,3))+' & '+str(round(wd_mid.mKs,3))+' & '+str(round(wd_hot.mKs,3))+' & '+str(round(ms_binary.mKs_combined,3))+' & '+str(round(k_wd_binary.mKs_combined,3))+' & '+str(round(m_wd_binary.mKs_combined,3))+'\n')
    t.write('\hline\n')
    t.write('\end{tabular}\n')
    t.write('\end{table}\n')
    
    t.close()
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        
        dm = float(raw_input('Please enter the distance modulus: '))
        
    else:
        dm = float(argv[1])
    
    predict_photometry(dm)