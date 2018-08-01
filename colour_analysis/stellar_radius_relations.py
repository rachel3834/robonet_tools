# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 13:56:20 2018

@author: rstreet
"""
import numpy as np

def calc_star_ang_radius_Adams2018(Q,sigQ,PQ,sigPQ,pqcolour,Lclass=None):
    """Function to calculate the limb-darkened angular diameter of a star, 
    given one of the photometric colours:
    (I-H), (I-K), (V-I), (V-H), (V-K), 
    
    using the expression:
    
    log(theta_Q=0) = log(theta_LD) + 0.2Q, 
    
    log(theta_Q=0) = Sum_n=0->N [ c_n (P-Q)**n ]
    where
    
    Q = apparent magnitude in a given band
    (P-Q) = the measured colour
    theta_LD = angular diameter, corrected for limb-darkening
    
    from Adams et al. (2018), MNRAS, 473, 3608. 
    Returned units are log10(micro-arcsec).
    """
    
    coeffs = fetch_coefficients_Adams2018(pqcolour,Lclass=Lclass)
    
    log_theta_Q0 = 0.0
    var_log_theta_Q0 = sigPQ * sigPQ
    
    for n in range(0,3,1):
        log_theta_Q0 += coeffs['c'][n] * PQ**n
        var_log_theta_Q0 += coeffs['sig_c'][n]*coeffs['sig_c'][n]
    
    log_theta_LD = log_theta_Q0 - 0.2 * Q
    
    sig_log_theta_LD = np.sqrt( var_log_theta_Q0 + ( sigQ*sigQ ) )
    
    return log_theta_LD, sig_log_theta_LD


def fetch_coefficients_Adams2018(pqcolour,Lclass=None):
    """Function to return the correct set of co-efficients for the
    stellar angular diameter calculation.
    pqcolour        str    one of {I-H, I-K, V-I, V-H, V-K}
    spectral_type   str    one of {None, dwarfs, subgiants, giants}
    """
    
    if Lclass == None:
        coeffs = {'I-H': {'c': [0.541, 0.133, 0.0],
                          'sig_c': [0.004,0.003,0.0],
                          'rms': 0.025,
                          'valid_range': [-0.042,1.935]},
                  'I-K': {'c': [0.528, 0.108, 0.0],
                          'sig_c': [0.005,0.003,0.0],
                          'rms': 0.020,
                          'valid_range': [0.019,2.159]},
                  'V-I': {'c': [0.542, 0.391, 0.0],
                          'sig_c': [0.006, 0.006, 0.0],
                          'rms': 0.028,
                          'valid_range': [-0.050,1.740]},
                  'V-H': {'c': [0.538, 0.074, 0.0],
                          'sig_c': [0.004, 0.002, 0.0],
                          'rms': 0.020,
                          'valid_range': [-0.052,3.615]},
                  'V-K': {'c': [0.529, 0.062, 0.0],
                          'sig_c': [0.004, 0.002, 0.0],
                          'rms': 0.021,
                          'valid_range': [-0.021,3.839]},
                  }
    
    elif Lclass == 'dwarfs' or Lclass == 'subgiants':
        coeffs = {'I-H': {'c': [0.529, 0.166, 0.0],
                          'sig_c': [0.007, 0.010, 0.0],
                          'rms': 0.031,
                          'valid_range': [-0.042,1.198]},
                  'I-K': {'c': [0.520, 0.118, 0.0], 
                          'sig_c': [0.007, 0.010, 0.0],
                          'rms': 0.023,
                          'valid_range': [0.019,1.309]},
                  'V-I': {'c': [0.542, 0.378, 0.0], 
                          'sig_c': [0.007, 0.011, 0.0],
                          'rms': 0.029,
                          'valid_range': [-0.050,1.160]},
                  'V-H': {'c': [0.534, 0.079, 0.0,],
                          'sig_c': [0.005, 0.003, 0.0],
                          'rms': 0.023,
                          'valid_range': [-0.052,3.447]},
                  'V-K': {'c': [0.523, 0.063, 0.0], 
                          'sig_c': [0.006, 0.005, 0.0],
                          'rms': 0.024,
                          'valid_range': [-0.021,2.209]}
                  }
                  
    elif Lclass == 'giants':
        coeffs = {'I-H': {'c': [0.523, 0.144, 0.0], 
                          'sig_c': [0.011, 0.007, 0.0],
                          'rms': 0.020,
                          'valid_range': [0.065,1.935]},
                  'I-K': {'c': [0.543, 0.098, 0.0],
                          'sig_c': [0.011, 0.007, 0.0],
                          'rms': 0.016,
                          'valid_range': [1.129,2.159]},
                  'V-I': {'c': [0.535, 0.490, -0.068],
                          'sig_c': [0.027, 0.046, 0.019],
                          'rms': 0.026,
                          'valid_range': [-0.010,1.740]},
                  'V-H': {'c': [0.532, 0.076, 0.0],
                          'sig_c': [0.009, 0.003, 0.0],
                          'rms': 0.016,
                          'valid_range': [0.055,3.615]},
                  'V-K': {'c': [0.562, 0.051, 0.0],
                          'sig_c': [0.009, 0.003, 0.0],
                          'rms': 0.019,
                          'valid_range': [2.049,3.839]}
                  }
    else:
        print('ERROR: unrecognised spectral type class '+Lclass)
        
    return coeffs[pqcolour]

def calc_star_ang_radius_Boyajian2014(colour,sig_colour,mag,sig_mag,pqcolour,FeH):
    """Function to calculate the angular diameter of a star from its (g-r) colour
    and metallicity.  
        using the expression:
    
    log(theta_mV=0) = log(theta_LD) + 0.2mV, 
    
    where:
    
    log(theta_m=0) = Sum_n=0->N [ a_n (X)**n ]
    
    where
    
    mV = star apparent magnitude in V
    (X) = star colour in SDSS passbands, one of {g-r, g-i}
    [Fe/H] = star metallicity
    
    and the RMS uncertainty on log(theta_mV) = 5.8%
    
    taken from Boyajian et al (2014) AJ, 147, 47.
    """
    
    coeffs = fetch_coefficients_Boyajian2014(pqcolour)
    
    log_theta_m0 = 0.0
    var_log_theta_m0 = sig_colour * sig_colour
    
    for n in range(0,5,1):
        log_theta_m0 += coeffs['a'][n] * colour**n
        var_log_theta_m0 += coeffs['sig_a'][n]*coeffs['sig_a'][n]
    

    log_theta_LD = log_theta_m0 - 0.2*mag
    
    sig_log_theta_LD = np.sqrt( var_log_theta_m0 + (sig_mag*sig_mag) )
    
    return log_theta_LD, sig_log_theta_LD

def fetch_coefficients_Boyajian2014(pqcolour):
    """Function to return the correct set of co-efficients for the
    stellar angular diameter calculation.
    pqcolour        str    one of {g-r, g-i}
    """
    
    coeffs = {'g-r': {'a': [0.66728, 0.58135, 0.88293, -1.41005, 0.67248],
                      'sig_a': [0.00203, 0.01180, 0.03470, 0.04331, 0.01736],
                      'rms': 0.097,
                      'valid_range': [-0.23,-1.40]},
              'g-i': {'a': [0.69174, 0.54346, -0.02149, 0.0, 0.0], 
                      'sig_a': [0.00125, 0.00266, 0.00097, 0.0, 0.0],
                      'rms': 0.092,
                      'valid_range': [-0.43,-2.78]}
             }
             
    return coeffs[pqcolour]
  
def scale_source_distance(Mv, theta_deg, err_theta_deg, DS):
    """Function to calculate the apparent angular size of the source at the 
    inferred source distance.  
    """
    
    Rsol = 6.96e10   # cm
    pc = 3.086e18   # cm
    
    # Calculating the distance at which a star of the given absolute V mag given
    # would have an apparent V mag of zero:
    DMv0 = 10**( ( -1.0 * Mv + 5.0 ) / 5.0 )
    print 'Distance mV=0: '+str(DMv0)+'pc.  Zero uncertainty as model value for Mv.'
    
    # Based on the distance, calculating the physical radius of the star:
    Rstar = DMv0 * math.sin(theta_deg*math.pi/180.0 / 2.0)
    Rstar_min = DMv0 * math.sin((theta_deg-err_theta_deg)*math.pi/180.0 / 2.0)
    Rstar_max = DMv0 * math.sin((theta_deg+err_theta_deg)*math.pi/180.0 / 2.0)
    sig_Rstar = ( ( Rstar_max - Rstar_min ) / 2.0 ) * pc / Rsol
    print 'Stellar radius '+str(Rstar * pc / Rsol)+'+/-'+str(sig_Rstar)+'Rsolar'
    
    # Calculating the apparent angular radius this star would have if
    # it were at the distance of the source as given:
    # Note no error is included for the distance since this value is expected
    # to be an estimate. 
    def calc_theta_S(Rstar, DS):
        sin_theta_S = Rstar / DS
        theta_S = math.asin( sin_theta_S ) * 180.0 / math.pi * 3600.0 * 1e6
        return sin_theta_S, theta_S
    (sin_theta_S, theta_S) = calc_theta_S(Rstar, DS)
    (sin_theta_S_min, theta_S_min) = calc_theta_S(Rstar_min, DS)
    (sin_theta_S_max, theta_S_max) = calc_theta_S(Rstar_max, DS)
    sig_theta_S = ( theta_S_max - theta_S_min ) / 2.0
    
    print 'Apparent angular source radius at D_S:'+str(theta_S)+'+/-'+str(sig_theta_S)+'mu-as'


def star_mass_radius_relation(teff, logg, FeH, log=None):
    """The relationship between mass and radius for dwarf and giant stars, 
    as derived by Torres et al. (2010), A&ARv, 18, 67."""
    
    b1 = 2.4427
    sig_b1 = 0.038
    b2 = 0.6679
    sig_b2 = 0.016
    b3 = 0.1771
    sig_b3 = 0.027
    b4 = 0.705 
    sig_b4 = 0.13
    b5 = -0.21415
    sig_b5 = 0.0075
    b6 = 0.02306
    sig_b6 = 0.0013
    b7 = 0.04173
    sig_b7 = 0.0082
    
    X = np.log10(teff) - 4.1
    
    log_R = b1 + b2*X + b3*X*X + b4*X*X*X + \
                                        b5*logg*logg + \
                                                b6*logg*logg*logg + b7*FeH
    
    sig_log_R = 0.014
    
    starR = 10**(log_R)
    sig_starR = (sig_log_R/log_R) * starR
    
    return starR, sig_starR
    
if __name__ == '__main__':
    
    VI_Sun = 0.702
    sigVI_Sun = 0.010
    
    (log_theta_LD, sig_log_theta_LD) = calc_star_ang_radius(14.0,0.001,VI_Sun,sigVI_Sun,'V-I',Lclass='dwarfs')
    
    theta_LD = 10**(log_theta_LD)
    
    print('Log_10(theta_LD) = '+str(log_theta_LD)+' +/- '+str(sig_log_theta_LD))
    print('Theta_LD = '+str(theta_LD)+'mas')
    