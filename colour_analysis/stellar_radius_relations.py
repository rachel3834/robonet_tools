# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 13:56:20 2018

@author: rstreet
"""
import numpy as np

def calc_star_ang_radius(Q,sigQ,PQ,sigPQ,pqcolour,Lclass=None):
    """Function to calculate the limb-darkened angular radius of a star, 
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
    """
    
    coeffs = fetch_coefficients(pqcolour,Lclass=Lclass)
    
    log_theta_Q0 = 0.0
    var_log_theta_Q0 = sigPQ * sigPQ
    
    for n in range(0,3,1):
        log_theta_Q0 += coeffs['c'][n] * PQ**n
        var_log_theta_Q0 += coeffs['sig_c'][n]*coeffs['sig_c'][n]
    
    log_theta_LD = log_theta_Q0 - 0.2 * Q
    
    sig_log_theta_LD = np.sqrt( var_log_theta_Q0 + ( sigQ*sigQ ) )
    
    return log_theta_LD, sig_log_theta_LD


def fetch_coefficients(pqcolour,Lclass=None):
    """Function to return the correct set of co-efficients for the
    stellar angular radius calculation.
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


if __name__ == '__main__':
    
    VI_Sun = 0.702
    sigVI_Sun = 0.010
    
    (log_theta_LD, sig_log_theta_LD) = calc_star_ang_radius(14.0,0.001,VI_Sun,sigVI_Sun,'V-I',Lclass='dwarfs')
    
    theta_LD = 10**(log_theta_LD)
    
    print('Log_10(theta_LD) = '+str(log_theta_LD)+' +/- '+str(sig_log_theta_LD))
    print('Theta_LD = '+str(theta_LD)+'rads = '+str(theta_LD*(180.0/np.pi))+'deg')
    