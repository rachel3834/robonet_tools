import pyLIMA

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
import cycler
from matplotlib.ticker import MaxNLocator
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show

#import PyQt5   # necessary to display interactive matplotlib plots
#matplotlib.use('QtAgg')

from pyLIMA.fits import LM_fit
from pyLIMA.fits import DE_fit
from pyLIMA.fits import TRF_fit
from pyLIMA.fits import MCMC_fit
from pyLIMA.models import PSPL_model
from pyLIMA.models import FSPL_model
from pyLIMA.models import USBL_model, pyLIMA_fancy_parameters
from pyLIMA.outputs import pyLIMA_plots
from pyLIMA.toolbox import fake_telescopes, plots
from pyLIMA.parallax import parallax

from pyLIMA import event
from pyLIMA import telescopes


your_event = event.Event(ra=262.75616,dec=-21.40123)   # setting no parallax for starters
your_event.name = 'Gaia21bsg'

data_1 = np.loadtxt('data/star_20957_Gaia21bsg_fs01_ip_reduced.dat')
telescope_1 = telescopes.Telescope(name='Gaia_20957_i',
                                  camera_filter = 'I',
                                  light_curve = data_1.astype(float),
                                  light_curve_names = ['time','mag','err_mag'],
                                  light_curve_units = ['JD','mag','mag'])

data_2 = np.loadtxt('data/star_50085_Gaia21bsg_gp_reduced.dat')
telescope_2 = telescopes.Telescope(name='Gaia__50085_g',
                                  camera_filter = 'G',
                                  light_curve = data_2.astype(float),
                                  light_curve_names = ['time','mag','err_mag'],
                                  light_curve_units = ['JD','mag','mag'])

data_3 = np.loadtxt('data/star_79874_Gaia21bsg_ip_reduced.dat')
telescope_3 = telescopes.Telescope(name='Gaia_79874_i',
                                  camera_filter = 'I',
                                  light_curve = data_3.astype(float),
                                  light_curve_names = ['time','mag','err_mag'],
                                  light_curve_units = ['JD','mag','mag'])

use_ztf = False
if use_ztf:
    data_4 = np.loadtxt('data/reduced_ztf.dat',delimiter=' ')   # Removing ZTF data for normalization issues
    telescope_4 = telescopes.Telescope(name='ZTF_r',
                                  camera_filter = 'R',
                                  light_curve = data_4.astype(float),
                                  light_curve_names = ['time','mag','err_mag'],
                                  light_curve_units = ['JD','mag','mag'])

data_5 = np.loadtxt('data/reduced_gaia.dat',delimiter=' ')
telescope_5 = telescopes.Telescope(name='Gaia_g',
                                  camera_filter = 'G',
                                  light_curve = data_5.astype(float),
                                  light_curve_names = ['time','mag','err_mag'],
                                  light_curve_units = ['JD', 'mag','mag'])

your_event.telescopes.append(telescope_1)
your_event.telescopes.append(telescope_2)
your_event.telescopes.append(telescope_3)
if use_ztf:
    your_event.telescopes.append(telescope_4)   # remove this as well
your_event.telescopes.append(telescope_5)

your_event.find_survey('Gaia')
your_event.check_event()


# PSPL Model
pspl = PSPL_model.PSPLmodel(your_event)
fit1 = LM_fit.LMfit(pspl)
fit1.fit_parameters
fit1.fit()

print('RESULTS PSPL model: ',fit1.fit_results['best_model'])
expected_pspl_results = [2459421.7367247217, -1643.1186207426726, 1016.7157834731744, -358953.74753607647, 367486.29517801385, -22010.230307457714, 22314.805814874573, -114441.63859215772, 115814.69267379036, -271085.94159024704, 272601.577490539]
assert( fit1.fit_results['best_model']==expected_pspl_results )

#guess_parameters = fit1.fit_results['best_model']

#fit2 = MCMC_fit.MCMCfit(pspl)
#fit2.model_parameters_guess = guess_parameters
#fit2.fit()

guess_parameters = [ 2.45933417e+06,    # t0
                    1.08350472e-02,     # u0
                    2.48698969e+00,     # log(tE)
                    -2.42567297e+00,    # log(rho)
                    2.41845537e-01,     # log(s)
                    -1.63468426e-01,    # log(q)
                    7.44880142e-01]     # alpha

fancy = pyLIMA_fancy_parameters.standard_fancy_parameters
usbl = USBL_model.USBLmodel(your_event, fancy_parameters=fancy, parallax=['None',2.45935387e+06])
fit2 = TRF_fit.TRFfit(usbl)
fit2.model_parameters_guess = guess_parameters
fit2.fit()
print('RESULTS USBL model: ',fit2.fit_results['best_model'])
expected_usbl_results = [2459334.3781281007, 0.011964425775738204, 2.435014189162299, -2.460628968443897, 0.2360312692735136, -0.18262860172435824, 0.7918451957533834, 1070.2723059375824, 3165.635781440189, 31.8340425838355, 126.22459101315538, 178.95999566233235, 448.28602898231225, 961.3383664886596, 355.9217885862539]
assert( fit2.fit_results['best_model']==expected_usbl_results )