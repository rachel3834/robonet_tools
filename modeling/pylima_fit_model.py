from os import path
import numpy as np
import argparse
from pyDANDIA import config_utils
from pyLIMA import event as pyEvent
from pyLIMA import telescopes
from pyLIMA import toolbox
from pyLIMA.fits import TRF_fit, DE_fit, MCMC_fit, LM_fit
from pyLIMA.fits import stats
from pyLIMA.models import PSPL_model, USBL_model, FSBL_model, pyLIMA_fancy_parameters
from pyLIMA.outputs import pyLIMA_plots
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import json
from datetime import datetime
import logging
import multiprocessing

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run_model_fit(config):
    logger.info('Starting to model ' + config['event_name'])

    # Create an event object to model
    event = create_event(config)

    # Load the lightcurve data and associate it with the event
    tel_list = load_lightcurves(config)
    for tel in tel_list:
        event.telescopes.append(tel)
    logger.info('Loaded lightcurve data')

    # Identify the lightcurve to be considered as the reference, and sanity check
    event.find_survey(config['reference_lightcurve_id'])
    event.check_event()

    # Create the selected type of model object
    model = create_model(config, event)

    # Create and configure the fitting method, and perform the fit
    fitter = create_model_fitter(config, model)
    if config['model_type'] in ['DE', 'MCMC']:
        pool = multiprocessing.Pool(processes=4)
        fitter.fit(computational_pool=pool)
    else:
        fitter.fit()

    # Extract and store model fit results
    model_params = gather_model_parameters(config, event, fitter)
    store_model_fit_params(config, model_params)

    pyLIMA_plots.plot_lightcurves(model, fitter.fit_results['best_model'])
    plt.show()

def create_event(config):

    s = SkyCoord(config['RA'], config['Dec'], frame='icrs', unit=(u.hourangle, u.deg))
    event = pyEvent.Event(ra=s.ra.deg, dec=s.dec.deg)
    event.name = config['event_name']

    return event

def load_lightcurves(config):

    tel_list = []
    for id, lc_file in config['data_files'].items():
        data = np.loadtxt(lc_file, comments='#', delimiter=' ')

        tel = telescopes.Telescope(name=id, camera_filter=id,
                               light_curve=data[:,[0,1,2]],
                               light_curve_names=['time', 'mag', 'err_mag'],
                               light_curve_units=['JD', 'mag', 'err_mag'])
        tel_list.append(tel)

    return tel_list

def create_model(config, event):

    if config['use_fancy_parameters']:
        fancy_set = pyLIMA_fancy_parameters.standard_fancy_parameters
    else:
        fancy_set = []

    if config['model_type'] == 'PSPL':
        model = PSPL_model.PSPLmodel(event, parallax=[config['parallax_type'], config['t0_par']])
    elif config['model_type'] == 'USBL':
        if config['use_fancy_parameters']:
            model = USBL_model.USBLmodel(event, fancy_parameters=fancy_set,
                                                parallax=[config['parallax_type'], config['t0_par']])
        else:
            model = USBL_model.USBLmodel(event, parallax=[config['parallax_type'], config['t0_par']])
    elif config['model_type'] == 'FSBL':
        model = FSBL_model.FSBLmodel(event, parallax=[config['parallax_type'], config['t0_par']])
    else:
        raise IOError('Configured model type ' + config['model_type'] + ' not recognised')
    model.define_model_parameters()

    if 'none' in str(config['parallax_type']).lower():
        logger.info('Applying a ' + config['model_type'] + ' model without parallax')
    else:
        logger.info('Applying a ' + config['model_type'] + ' model with parallax '
                        + config['parallax_type'] + ' and t0_par=' + str(config['t0_par']))

    return model

def create_model_fitter(config, model):

    if config['fit_type'] == 'LM':
        fitter = LM_fit.LMfit(model)
    elif config['fit_type'] == 'TRF':
        fitter = TRF_fit.TRFfit(model, loss_function=config['loss_function'])
    elif config['fit_type'] == 'DE':
        if config['use_fancy_parameters']:
            fitter = DE_fit.DEfit(model, DE_population_size = 20, max_iteration = 10000,
                                                display_progress = True,)
        else:
            fitter = DE_fit.DEfit(model, DE_population_size = 20, max_iteration = 10000,
                                                display_progress = True)
    elif config['fit_type'] == 'MCMC':
        fitter = MCMC_fit.MCMCfit(model)
    else:
        raise IOError('Configured fit type ' + config['fit_type'] + ' not recognised')
    logger.info('Using the ' + config['fit_type'] + ' method')

    if config['use_boundaries']:
        fitter.fit_parameters["t0"][1] = [config['t0_min'], config['t0_max']]
        fitter.fit_parameters["tE"][1] = [config['tE_min'], config['tE_max']]
        fitter.fit_parameters["u0"][1] = [config['u0_min'], config['u0_max']]

    if config['use_initial_guess']:
        fitter.model_parameters_guess = config['guess']

    return fitter

def load_config():

    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help='Path to configuration file')
    args = parser.parse_args()

    config = config_utils.read_config(args.config_file)

    # Parse modeling-specific parameters
    boolean_keys = ['use_boundaries', 'use_initial_guess', 'use_fancy_parameters']
    for key in boolean_keys:
        if 'true' in str(config[key]).lower():
            config[key] = True
        else:
            config[key] = False

    float_keys = ['t0_min', 't0_max', 'u0_min', 'u0_max', 'tE_min', 'tE_max', 't0_par']
    for key in float_keys:
        config[key] = float(config[key])

    parse_initial_guess(config)

    return config

def parse_initial_guess(config):

    if config['model_type'] == 'PSPL' and 'none' in str(config['parallax_type']).lower():
        config['guess'] = [config['initial_guess_parameters']['t0'],
                           config['initial_guess_parameters']['u0'],
                           config['initial_guess_parameters']['tE']]

    elif config['model_type'] == 'USBL' and 'none' in str(config['parallax_type']).lower():
        config['guess'] = [config['initial_guess_parameters']['t0'],
                           config['initial_guess_parameters']['u0'],
                           config['initial_guess_parameters']['tE'],
                           config['initial_guess_parameters']['rho'],
                           config['initial_guess_parameters']['s'],
                           config['initial_guess_parameters']['q'],
                           config['initial_guess_parameters']['alpha']]

        # If fancy parameters are selected, order is: [t0, u0, log_tE, log_rho, log_s, log_q, alpha]
        if config['use_fancy_parameters']:
            config['guess'] = [config['initial_guess_parameters']['t0'],
                               config['initial_guess_parameters']['u0'],
                               np.log10(config['initial_guess_parameters']['tE']),
                               np.log10(config['initial_guess_parameters']['rho']),
                               np.log10(config['initial_guess_parameters']['s']),
                               np.log10(config['initial_guess_parameters']['q']),
                               config['initial_guess_parameters']['alpha']]

    return config

def gather_model_parameters(config, event, fitter):
    """
    Function to gather the parameters of a PyLIMA fitted model into a dictionary for easier handling.
    """

    # PyLIMA model objects store the fitted values of the model parameters in the fit_results attribute,
    # which is a list of the values pertaining to the model used for the fit.  Since this model can have a
    # variable number of parameters depending on which type of model is used, we use the fit object's built-in
    # list of key indices
    param_keys = list(fitter.fit_parameters.keys())

    model_params = {}

    for i, key in enumerate(param_keys):
        if key in ['t0' 'tE']:
            ndp = 3
        else:
            ndp = 5
        model_params[key] = np.around(fitter.fit_results["best_model"][i], ndp)
        if config['fit_type'] in ['TRF']:
            model_params[key+'_error'] = np.around(np.sqrt(fitter.fit_results["covariance_matrix"][i,i]), ndp)

    # model_params['chi2'] = np.around(fitter.fit_results["best_model"][-1], 3)
    # Reporting actual chi2 instead value of the loss function
    (chi2, pyLIMA_parameters) = fitter.model_chi2(fitter.fit_results["best_model"])
    model_params['chi2'] = np.around(chi2, 3)

    # If the model did not include parallax, zero those parameters
    if 'piEN' not in param_keys:
        model_params['piEN'] = 0.0
        model_params['piEN_error'] = 0.0
        model_params['piEE'] = 0.0
        model_params['piEE_error'] = 0.0

    # Calculate the reduced chi2
    ndata = 0
    for i,tel in enumerate(event.telescopes):
        ndata += len(tel.lightcurve_magnitude)
    model_params['red_chi2'] = np.around(model_params['chi2'] / float(ndata - len(param_keys)),3)

    if config['fit_type'] in ['TRF']:
        model_params['Fit_covariance'] = fitter.fit_results["covariance_matrix"]

    # Parse the fit parameters because some of these items contain
    # array values that JSON can't serialize
    fit_params = []
    for key, value in fitter.fit_parameters.items():
        print('1 ',key, value)
        if type(value[1]) == type(np.zeros(1)):
            value = [key, [value[0], value[1].tolist()]]
        else:
            value = [key, [value[0], value[1]]]
        print('2 ',key, value)
        fit_params.append(value)
    model_params['fit_parameters'] = fit_params

    # Calculate fit statistics
    # The fitter.model_residuals returns photometric and astrometric residuals as a dictionary
    # while the photometric residuals provides a list of arrays consisting of the
    # photometric residuals, photometric errors, and error_flux
    try:
        res = fitter.model_residuals(fitter.fit_results['best_model'])
        sw_test = stats.normal_Shapiro_Wilk(
            (np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])))
        model_params['SW_test'] = np.around(sw_test[0],3)
        ad_test = stats.normal_Anderson_Darling(
            (np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])))
        model_params['AD_test'] = np.around(ad_test[0],3)
        ks_test = stats.normal_Kolmogorov_Smirnov(
            (np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])))
        model_params['KS_test'] = np.around(ks_test[0],3)
        model_params['chi2_dof'] = np.sum((np.ravel(res[0]['photometry'][0]) / np.ravel(res[1]['photometry'][0])) ** 2) / (
                len(np.ravel(res[0]['photometry'][0])) - 5)
    except:
        model_params['SW_test'] = np.nan
        model_params['AD_test'] = np.nan
        model_params['KS_test'] = np.nan
        model_params['chi2_dof'] = np.nan

    if config['fit_type'] == 'DE':
        output_file = config['model_type'] + '_' + config['fit_type'] + '_population.npy'
        np.save(path.join(config['output_dir'], output_file), fitter.fit_results['DE_population'])
    elif config['fit_type'] == 'MCMC':
        output_file = config['model_type'] + '_' + config['fit_type'] + '_chains.npy'
        np.save(path.join(config['output_dir'], output_file), fitter.fit_results['MCMC_chains'])

    print(model_params)

    return model_params

def store_model_fit_params(config, model_params):

    # Serialize the convarience array first, since JSON can't handle this otherwise
    if config['fit_type'] in ['TRF']:
        model_params['Fit_covariance'] = model_params['Fit_covariance'].tolist()

    # Add the input parameters for the record:
    for key, value in config.items():
        model_params[key] = value

    # Add a timestamp
    now = datetime.utcnow()
    model_params['timestamp'] = now.strftime('%Y-%m-%dT%H:%M:%S') + 'UTC'

    # Convert to JSON:
    json_object = json.dumps(model_params, indent=4)

    # Output to file
    file_path = path.join(config['output_dir'], config['model_type'] + '_fit.json')
    with open(file_path, "w") as outfile:
        outfile.write(json_object)

def flux_to_mag(flux):

    ZP_pyLIMA = 27.4
    magnitude = ZP_pyLIMA - 2.5 * np.log10(flux)
    return magnitude

def fluxerror_to_magerror(flux, flux_error):

    mag_err = (2.5 / np.log(10.0)) * flux_error / flux
    return mag_err

def mag_to_flux(mag):

    ZP_pyLIMA = 27.4
    flux = 10**((mag - ZP_pyLIMA) / -2.5)

    return flux

if __name__ == '__main__':
    config = load_config()
    run_model_fit(config)


