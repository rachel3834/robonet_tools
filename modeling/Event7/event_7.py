from pyLIMA.models import USBL_model
from pyLIMA.fits import TRF_fit
from pyLIMA.fits import MCMC_fit
from pyLIMA.models import PSBL_model
import multiprocessing as mul
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA.outputs import pyLIMA_plots
from pyLIMA.outputs import file_outputs
import pandas as pd
import pygtc

path_save = '/data/software/robonet_tools/modeling/Event_7/'
path_ephemerides = '/data/software/robonet_tools/modeling/Event_7/james_webb.txt'
path_model = '/data/software/robonet_tools/modeling/Event_7'

def bineadora(evento,fil, n):
    evento_sorted = evento[evento['band']==fil].sort_values(by='mjd')
    fs = evento_sorted['mjd'].values[0]
    tin = evento_sorted['mjd'].values[0]
    tfin = evento_sorted['mjd'].values[-1]
    bins = np.arange(tin,tfin+n,n)
    evento_sorted['binned'] = pd.cut(evento_sorted['mjd'], bins,labels=False)
    time = []
    mags = []
    errmags = []
    for i in range(len(bins)):
        if not len(evento_sorted[evento_sorted['binned']==i])==0:
            time.append(np.mean(evento_sorted[evento_sorted['binned']==i]['mjd']))
            mags.append(np.mean(evento_sorted[evento_sorted['binned']==i]['mag']))
            errmags.append(np.sqrt(sum(np.array(evento_sorted[evento_sorted['binned']==i]['magerr'])**2))/len(evento_sorted[evento_sorted['binned']==i]['magerr']))
    df0 = pd.DataFrame({
        'band': [fil]*len(errmags),
        'mjd': time,
        'mag': mags,
        'magerr': errmags})
    return df0

def chi2(mag,err,c):
    l = []
    for i in range(len(mag)):
        l.append((mag[i]-c)**2/err[i]**2)
    return np.sum(l)/len(mag)

def filtros(file_name,tit, N,binning):
    '''
    function with binnig
    Falta el filtro de eliminar las curvas que no tienen bandas de Rubin
    '''
    evento = pd.read_csv(file_name , sep = ',' , decimal = '.', skiprows = 20)
    params = pd.read_csv(file_name , sep = ':' , decimal = '.', skiprows = 0)
    
    filtercolor = {'w':'b','u':'c', 'g':'g', 'r':'y', 'i':'r', 'z':'m', 'y':'k'}
    mag_sat = {'w':14.8, 'u':14.7, 'g': 15.7, 'r': 15.8, 'i': 15.8, 'z': 15.3, 'y': 13.9}
    
    dict_params = {}
    for i in range(len(params[0:18].values[:,1])):
        dict_params[params[0:18].values[:,0][i]] = params[0:18].values[:,1][i]
#     print(dict_params)

    t0 = dict_params['t0']
    dict_params['te'] = dict_params.pop('tE')
    te = dict_params['te']
    
    #------------SI TIENE MAG MENOS DE LA MAG_SAT Y SI ESTA MAS DE 1 SIGMA----------------
    evento.loc[evento['band'] == 'w', 'm5'] = 25.9
    evento['mag_sat'] = evento['band'].map(mag_sat)
    criteria = evento['m5'] > evento['mag'] + 1*evento['magerr']
    criteria2 = evento['mag']>evento['mag_sat']
    filtered_evento = evento[criteria&criteria2]
    evento = filtered_evento
    
    #----SI TIENE MENOS DE 4 PUNTOS EN UNA BANDA, LA TIRO----------------    
    value_counts = evento['band'].value_counts()   #cuenta el número de filas(datos) en cada banda
    # Get the bands with less than 4 values
    classes_to_drop = value_counts[value_counts < 6].index.tolist()
    # Drop the classes with less than 4 values from the DataFrame
    evento = evento[~evento['band'].isin(classes_to_drop)]
    evento = evento.groupby('band').filter(lambda group: len(group) >= 10) #chequea que existan al menos 10 puntos por banda

    #----------------Separación para el bineado-----------------------------------
    only_event = evento.loc[(evento['mjd'] < t0+3.5*te) & (evento['mjd'] > t0-3.5*te)]
    only_event = only_event[['band','mjd','mag','magerr']]
    interval = [t0-5*te,t0+5*te]
    no_event = evento[~evento['mjd'].between(interval[0], interval[1])]
    #---------------Si No hay ni una banda de Rubin, tira el evento-------------------------------
    is_only_w = (evento['band'] == 'w').all()

    # Print the result
    if not is_only_w:
        crit_1 = {}
        if binning:
            p = 0
            tot_dot = []
            bin_dot = []
            for fil in ('w','u','g','r','i','z','y'):
                if fil in evento['band'].values:
                    if fil in no_event['band'].values:
                        bineado = bineadora(no_event,fil, N)
                        df = pd.concat([bineado, only_event[only_event['band']==fil]],ignore_index=True)            
                        mjd = df[df['band']==fil][['mjd','mag','magerr']].values[:,0]
                        mag = df[df['band']==fil][['mjd','mag','magerr']].values[:,1]
                        magerr = df[df['band']==fil][['mjd','mag','magerr']].values[:,2]
                        crit_1[fil] = np.c_[mjd,mag,magerr]
                        plt.errorbar(mjd,mag,magerr,linestyle=' ',color=filtercolor[fil],marker='.',alpha=1,capsize=2,label=str(fil)+': '+str(len(mjd)))
                        plt.legend()
                else:
                    crit_1[fil] = np.array([])
                p=p+1
            plt.gca().invert_yaxis()
            plt.xlabel('days [mjd]')
            plt.ylabel('magnitude')
            plt.axvspan(t0-te,t0+te, color='blue',alpha=0.2)
            plt.tight_layout()
            plt.show()
            return crit_1, dict_params
        else:
            chis=[]
            for fil in ('w','u','g','r','i','z','y'):
                if fil in evento['band'].values:
                        df = evento
                        mjd = df[df['band']==fil][['mjd','mag','magerr']].values[:,0]
                        mag = df[df['band']==fil][['mjd','mag','magerr']].values[:,1]
                        magerr = df[df['band']==fil][['mjd','mag','magerr']].values[:,2]
                        crit_1[fil] = np.c_[mjd,mag,magerr]
                        mag_baseline = np.mean(no_event['mag'][no_event['band']==fil])
                else:
                    crit_1[fil] = np.array([])
            is_value_greater_than_2 = False
            for value in chis:
                if value > 2:
                    is_value_greater_than_2 = True
                    break  # No need to continue checking once we find one value greater than 2

            if is_value_greater_than_2:
                for fil in ('w','u','g','r','i','z','y'):
                    if fil in evento['band'].values:
                            df = evento
                            mjd = df[df['band']==fil][['mjd','mag','magerr']].values[:,0]
                            mag = df[df['band']==fil][['mjd','mag','magerr']].values[:,1]
                            magerr = df[df['band']==fil][['mjd','mag','magerr']].values[:,2]
                            crit_1[fil] = np.c_[mjd,mag,magerr]

            return crit_1, dict_params
    else:
        return 0,0


from pyLIMA.fits import DE_fit

def fit_event_model(n,event_params,algo, model_type, wfirst_lc, lsst_u,lsst_g,lsst_r,lsst_i,lsst_z,lsst_y):
    tlsst = 60350.38482057137+2400000.5
    e = event.Event()
    e.name = 'Event_RR_'+str(int(n))
    e.ra = 270
    e.dec = -30
    tel_list = []

# Add a PyLIMA telescope object to the event with the Gaia lightcurve
    tel1 = telescopes.Telescope(name='Roman', camera_filter='W149',
                                     light_curve=wfirst_lc,
                                     light_curve_names = ['time','mag','err_mag'],
                                     light_curve_units = ['JD','mag','mag'],
                                     location='Space')
    
    ephemerides = np.loadtxt(path_ephemerides)
    ephemerides[:,0] = ephemerides[:,0]
    ephemerides[:,3] *=  60*300000/150000000
    deltaT = tlsst-ephemerides[:,0][0]
    ephemerides[:,0] = ephemerides[:,0]+deltaT
    tel1.spacecraft_positions ={'astrometry':[],'photometry':ephemerides}
    e.telescopes.append(tel1)
    
    tel_list.append('Roman')
        
    
    if len(lsst_u)!=0:    
        # Add a PyLIMA telescope object to the event with the LCO lightcurve
        tel2 = telescopes.Telescope(name='Rubin_u', camera_filter='u',
                                         light_curve=lsst_u,
                                         light_curve_names = ['time','mag','err_mag'],
                                         light_curve_units = ['JD','mag','mag'],
                                         location='Earth')#,

        e.telescopes.append(tel2)
        tel_list.append('Rubin_u')
    if len(lsst_g)!=0:    
    # Add a PyLIMA telescope object to the event with the LCO lightcurve
        tel3 = telescopes.Telescope(name='Rubin_g', camera_filter='g',
                                         light_curve=lsst_g,
                                         light_curve_names = ['time','mag','err_mag'],
                                         light_curve_units = ['JD','mag','mag'],
                                         location='Earth')#,
#                                          light_curve_magnitude_dictionnary={'time': 0, 'mag': 1, 'err_mag': 2},
#                                          clean_the_lightcurve=False)

        e.telescopes.append(tel3)
        tel_list.append('Rubin_g')
    if len(lsst_r)!=0:    
        # Add a PyLIMA telescope object to the event with the LCO lightcurve
        tel4 = telescopes.Telescope(name='Rubin_r', camera_filter='r',
                                         light_curve=lsst_r,
                                         light_curve_names = ['time','mag','err_mag'],
                                         light_curve_units = ['JD','mag','mag'],
                                         location='Earth')#,
#                                          light_curve_magnitude_dictionnary={'time': 0, 'mag': 1, 'err_mag': 2},
#                                          clean_the_lightcurve=False)

        e.telescopes.append(tel4)
        tel_list.append('Rubin_r')
    if len(lsst_i)!=0:    
        # Add a PyLIMA telescope object to the event with the LCO lightcurve
        tel5 = telescopes.Telescope(name='Rubin_i', camera_filter='i',
                                         light_curve=lsst_i,
                                         light_curve_names = ['time','mag','err_mag'],
                                         light_curve_units = ['JD','mag','mag'],
                                         location='Earth')#,


        e.telescopes.append(tel5)
        tel_list.append('Rubin_i')
    if len(lsst_z)!=0:    

        # Add a PyLIMA telescope object to the event with the LCO lightcurve
        tel6 = telescopes.Telescope(name='Rubin_z', camera_filter='z',
                                         light_curve=lsst_z,
                                         light_curve_names = ['time','mag','err_mag'],
                                         light_curve_units = ['JD','mag','mag'],
                                         location='Earth')

        e.telescopes.append(tel6)
        tel_list.append('Rubin_z')
    if len(lsst_y)!=0:    
        # Add a PyLIMA telescope object to the event with the LCO lightcurve
        tel7 = telescopes.Telescope(name='Rubin_y', camera_filter='y',
                                         light_curve = lsst_y,
                                         light_curve_names = ['time','mag','err_mag'],
                                         light_curve_units = ['JD','mag','mag'],
                                         location='Earth')

        e.telescopes.append(tel7)
        tel_list.append('Rubin_y')
    
    if len(e.telescopes) < 2:
        return 0,0,0
        
    else:

        e.check_event()

        psbl = PSBL_model.PSBLmodel(e, parallax=['Full', event_params['t0']])
        
        # Give the model initial guess values somewhere near their actual values so that the fit doesn't take all day
        lensing_parameters = [float(event_params['t0']), float(event_params['u0']), float(event_params['te']), 
                                  float(event_params['s']),float(event_params['q']),float(event_params['alpha']), float(event_params['piEN']), float(event_params['piEE'])]


        if algo == 'TRF':
            fit_2 = TRF_fit.TRFfit(psbl)
            fit_2.model_parameters_guess = lensing_parameters 
            rango = 1

            fit_2.fit_parameters['t0'][1] = [float(event_params['t0'])-abs(float(event_params['t0']))*rango,float(event_params['t0'])+abs(float(event_params['t0']))*rango] # t0 limits
            fit_2.fit_parameters['u0'][1] = [float(event_params['u0'])-abs(float(event_params['u0']))*rango,float(event_params['u0'])+abs(float(event_params['u0']))*rango] # u0 limits
            fit_2.fit_parameters['tE'][1] = [float(event_params['te'])-abs(float(event_params['te']))*rango,float(event_params['te'])+abs(float(event_params['te']))*rango] # logtE limits in days
            fit_2.fit_parameters['separation'][1] = [float(event_params['s'])-abs(float(event_params['s']))*rango,float(event_params['s'])+abs(float(event_params['s']))*rango] # logs limits
            fit_2.fit_parameters['mass_ratio'][1] = [float(event_params['q'])-abs(float(event_params['q']))*rango,float(event_params['q'])+abs(float(event_params['q']))*rango] # logq limits
            fit_2.fit_parameters['alpha'][1] = [float(event_params['alpha'])-abs(float(event_params['alpha']))*rango,float(event_params['alpha'])+abs(float(event_params['alpha']))*rango] # alpha limits (in radians)
            fit_2.fit_parameters['piEE'][1] = [float(event_params['piEE'])-abs(float(event_params['piEE']))*rango,float(event_params['piEE'])+abs(float(event_params['piEE']))*rango]
            fit_2.fit_parameters['piEN'][1]= [float(event_params['piEN'])-abs(float(event_params['piEN']))*rango,float(event_params['piEN'])+abs(float(event_params['piEN']))*rango]

            pool = mul.Pool(processes = 16)
            fit_2.fit()

            cov = fit_2.fit_results["covariance_matrix"]
            best_fit = np.array(fit_2.fit_results["best_model"])
            chichi = fit_2.fit_results["chi2"]
            time_fit = fit_2.fit_results['fit_time']
            true_values =  np.array(event_params)
            results = {'covariance_matrix':cov, 'best_model':best_fit, 'chi2':chichi, 'fit_time': time_fit, 'true_values':true_values}


            np.save(path_save+e.name+'_trf.npy', results)
            return fit_2, e, tel_list

            
        elif algo == 'MCMC':
            
            fit_2 = MCMC_fit.MCMCfit(psbl,MCMC_links=100000)

            fit_2.model_parameters_guess = lensing_parameters 
            rango = 1

            fit_2.fit_parameters['t0'][1] = [float(event_params['t0'])-abs(float(event_params['t0']))*rango,float(event_params['t0'])+abs(float(event_params['t0']))*rango] # t0 limits
            fit_2.fit_parameters['u0'][1] = [float(event_params['u0'])-abs(float(event_params['u0']))*rango,float(event_params['u0'])+abs(float(event_params['u0']))*rango] # u0 limits
            fit_2.fit_parameters['tE'][1] = [float(event_params['te'])-abs(float(event_params['te']))*rango,float(event_params['te'])+abs(float(event_params['te']))*rango] # logtE limits in days
            fit_2.fit_parameters['separation'][1] =[0,1]# [float(event_params['s'])-abs(float(event_params['s']))*rango,float(event_params['s'])+abs(float(event_params['s']))*rango] # logs limits
            fit_2.fit_parameters['mass_ratio'][1] = [float(event_params['q'])-abs(float(event_params['q']))*rango,float(event_params['q'])+abs(float(event_params['q']))*rango] # logq limits
            fit_2.fit_parameters['alpha'][1] = [float(event_params['alpha'])-abs(float(event_params['alpha']))*rango,float(event_params['alpha'])+abs(float(event_params['alpha']))*rango] # alpha limits (in radians)
            fit_2.fit_parameters['piEE'][1] = [float(event_params['piEE'])-abs(float(event_params['piEE']))*rango,float(event_params['piEE'])+abs(float(event_params['piEE']))*rango]
            fit_2.fit_parameters['piEN'][1]= [float(event_params['piEN'])-abs(float(event_params['piEN']))*rango,float(event_params['piEN'])+abs(float(event_params['piEN']))*rango]

            pool = mul.Pool(processes = 32)
            fit_2.fit(computational_pool = pool)

            chains = fit_2.fit_results["MCMC_chains"]
            chains_fluxes = fit_2.fit_results["MCMC_chains_with_fluxes"]
            ln_likelihood =fit_2.fit_results['ln(likelihood)'] 
            best_fit = np.array(fit_2.fit_results["best_model"])
            time_fit = fit_2.fit_results['fit_time']
            true_values =  np.array(event_params)
            results = {'MCMC_chains':chains, "MCMC_chains_with_fluxes" : chains_fluxes, 'ln_likelihood':ln_likelihood,'best_model':best_fit, 'fit_time': time_fit, 'true_values':true_values}


            np.save(path_save+e.name+'_mcmc.npy', results)
            return fit_2, e, tel_list

        elif algo == 'DE':

            fit_2 = DE_fit.DEfit(psbl)

            fit_2.model_parameters_guess = lensing_parameters 
            rango = 1

            fit_2.fit_parameters['t0'][1] = [float(event_params['t0'])-abs(float(event_params['t0']))*rango,float(event_params['t0'])+abs(float(event_params['t0']))*rango] # t0 limits
            fit_2.fit_parameters['u0'][1] = [float(event_params['u0'])-abs(float(event_params['u0']))*rango,float(event_params['u0'])+abs(float(event_params['u0']))*rango] # u0 limits
            fit_2.fit_parameters['tE'][1] = [float(event_params['te'])-abs(float(event_params['te']))*rango,float(event_params['te'])+abs(float(event_params['te']))*rango] # logtE limits in days
            fit_2.fit_parameters['separation'][1] = [float(event_params['s'])-abs(float(event_params['s']))*rango,float(event_params['s'])+abs(float(event_params['s']))*rango] # logs limits
            fit_2.fit_parameters['mass_ratio'][1] = [float(event_params['q'])-abs(float(event_params['q']))*rango,float(event_params['q'])+abs(float(event_params['q']))*rango] # logq limits
            fit_2.fit_parameters['alpha'][1] = [float(event_params['alpha'])-abs(float(event_params['alpha']))*rango,float(event_params['alpha'])+abs(float(event_params['alpha']))*rango] # alpha limits (in radians)
            fit_2.fit_parameters['piEE'][1] = [float(event_params['piEE'])-abs(float(event_params['piEE']))*rango,float(event_params['piEE'])+abs(float(event_params['piEE']))*rango]
            fit_2.fit_parameters['piEN'][1]= [float(event_params['piEN'])-abs(float(event_params['piEN']))*rango,float(event_params['piEN'])+abs(float(event_params['piEN']))*rango]

            pool = mul.Pool(processes = 32)
            fit_2.fit(computational_pool = pool)

            de_population = fit_2.fit_results['DE_population']
            ln_likelihood =fit_2.fit_results['-(ln_likelihood)'] 
            best_fit = np.array(fit_2.fit_results["best_model"])
            time_fit = fit_2.fit_results['fit_time']
            true_values =  np.array(event_params)

            results = {'DE_population':de_population, 'ln_likelihood':ln_likelihood,'best_model':best_fit, 'fit_time': time_fit, 'true_values':true_values}


            np.save(path_save+e.name+'_de.npy', results)
            return fit_2, e, tel_list


def fit_event_roman(n,event_params, algo,model_type, wfirst_lc):
    tlsst = 60350.38482057137+ 2400000.5
    e = event.Event(ra=270, dec=-30)
    
    e.name = 'Event_Roman_'+str(int(n))
#     e.ra = 270
#     e.dec = -30
    tel_list = []
    # Add a PyLIMA telescope object to the event with the Gaia lightcurve
    tel1 = telescopes.Telescope(name='Roman', camera_filter='W149',
#                                      spacecraft_name = 'Gaia',
                                     light_curve=wfirst_lc,
                                     light_curve_names = ['time','mag','err_mag'],
                                     light_curve_units = ['JD','mag','mag'],
                                     location='Space')
    # For spacecraft parallax, need to append the spacecraft_positions table here
    
    ephemerides = np.loadtxt(path_ephemerides)
    ephemerides[:,0] = ephemerides[:,0]
    ephemerides[:,3] *=  60*300000/150000000
    deltaT = tlsst-ephemerides[:,0][0]
    ephemerides[:,0] = ephemerides[:,0]+deltaT
    tel1.spacecraft_positions ={'astrometry':[],'photometry':ephemerides}
#     print()
    e.telescopes.append(tel1)
    tel_list.append('Roman')
    

    # Identify which dataset to use as the baseline:
    
    e.check_event()

    psbl = PSBL_model.PSBLmodel(e, parallax=['Full', event_params['t0']])


    # Give the model initial guess values somewhere near their actual values so that the fit doesn't take all day
    lensing_parameters = [float(event_params['t0']), float(event_params['u0']), float(event_params['te']), 
                              float(event_params['s']),float(event_params['q']),float(event_params['alpha']), float(event_params['piEN']), float(event_params['piEE'])]
    
    e.check_event()

    psbl = PSBL_model.PSBLmodel(e, parallax=['Full', event_params['t0']])

    # Give the model initial guess values somewhere near their actual values so that the fit doesn't take all day
    lensing_parameters = [float(event_params['t0']), float(event_params['u0']), float(event_params['te']), 
                              float(event_params['s']),float(event_params['q']),float(event_params['alpha']), float(event_params['piEN']), float(event_params['piEE'])]


    if algo == 'TRF':
        fit_2 = TRF_fit.TRFfit(psbl)
        fit_2.model_parameters_guess = lensing_parameters 
        rango = 1

        fit_2.fit_parameters['t0'][1] = [float(event_params['t0'])-abs(float(event_params['t0']))*rango,float(event_params['t0'])+abs(float(event_params['t0']))*rango] # t0 limits
        fit_2.fit_parameters['u0'][1] = [float(event_params['u0'])-abs(float(event_params['u0']))*rango,float(event_params['u0'])+abs(float(event_params['u0']))*rango] # u0 limits
        fit_2.fit_parameters['tE'][1] = [float(event_params['te'])-abs(float(event_params['te']))*rango,float(event_params['te'])+abs(float(event_params['te']))*rango] # logtE limits in days
        fit_2.fit_parameters['separation'][1] = [float(event_params['s'])-abs(float(event_params['s']))*rango,float(event_params['s'])+abs(float(event_params['s']))*rango] # logs limits
        fit_2.fit_parameters['mass_ratio'][1] = [float(event_params['q'])-abs(float(event_params['q']))*rango,float(event_params['q'])+abs(float(event_params['q']))*rango] # logq limits
        fit_2.fit_parameters['alpha'][1] = [float(event_params['alpha'])-abs(float(event_params['alpha']))*rango,float(event_params['alpha'])+abs(float(event_params['alpha']))*rango] # alpha limits (in radians)
        fit_2.fit_parameters['piEE'][1] = [float(event_params['piEE'])-abs(float(event_params['piEE']))*rango,float(event_params['piEE'])+abs(float(event_params['piEE']))*rango]
        fit_2.fit_parameters['piEN'][1]= [float(event_params['piEN'])-abs(float(event_params['piEN']))*rango,float(event_params['piEN'])+abs(float(event_params['piEN']))*rango]

        pool = mul.Pool(processes = 16)
        fit_2.fit()

        cov = fit_2.fit_results["covariance_matrix"]
        best_fit = np.array(fit_2.fit_results["best_model"])
        chichi = fit_2.fit_results["chi2"]
        time_fit = fit_2.fit_results['fit_time']
        true_values =  np.array(event_params)
        results = {'covariance_matrix':cov, 'best_model':best_fit, 'chi2':chichi, 'fit_time': time_fit, 'true_values':true_values}


        np.save(path_save+e.name+'_trf.npy', results)
        return fit_2, e, tel_list


    elif algo == 'MCMC':

        fit_2 = MCMC_fit.MCMCfit(psbl,MCMC_links=100000)

        fit_2.model_parameters_guess = lensing_parameters 
        rango = 1

        fit_2.fit_parameters['t0'][1] = [float(event_params['t0'])-abs(float(event_params['t0']))*rango,float(event_params['t0'])+abs(float(event_params['t0']))*rango] # t0 limits
        fit_2.fit_parameters['u0'][1] = [float(event_params['u0'])-abs(float(event_params['u0']))*rango,float(event_params['u0'])+abs(float(event_params['u0']))*rango] # u0 limits
        fit_2.fit_parameters['tE'][1] = [float(event_params['te'])-abs(float(event_params['te']))*rango,float(event_params['te'])+abs(float(event_params['te']))*rango] # logtE limits in days
        fit_2.fit_parameters['separation'][1] = [float(event_params['s'])-abs(float(event_params['s']))*rango,float(event_params['s'])+abs(float(event_params['s']))*rango] # logs limits
        fit_2.fit_parameters['mass_ratio'][1] = [float(event_params['q'])-abs(float(event_params['q']))*rango,float(event_params['q'])+abs(float(event_params['q']))*rango] # logq limits
        fit_2.fit_parameters['alpha'][1] = [float(event_params['alpha'])-abs(float(event_params['alpha']))*rango,float(event_params['alpha'])+abs(float(event_params['alpha']))*rango] # alpha limits (in radians)
        fit_2.fit_parameters['piEE'][1] = [float(event_params['piEE'])-abs(float(event_params['piEE']))*rango,float(event_params['piEE'])+abs(float(event_params['piEE']))*rango]
        fit_2.fit_parameters['piEN'][1]= [float(event_params['piEN'])-abs(float(event_params['piEN']))*rango,float(event_params['piEN'])+abs(float(event_params['piEN']))*rango]

        pool = mul.Pool(processes = 16)
        fit_2.fit(computational_pool = pool)

        chains = fit_2.fit_results["MCMC_chains"]
        chains_fluxes = fit_2.fit_results["MCMC_chains_with_fluxes"]
        ln_likelihood =fit_2.fit_results['ln(likelihood)'] 
        best_fit = np.array(fit_2.fit_results["best_model"])
        time_fit = fit_2.fit_results['fit_time']
        true_values =  np.array(event_params)
        results = {'MCMC_chains':chains, "MCMC_chains_with_fluxes" : chains_fluxes, 'ln_likelihood':ln_likelihood,'best_model':best_fit, 'fit_time': time_fit, 'true_values':true_values}

        np.save(path_save+e.name+'_mcmc.npy', results)

        return fit_2, e, tel_list
    
    elif algo == 'DE':

        fit_2 = DE_fit.DEfit(psbl)

        fit_2.model_parameters_guess = lensing_parameters 
        rango = 1

        fit_2.fit_parameters['t0'][1] = [float(event_params['t0'])-abs(float(event_params['t0']))*rango,float(event_params['t0'])+abs(float(event_params['t0']))*rango] # t0 limits
        fit_2.fit_parameters['u0'][1] = [float(event_params['u0'])-abs(float(event_params['u0']))*rango,float(event_params['u0'])+abs(float(event_params['u0']))*rango] # u0 limits
        fit_2.fit_parameters['tE'][1] = [float(event_params['te'])-abs(float(event_params['te']))*rango,float(event_params['te'])+abs(float(event_params['te']))*rango] # logtE limits in days
        fit_2.fit_parameters['separation'][1] = [float(event_params['s'])-abs(float(event_params['s']))*rango,float(event_params['s'])+abs(float(event_params['s']))*rango] # logs limits
        fit_2.fit_parameters['mass_ratio'][1] = [float(event_params['q'])-abs(float(event_params['q']))*rango,float(event_params['q'])+abs(float(event_params['q']))*rango] # logq limits
        fit_2.fit_parameters['alpha'][1] = [float(event_params['alpha'])-abs(float(event_params['alpha']))*rango,float(event_params['alpha'])+abs(float(event_params['alpha']))*rango] # alpha limits (in radians)
        fit_2.fit_parameters['piEE'][1] = [float(event_params['piEE'])-abs(float(event_params['piEE']))*rango,float(event_params['piEE'])+abs(float(event_params['piEE']))*rango]
        fit_2.fit_parameters['piEN'][1]= [float(event_params['piEN'])-abs(float(event_params['piEN']))*rango,float(event_params['piEN'])+abs(float(event_params['piEN']))*rango]

        pool = mul.Pool(processes = 32)
        fit_2.fit(computational_pool = pool)

        de_population = fit_2.fit_results['DE_population']
        ln_likelihood =fit_2.fit_results['-(ln_likelihood)'] 
        best_fit = np.array(fit_2.fit_results["best_model"])
        time_fit = fit_2.fit_results['fit_time']
        true_values =  np.array(event_params)
        
        results = {'DE_population':de_population, 'ln_likelihood':ln_likelihood,'best_model':best_fit, 'fit_time': time_fit, 'true_values':true_values}


        np.save(path_save+e.name+'_de.npy', results)
        return fit_2, e, tel_list
    
def fit_light_curve(file_name,binning,algo):
    '''
    This function take the name of the file that contains the light curves
    if pass all the selection criteria then fit the event and save the true parameters
    with the fited ones and his errors
    '''

    modeling_results = []
    events = []
    tels = []

    filtercolor = {'w':'b','u':'c', 'g':'g', 'r':'y', 'i':'r', 'z':'m', 'y':'k'}

    light_curve,PARAMS = filtros(file_name,1,100,binning)
    if not light_curve == 0:
        array_w149, array_u, array_g, array_r, array_i, array_z, array_y = light_curve['w'],light_curve['u'],light_curve['g'],light_curve['r'],light_curve['i'],light_curve['z'],light_curve['y']
        t0 = PARAMS['t0']
        te = PARAMS['te']
        model_type = 'PSBL'

        fit_2, e, tel_list = fit_event_model(PARAMS['Source'],PARAMS,algo, model_type, array_w149, array_u, array_g, array_r, array_i, array_z, array_y)
        fit_2, e, tel_list = fit_event_roman(PARAMS['Source'],PARAMS,algo, model_type, array_w149)
        return fit_2, e, tel_list
    else:
        print('This event is not selected and can´t be fitted.')
        return 0,0,0


i=int(7)

fit_light_curve(path_model+f'Event_{i}.txt', False, 'MCMC')


