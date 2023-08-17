import numpy as np
import pytest

from pyLIMA.fits import LM_fit
from pyLIMA.fits import DE_fit
from pyLIMA.fits import TRF_fit
from pyLIMA.fits import MCMC_fit
from pyLIMA.models import PSPL_model
from pyLIMA.models import FSPL_model
from pyLIMA.models import USBL_model
from pyLIMA import event
from pyLIMA import telescopes

def EventFactory(name, ra, dec):

    # Instantiate an Event object
    e = event.Event(ra=ra,dec=dec)
    e.name = name

    # Generate a simulated lightcurve, and instantiate a corresponding
    # Telescope object
    # NOTE: This does not yet simulate an actual event lightcurve
    ndata = 100
    lc = np.zeros((ndata,3))
    lc[:,0] = np.linspace(2459000.0, 2459000.0+ndata, ndata)
    lc[:,1].fill(14.5)
    lc[:,2].fill(0.002)
    telescope1 = telescopes.Telescope(name=name+'_lc_1',
                                      camera_filter = 'I',
                                      light_curve = lc,
                                      light_curve_names = ['time','mag','err_mag'],
                                      light_curve_units = ['JD','mag','mag'])
    e.telescopes.append(telescope1)

    # Identify the lightcurve to be used as the survey reference, and
    # use PyLIMA's built in event verification method
    e.find_survey(name+'_lc_1')
    e.check_event()

    return e

def FitFactory(e, model_type='PSPL', fit_method='TRF'):
    """Function to generate simulated Fit objects, for event e"""

    if model_type == 'PSPL':
        m = PSPL_model.PSPLmodel(e)
    else:
        raise IOError('Unsupported model type ('+model_type+') requested')

    if fit_method == 'TRF':
        f = TRF_fit.TRFfit(m)
    else:
        raise IOError('Unsupported fit method ('+fit_method+') requested')

    return f

@pytest.mark.parametrize(
    "test, expected",
    [
        ({
          'name': 'Test1',
          'ra': 262.75616,
          'dec': -21.40123,
          'model_type': 'PSPL',
          'fit_method': 'TRF',
          't0': 2459334.0,
          'u0': 0.01,
          'tE': 30.0
          },
        False),
    ])
def test_verify_fit_config(test, expected):
    """Test the function to verify that a model fit is properly configured"""

    from pylima_verify import verify_fit_config

    e = EventFactory(test['name'], test['ra'], test['dec'])
    f = FitFactory(
                    e,
                    model_type=test['model_type'],
                    fit_method=test['fit_method']
                    )
    f.model_parameters_guess = [test['t0'], test['u0'], test['tE']]

    (status, message) = verify_fit_config(f)

    assert(status, expected)
