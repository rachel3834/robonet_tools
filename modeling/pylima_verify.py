
# These declarations allow us to configure ranges for the key model parameters
# where the values are physically meaningful.  Values outside these ranges
# should raise a red flag
EXPECTATION_VALUES = {
            # Expected range of t0 should lie within the JD range of the
            # Rubin LSST survey
            't0': [2460796.5, 2464448.500]
}

# Since the order of the parameters depends on the model_type, we use a
# parameter_index to find the right list index for each parameter in the
# model_parameters_guess list
PARAMETER_INDEX = {
            't0': 0
}

def verify_fit_config(f):
    """Function to verify the configuration of a PyLIMA model fit"""

    status = True
    message = 'OK'

    for key, range in EXPECTATION_VALUES.items():
        # Extract the configured value from the list of initial guesses for
        # the parameter set.
        config_value = f.model_parameters_guess[PARAMETER_INDEX[key]]

        # Check that the configured value lies within the expected range
        if config_value < range[0] or config_value > range[1]:
            status = False
            message = key+' is outside expected range '+str(range[0])\
                        +' to '+str(range[1])

    return status, message
