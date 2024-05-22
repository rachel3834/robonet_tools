import fieldphot_to_ipaclc
import crossmatch_to_ipactable
import pytest
import numpy as np
from pyDANDIA import crossmatch
from astropy.table import Table, Column
from os import path, remove
from astropy.io import fits

def mock_field_photometry(nstars, nimages):
    """
    Function to create a mock field photometry data array, in PyDANDIA format, for unit testing purposes.
    """

    tzero = 2560400.0
    ap_radius = 5.0
    mean_mags = np.linspace(15.0, 22.0, nstars)
    ZP = 25.0

    phot_data = np.zeros((nstars, nimages, 17))
    for i in range(0,nimages,1):
        phot_data[:,i,0].fill(tzero + i)                                        # HJD
        phot_data[:,i,14] = np.random.uniform(low=0.0, high=100.0, size=nstars)  # Delta x
        phot_data[:,i,15] = np.random.uniform(low=0.0, high=100.0, size=nstars)  # Delta y
        phot_data[:,i,9].fill(1.0)                                             # Phot scale factor
        phot_data[:,i,10].fill(1.0)                                             # Phot scale factor uncertainty


    for j in range(0,nstars,1):
        phot_data[j,:,1] = np.random.normal(loc=mean_mags[j], scale=1.0, size=nimages) # Inst mag
        phot_data[j,:,2] = abs(np.random.normal(loc=0.01, scale=1.0, size=nimages))         # Inst mag error
        phot_data[j,:,3] = np.random.normal(loc=mean_mags[j], scale=1.0, size=nimages) # Calib mag
        phot_data[j,:,4] = abs(np.random.normal(loc=0.01, scale=1.0, size=nimages))         # Calib mag error
        phot_data[j,:,5] = np.random.normal(loc=mean_mags[j], scale=1.0, size=nimages) # Normalized mag
        phot_data[j,:,6] = abs(np.random.normal(loc=0.01, scale=1.0, size=nimages))         # Normalized mag error
        phot_data[j,:,7] = np.random.normal(loc=mean_mags[j], scale=1.0, size=nimages) # Corrected mag
        phot_data[j,:,8] = abs(np.random.normal(loc=0.01, scale=1.0, size=nimages))         # Corrected mag error
    # Column 16 is the qc_flag, left at zero for no error

    return phot_data

def mock_xmatch(nstars, nimages):
    """
    Function to create a mock crossmatch object for unittesting purposes.
    """
    tzero = 2560400.0

    xmatch = crossmatch.CrossMatchTable()
    xmatch.gaia_dr = 'Gaia_EDR3'

    # Create table structure and add datasets
    params = {'datasets': {
        'ROME-FIELD-01_lsc-doma-1m0-05-fa15_ip': ['primary_ref', '/data/testing_path/ROME-FIELD-01_lsc-doma-1m0-05-fa15_ip', None],
        'ROME-FIELD-01_lsc-doma-1m0-05-fa15_rp': ['no_ref', '/data/testing_path/ROME-FIELD-01_lsc-doma-1m0-05-fa15_rp', None],
        'ROME-FIELD-01_lsc-doma-1m0-05-fa15_gp': ['no_ref', '/data/testing_path/ROME-FIELD-01_lsc-doma-1m0-05-fa15_gp', None],
    },
            'primary_ref': 'ROME-FIELD-01_lsc-doma-1m0-05-fa15_ip'}
    xmatch.create(params)

    # Populate the field index
    ncol = len(xmatch.field_index.colnames)
    for star in range(0,nstars,1):
        xmatch.field_index.add_row(
            [
                star + 1,                                           # Field ID
                np.random.normal(loc=267.0, scale=5.0, size=1),     # RA
                np.random.normal(loc=-29.5, scale=1.0, size=1),     # Dec
                np.random.randint(1, high=5, size=1),               # Quadrant
                star + 1,                                           # Quadrant ID mimicing field ID
                str(4056436121079690000 + star),                    # Gaia ID
                star + 1                                            # Star ID in first dataset
            ] + ([0] * (ncol - 7))                                  # Star ID in other datasets
        )


    # Initialize the stars table
    xmatch.init_stars_table()

    # Populate the images table
    rints = np.random.randint(0, high=len(params['datasets'].keys()), size=nimages)
    datasets = list(params['datasets'].keys())
    filters = ['gp', 'rp', 'ip']
    for i in range(0,nimages,1):
        xmatch.images.add_row(
            [
                i,                      # index
                'test_file',            # filename
                datasets[rints[i]],     # dataset_code
                filters[rints[i]],      # filter
                tzero + i,              # HJD
                '2025-05-22T00:00:00',  # Datetime
                300.0,                  # Exposure
                '18:00:00.0',           # RA
                '-27:30:00.0',          # Dec
                20.0,                   # Moon angular separation
                0.5,                    # Moon fraction
                np.random.uniform(low=1.0, high=1.8, size=1),   # Airmass
                0.0, 0.0,               # Sigma x, y
                2000.0,                 # Sky background
                2000.0,                 # Median sky
                4.0,                    # FWHM
                0.0,                    # Corr_xy
                nstars,                 # Nstars
            ] + [0.0]*29
        )

    return xmatch

def extract_test_lc(phot_data, xmatch, quad_idx):
    """
    Function to extract the expected lightcurve data from the synthesized test dataset
    """

    lc = {}
    datacount = {}

    # Build indices of which images are in which filter
    filter_set = {'gp': None, 'rp': None, 'ip': None}
    for f in filter_set.keys():
        filter_set[f] = np.where(xmatch.images['filter'] == f)[0]
        shortcodes = [xmatch.get_dataset_shortcode(x) for x in xmatch.images['dataset_code'][filter_set[f]]]

        lc[f] = Table([
                Column(name='HJD', data=phot_data[quad_idx,filter_set[f],0]),
                Column(name='inst_mag', data=phot_data[quad_idx,filter_set[f],1]),
                Column(name='inst_mag_error', data=phot_data[quad_idx,filter_set[f],2]),
                Column(name='calib_mag', data=phot_data[quad_idx,filter_set[f],3]),
                Column(name='calib_mag_error', data=phot_data[quad_idx,filter_set[f],4]),
                Column(name='norm_mag', data=phot_data[quad_idx,filter_set[f],7]),
                Column(name='norm_mag_error', data=phot_data[quad_idx,filter_set[f],8]),
                Column(name='qc_flag', data=phot_data[quad_idx,filter_set[f],16]),
                Column(name='dataset', data=shortcodes),
                Column(name='airmass', data=xmatch.images['airmass'][filter_set[f]]),
                Column(name='moon_frac', data=xmatch.images['moon_fraction'][filter_set[f]]),
                Column(name='moon_sep', data=xmatch.images['moon_ang_separation'][filter_set[f]]),
                Column(name='sky_bkgd', data=xmatch.images['sky'][filter_set[f]]),
                Column(name='fwhm', data=xmatch.images['fwhm'][filter_set[f]]),
            ])

        lc[f].sort(['HJD'])
        datacount[f] = len(lc[f])

    return lc, datacount

def test_fetch_photometry_for_star():
    # Synthesize a small test dataset in the right format
    phot_data = mock_field_photometry(100, 10)
    xmatch = mock_xmatch(100, 10)

    # Select a star for testing and generate the expected output
    quad_idx = 10
    test_lc, test_datacount = extract_test_lc(phot_data, xmatch, quad_idx)

    # Apply function to be tested and compare output
    lc, datacount = fieldphot_to_ipaclc.fetch_photometry_for_star(quad_idx, xmatch, phot_data)

    # Define tolerances for comparing floating point values in the lightcurves:
    tols = {
        'HJD': 1.0,
        'inst_mag': 1e-4,
        'inst_mag_error': 1e-4,
        'calib_mag': 1e-4,
        'calib_mag_error': 1e-4,
        'norm_mag': 1e-4,
        'norm_mag_error': 1e-4,
        'qc_flag': 1,
        'airmass': 1e-2,
        'moon_frac': 0.1,
        'moon_sep': 0.1,
        'sky_bkgd': 1.0,
        'fwhm': 1e-3,
    }

    exact_cols = ['dataset']

    assert(lc.keys() == test_lc.keys())
    for f, data in test_lc.items():
        print('TEST DATA: ', f, data)
        print('RETURNED DATA: ',f, lc[f])
        for col, rtol in tols.items():
            np.testing.assert_allclose(lc[f][col], data[col], rtol=rtol)
        for col in exact_cols:
            print('COL: ',col, data[col])
            assert((lc[f][col] == data[col]).all())
        assert(test_datacount[f] == datacount[f])

class CodeArguments():
    def __init__(self):
        self.field_name = 'ROME-FIELD-01'
        self.qid = 1
        self.output_dir = './test_output'

def test_output_to_ipac_lightcurve():

    # Simulate code arguments:
    args = CodeArguments()

    # Synthesize a small test dataset in the right format
    phot_data = mock_field_photometry(100, 10)
    xmatch = mock_xmatch(100, 10)
    star_field_id = 50
    star_field_idx = star_field_id - 1
    str_field_id = crossmatch_to_ipactable.bin_num_string(star_field_id)

    quad_id = 10
    test_lc, test_datacount = extract_test_lc(phot_data, xmatch, quad_id-1)

    expected_header_keys = {
        'NAME': args.field_name + '_' + str_field_id,
        'FIELD': args.field_name,
        'FIELD_ID': star_field_id,
        'QUADRANT': args.qid,
        'RA': xmatch.stars['ra'][star_field_idx],
        'DEC': xmatch.stars['dec'][star_field_idx],
        'GAIA_ID': xmatch.stars['gaia_source_id'][star_field_idx],
        'GAIACAT': 'Gaia_EDR3',
        'NDATA_G': len(test_lc['gp']),
        'NDATA_R': len(test_lc['rp']),
        'NDATA_I': len(test_lc['ip']),
    }
    expected_table_cols = [
        'HJD',
        'inst_mag', 'inst_mag_error',
        'calib_mag', 'calib_mag_error',
        'norm_mag', 'norm_mag_error',
        'qc_flag',
        'dataset',
        'airmass', 'moon_frac', 'moon_sep',
        'sky_bkgd', 'fwhm'
    ]
    hdu_map = {
        'LIGHTCURVE_SDSS_G': 'gp',
        'LIGHTCURVE_SDSS_R': 'rp',
        'LIGHTCURVE_SDSS_I': 'ip',
    }

    file_path = fieldphot_to_ipaclc.output_to_ipac_lightcurve(args, star_field_id, quad_id, xmatch, test_lc)

    assert(path.isfile(file_path))
    with fits.open(file_path) as hdul:
        # Check the primary header data is correct
        for key, value in expected_header_keys.items():
            assert(hdul[0].header[key] == value)

        # Check we have 3 tables plus the primary data unit, and that
        # each table has the right columns and number of entries
        assert(len(hdul) == 4)
        for hdu in hdul[1:]:
            for col in hdu.data.columns:
                assert(col.name in expected_table_cols)
            assert(len(hdu.data) == len(test_lc[hdu_map[hdu.name]]))

    # Force a test where one of the lightcurves is zero-length
    remove(file_path)
    test_lc['rp'] = Table([
                Column(name='HJD', data=np.array([])),
                Column(name='inst_mag', data=np.array([])),
                Column(name='inst_mag_error', data=np.array([])),
                Column(name='calib_mag', data=np.array([])),
                Column(name='calib_mag_error', data=np.array([])),
                Column(name='norm_mag', data=np.array([])),
                Column(name='norm_mag_error', data=np.array([])),
                Column(name='qc_flag', data=np.array([])),
                Column(name='dataset', data=np.array([])),
                Column(name='airmass', data=np.array([])),
                Column(name='moon_frac', data=np.array([])),
                Column(name='moon_sep', data=np.array([])),
                Column(name='sky_bkgd', data=np.array([])),
                Column(name='fwhm', data=np.array([])),
            ])
    file_path = fieldphot_to_ipaclc.output_to_ipac_lightcurve(args, star_field_id, quad_id, xmatch, test_lc)

    with fits.open(file_path) as hdul:
        for hdu in hdul[1:]:
            if hdu.name == 'LIGHTCURVE_SDSS_R':
                for col in hdu.data.columns:
                    assert (col.name in expected_table_cols)
                assert (len(hdu.data) == 0)
    remove(file_path)