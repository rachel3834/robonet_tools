import create_known_variable_table
import pytest

@pytest.mark.parametrize("test_input,catalog_type,expected", [
    ( {"OGLE-2017-BLG-0171": {
        "target_data": {
            "ra": 268.0280833333333,
            "dec": -30.179333333333336,
            "baseline_mag": 19.19,
            "spitzer_target": 'true',
            "MOA_alert_ID": "MOA-XXXX-BLG-123",
            "Gaia_alert_ID": 'null',
            "Gaia_alert_ra": 'null',
            "Gaia_alert_dec": 'null',
            "Gaia_alert_class": 'null',
            "Gaia_alert_comment": 'null',
            "ATel": 'null',
            "t0": 2457873.24,
            "tE": 144.8,
            "u0": 0.748
        },
        "search_radius_deg": 0.0005555555555555556,
        "rome_stars": [
            {
                "field_id": 177462,
                "ra": 268.0281140567904,
                "dec": -30.17935574282786,
                "quadrant": 2,
                "separation_deg": 3.4750110057902514e-05,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 0.0,
                "cal_i_mag_lsc_doma": 17.21795334366742
            },
            {
                "field_id": 177513,
                "ra": 268.0276490204527,
                "dec": -30.179511001820348,
                "quadrant": 2,
                "separation_deg": 0.00041536053827886856,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 18.259769153914483,
                "cal_i_mag_lsc_doma": 17.302972135689657
            }
        ]
    },
    "OGLE-2017-BLG-0196": {
        "target_data": {
            "ra": 267.89750000000004,
            "dec": -30.18041666666667,
            "baseline_mag": 19.565,
            "spitzer_target": 'false',
            "Gaia_alert_ID": 'null',
            "Gaia_alert_ra": 'null',
            "Gaia_alert_dec": 'null',
            "Gaia_alert_class": 'null',
            "Gaia_alert_comment": 'null',
            "ATel": 'null',
            "t0": 2457825.923,
            "tE": 28.76,
            "u0": 0.273,
            "KMTNet_alert_ID": "KMT-2017-BLG-0042",
            "KMTNet_ra": 267.89754166666665,
            "KMTNet_dec": -30.18043888888889,
            "KMTNet_t0": 2457825.81289,
            "KMTNet_tE": 26.83,
            "KMTNet_u0": 0.283
        },
        "search_radius_deg": 0.0005555555555555556,
        "rome_stars": [
            {
                "field_id": 177749,
                "ra": 267.8974351463135,
                "dec": -30.18038335746541,
                "quadrant": 2,
                "separation_deg": 6.52113020547056e-05,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 0.0,
                "cal_i_mag_lsc_doma": 17.749490268515952
            }
        ]
    }},
      'event',
        {
            177462: {
                    'ogle_event_id': 'OGLE-2017-BLG-0171',
                    'ogle_variable_id': None,
                    'moa_event_id': 'MOA-XXXX-BLG-123',
                    'kmtnet_event_id': None,
                    'spitzer_event': 'true',
                    'vvv_variable_id': None
        },
            177749: {
                'ogle_event_id': 'OGLE-2017-BLG-0196',
                'ogle_variable_id': None,
                'moa_event_id': None,
                'kmtnet_event_id': 'KMT-2017-BLG-0042',
                'spitzer_event': 'false',
                'vvv_variable_id': None
         }
        }
    ),
        (
        {"OGLE-BLG-ELL-008964": {
        "target_data": {
            "ra": 268.05783333333335,
            "dec": -30.02322222222222,
            "ogle_class": "ELL",
            "ogle_subclass": "ELL",
            "rome_field": 1,
            "VVV_name":'null',
            "VVV_class":'null',
            "VVV_period":'null',
            "Gaia_EDR3_ID":'null',
            "Gaia_alert_ID":'null',
            "Gaia_alert_ra":'null',
            "Gaia_alert_dec":'null',
            "Gaia_alert_class":'null',
            "Gaia_alert_comment":'null',
            "ATel": 'null'
        },
        "search_radius_deg": 0.0005555555555555556,
        "rome_stars": [
            {
                "field_id": 94690,
                "ra": 268.0575083843225,
                "dec": -30.022759487610063,
                "quadrant": 3,
                "separation_deg": 0.0005415537957165766,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 18.28222371210053,
                "cal_i_mag_lsc_doma": 17.15314500046374
            },
            {
                "field_id": 95042,
                "ra": 268.0578780379656,
                "dec": -30.0232492608699,
                "quadrant": 3,
                "separation_deg": 4.7215087879266345e-05,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 16.847932155451147,
                "cal_i_mag_lsc_doma": 14.898636594254745
            },
            {
                "field_id": 231185,
                "ra": 268.057448185379,
                "dec": -30.02348035424348,
                "quadrant": 3,
                "separation_deg": 0.00042170366519753133,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 17.40298638185171,
                "cal_i_mag_lsc_doma": 0.0
            },
            {
                "field_id": 361359,
                "ra": 268.05807867916633,
                "dec": -30.02360280036499,
                "quadrant": 3,
                "separation_deg": 0.0004358489983904456,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 0.0,
                "cal_i_mag_lsc_doma": 0.0
            }
        ]
    },
    " 370408": {
        "target_data": {
            "ra": 267.63882833223,
            "dec": -30.23672897454,
            "ogle_class": 'null',
            "ogle_subclass": 'null',
            "VVV_name": " 370408",
            "VVV_ra": 267.63882833223,
            "VVV_dec": -30.23672897454,
            "VVV_class": "EA/EB     ",
            "VVV_period": "   0.153427835195354773",
            "Gaia_EDR3_ID": "                   ",
            "Gaia_alert_ID":'null',
            "Gaia_alert_ra":'null',
            "Gaia_alert_dec":'null',
            "Gaia_alert_class":'null',
            "Gaia_alert_comment":'null',
            "ATel": 'null'
        },
        "search_radius_deg": 0.0005555555555555556,
        "rome_stars": [
            {
                "field_id": 207908,
                "ra": 267.63850549619366,
                "dec": -30.236511058376603,
                "quadrant": 1,
                "separation_deg": 0.00035395076563564815,
                "cal_g_mag_lsc_doma": 18.831220404942847,
                "cal_r_mag_lsc_doma": 18.159368481531637,
                "cal_i_mag_lsc_doma": 17.396697545559718
            },
            {
                "field_id": 280360,
                "ra": 267.6391630247741,
                "dec": -30.236790081497272,
                "quadrant": 1,
                "separation_deg": 0.0002955445128174325,
                "cal_g_mag_lsc_doma": 0.0,
                "cal_r_mag_lsc_doma": 0.0,
                "cal_i_mag_lsc_doma": 0.0
            }
        ]
    }
        },
        'variable',
        {
            95042: {
                    'ogle_event_id': None,
                    'ogle_variable_id': 'OGLE-BLG-ELL-008964',
                    'moa_event_id': None,
                    'kmtnet_event_id': None,
                    'spitzer_event': 'false',
                    'vvv_variable_id': None
        },
            280360: {
                    'ogle_event_id': None,
                    'ogle_variable_id': None,
                    'moa_event_id': None,
                    'kmtnet_event_id': None,
                    'spitzer_event': 'false',
                    'vvv_variable_id': ' 370408'
            }
        }
    )
])
def test_extract_rome_star_matches(test_input, catalog_type, expected):
    LUT = {}
    LUT = create_known_variable_table.extract_rome_star_matches(LUT, test_input, catalog_type=catalog_type)

    assert(LUT == expected)
