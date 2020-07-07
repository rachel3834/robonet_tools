import data_download
import log_utils
import framelist_utils

def test_is_frame_calibration_data():

    status = data_download.is_frame_calibration_data('cpt1m003-fa12-20200607-0123-e91.fits')
    assert status == False
    status = data_download.is_frame_calibration_data('cpt1m003-fa12-20200607-skyflat-bin2x2-B.fits')
    assert status == True
    status = data_download.is_frame_calibration_data('cpt1m003-fa12-20200607-dark-bin1x1.fits')
    assert status == True
    status = data_download.is_frame_calibration_data('cpt1m003-fa12-20200607-bias-bin1x1.fits')
    assert status == True

def test_build_frame_list():

    proposal = 'KEY2020B-003'
    config = {'log_dir': '.', 'proposal_ids': '[KEY2020B-003, LCO2020A-002]'}
    log = log_utils.start_day_log(config,'test_data_download')

    default_params = { 'url': 'https://archive.lco.global',
                        'DATE_OBS': '2020-07-06 12:30:00',
                        'PROPID': 'KEY2020B-003',
                        'INSTRUME': 'fa12',
                        'OBJECT': 'TEST',
                        'SITEID': 'cpt',
                        'TELID': '1m003',
                        'EXPTIME': 30.0,
                        'FILTER': 'gp',
                        'REQNUM': 2168224 }

    frame1 = default_params.copy()
    frame1['filename'] = 'cpt1m003-fa12-20200607-0123-e91.fits'

    frame2 = default_params.copy()
    frame2['filename'] = 'cpt1m003-fa12-20200607-skyflat-bin2x2-B.fits'

    frame3 = default_params.copy()
    frame3['filename'] = 'cpt1m003-fa12-20200607-dark-bin1x1.fits'

    frame4 = default_params.copy()
    frame4['filename'] = 'cpt1m003-fa12-20200607-bias-bin1x1.fits'

    query_results = {'results': [ frame1, frame2, frame3, frame4 ]}

    new_frames = []

    new_frames = framelist_utils.build_frame_list(config, query_results, proposal, new_frames, log)

    assert len(new_frames) == 1

    log_utils.close_log(log)
