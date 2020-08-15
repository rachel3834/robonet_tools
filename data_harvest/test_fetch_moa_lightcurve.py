import fetch_moa_lightcurve
import requests
from os import path

def test_scrap_moa_event_page():

    event_url = 'https://www.massey.ac.nz/~iabond/moa/alert2019/display.php?id=gb12-R-1-110418'

    (result,soup) = fetch_moa_lightcurve.scrap_moa_event_page(event_url)

    assert type(soup.get_text()) == type('a')

def test_find_moa_calibration():

    event_url = 'https://www.massey.ac.nz/~iabond/moa/alert2019/display.php?id=gb12-R-1-110418'

    (result,soup) = fetch_moa_lightcurve.scrap_moa_event_page(event_url)

    (ZP, flux_0) = fetch_moa_lightcurve.find_moa_calibration(soup)

    assert type(ZP) == type(0.0)
    assert type(flux_0) == type(0.0)
    assert ZP > 0.0
    assert flux_0 > 0.0

def test_find_moa_lc_url():

    event_url = 'https://www.massey.ac.nz/~iabond/moa/alert2019/display.php?id=gb12-R-1-110418'
    event_year = '2019'

    (result,soup) = fetch_moa_lightcurve.scrap_moa_event_page(event_url)

    lc_url = fetch_moa_lightcurve.find_moa_lc_url(soup,event_year)

    assert 'https://' in lc_url

def test_download_moa_lc():

    lc_url = 'https://www.massey.ac.nz/~iabond/moa/alert2019/fetchtxt.php?path=/moa/ephot/phot-gb12-R-1-110418.dat'

    output_file = fetch_moa_lightcurve.download_moa_lightcurve(lc_url, '.''')

    assert path.isfile(output_file)
    
if __name__ == '__main__':

    test_scrap_moa_event_page()
    test_find_moa_calibration()
    test_find_moa_lc_url()
    test_download_moa_lc()
