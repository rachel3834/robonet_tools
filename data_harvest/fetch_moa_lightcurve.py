from os import path
from sys import argv
import lightcurve_utils
from bs4 import BeautifulSoup
import requests

MOA_BASE_URL = 'https://www.massey.ac.nz/~iabond/moa/'

def get_moa_data(moa_event_id, event_year, output_file):
    """Function to scrap the information for a specific MOA event from their
    website, download their survey data and parse their lightcurve into
    TOM-compatible format"""

    event_url = path.join(MOA_BASE_URL, 'alert'+event_year, 'display.php?id='+moa_event_id)

    (page_response, soup) = scrap_moa_event_page(event_url)

    (ZP, flux_0) = find_moa_calibration(soup)

    lc_url = find_moa_lc_url(soup, event_year)

    local_file = download_moa_lightcurve(lc_url, path.dirname(output_file))
    print(local_file)

    parse_moa_lightcurve(local_file, output_file)

def scrap_moa_event_page(event_url):

    result = requests.get(event_url)

    if result.status_code == 200:
        soup = BeautifulSoup(result.text, 'html.parser')
    else:
        soup = None
        raise IOError('Warning: Got status='+str(result.status_code)+' from MOA for this event')

    return result, soup

def find_moa_calibration(soup):

    plist = soup.find_all('p')

    flux_0 = 0.0
    ZP = 27.40
    for entry in plist:
        if 'Calibration' in entry.text:
            substr = extract_text_between_strings(entry.text, '&#916F + ', ')')
            try:
                flux_0 = float(substr)
            except ValueError:
                flux_0 = 0.0

            substr = extract_text_between_strings(entry.text, 'I = ', '-')
            try:
                ZP = float(substr)
            except ValueError:
                ZP = 0.0

    return ZP, flux_0

def find_moa_lc_url(soup, event_year):

    plist = soup.find_all('p')

    lc_url = None
    for entry in plist:
        if 'fetchtxt.php' in entry.text:
            substr = extract_text_between_strings(entry.text, 'fetchtxt.php?path=', '>')

            lc_url = path.join(MOA_BASE_URL, 'alert'+event_year, 'fetchtxt.php?path=', substr)

    return lc_url

def extract_text_between_strings(full_line, str1, str2):
    """Function to isolate a length of a text string in full_line between a
    starting string str1 and and ending string str2"""

    i = full_line.find(str1) + len(str1)
    ii = full_line[i:].find(str2)

    return full_line[i:i+ii]

def download_moa_lightcurve(lc_url, output_dir):

    output_file = path.join(output_dir, path.basename(lc_url))

    result = requests.get(lc_url)

    if result.status_code == 200:
        soup = BeautifulSoup(result.text, 'html.parser')

        prelist = soup.find_all('pre')

        if len(prelist) == 1:
            f = open(output_file,'w')
            f.write(prelist[0].text)
            f.close()

        elif len(prelist) == 0:
            raise IOError('MOA lightcurve page appears to be empty')

        elif len(prelist) > 1:
            raise IOError('MOA lightcurve page ill-formatted')

    else:
        raise IOError('Unable to download lightcurve')

    return output_file

def parse_moa_lightcurve(input_file, output_file, ZP=27.4, flux_0=0.0):
    """Function to convert a lightcurve in MOAs standard format to a CSV
    format that can be uploaded to MOP"""

    lc = lightcurve_utils.read_moa_lightcurve(input_file, ZP=ZP, flux_0=flux_0)

    lightcurve_utils.output_tom_csv(lc, output_file)

if __name__ == '__main__':

    if len(argv) < 3:
        moa_event_id = input('Please enter MOA *internal* event ID: ')
        event_year = input('Please enter the year of detection: ')
        output_file = input('Please enter the path for the output file: ')
    else:
        moa_event_id = argv[1]
        event_year = argv[2]
        output_file = argv[3]

    get_moa_data(moa_event_id, event_year, output_file)
