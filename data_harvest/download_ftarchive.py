from sys import argv
from os import path
import requests
from datetime import datetime, timedelta

FTARCHIVEURL = 'http://ftarchive.lco.global/'

if len(argv) == 1:
    start_date = input('Please enter the UTC start date that you wish to access [YYYMMDD]:' )
    end_date = input('Please enter the UTC end date that you wish to access [YYYMMDD]:' )
    tel = input('Please enter the telescope whose data you need [FTN or FTS]: ')
    data_dir = input('Please enter a directory path for the downloaded FITS files: ')
else:
    start_date = argv[1]
    end_date = argv[2]
    tel = argv[3]
    data_dir = argv[4]

if 'FTN' in str(tel).upper():
    tel_extn = '#ft1'
elif 'FTS' in str(tel).upper():
    tel_extn = '#ft2'
else:
    raise IOError('Unrecognized telescope code given')

ts_start = datetime.strptime(start_date, "%Y%m%d")
ts_end = datetime.strptime(end_date, "%Y%m%d")
n_days = (ts_end-ts_start).days

date_range = [ts_start + timedelta(days=x) for x in range(n_days)]

for ts in date_range:
    url = path.join(FTARCHIVEURL,tel_extn, ts.strftime("%Y%m%d"))
    print(url)
    result = requests.get(url)
    print(result.text)
