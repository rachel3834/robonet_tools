import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u

parser = argparse.ArgumentParser()
parser.add_argument('ra', help='RA as sexigesimal colon-separated or decimal degrees')
parser.add_argument('dec', help='Dec as sexigesimal colon-separated or decimal degrees')
args = parser.parse_args()

if ':' in args.ra:
    s = SkyCoord(args.ra, args.dec, frame='icrs', unit=(u.hourangle, u.deg))
else:
    s = SkyCoord(args.ra, args.dec, frame='icrs', unit=(u.deg, u.deg))

sgal = s.transform_to('galactic')
print(sgal)