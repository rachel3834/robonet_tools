from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def load_KMTNet_fields():
    """
    Returns a numpy array with the vertices
    describing the KMTNet zone polygon.
    """
    fields = np.array([[264.00, -37.40],
                       [270.50, -37.40],
                       [270.50, -33.75],
                       [272.50, -33.75],
                       [272.50, -32.00],
                       [275.50, -32.00],
                       [275.50, -25.30],
                       [275.60, -25.30],
                       [275.60, -21.90],
                       [272.00, -21.90],
                       [272.00, -23.00],
                       [270.40, -23.00],
                       [270.40, -20.50],
                       [264.50, -20.50],
                       [264.50, -22.70],
                       [262.00, -22.70],
                       [262.00, -26.25],
                       [260.50, -26.25],
                       [260.50, -31.40],
                       [262.00, -31.40],
                       [262.00, -36.00],
                       [264.00, -36.00]])
    return fields


def event_not_in_OMEGA_II(ra, dec, KMTNet_fields):
    """
    This function checks if the event is within the KMTNet fields.
    If it is, the event has to be rejected and not followed by OMEGA-II.

    :param ra: Right Ascention of the event.
    :param dec: Declination of the event.
    :param KMTNet_fields: A numpy array that contains a series of
                          points describing roughly
    :return: Boolean value if the event is not ok for OMEGA II.
    """

    exclusion_zone = Polygon(zip(KMTNet_fields[:, 0], KMTNet_fields[:, 1]))

    not_in_OMEGA_II = False

    point = Point(ra, dec)
    if (exclusion_zone.contains(point)):
        not_in_OMEGA_II = True

    return not_in_OMEGA_II

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('ra', type=str, help='RA in sexigesimal')
    parser.add_argument('dec', type=str, help='DEC in sexigesimal, substituting n for -ve')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    print(args)
    s = SkyCoord(args.ra, args.dec.replace('n','-'), frame='icrs', unit=(u.hourangle, u.deg))
    kmtnet_fields = load_KMTNet_fields()
    event_status = event_not_in_OMEGA_II(s.ra.deg, s.dec.deg, kmtnet_fields)
    print('Event not in OMEGA-II: ',event_status)