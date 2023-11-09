from os import path
import argparse
import matplotlib.pyplot as plt
import numpy as np
import plot_phase_curve

def plot_lightcurve(args):
    dt = 2450000.0

    data_gp = plot_phase_curve.load_data(args.input_file_gp)
    data_rp = plot_phase_curve.load_data(args.input_file_rp)
    data_ip = plot_phase_curve.load_data(args.input_file_ip)

    fig = plt.figure(1, (10, 10))
    plt.rcParams.update({'font.size': 18})

    plt.errorbar(data_gp[:,0]-dt, data_gp[:,1], yerr=data_gp[:,2],
                 marker='v', fmt='mv', label='g-band')
    plt.errorbar(data_rp[:,0]-dt, data_rp[:,1], yerr=data_rp[:,2],
                 fmt='bs', label='r-band')
    plt.errorbar(data_ip[:,0]-dt, data_ip[:,1], yerr=data_ip[:,2],
                 fmt='r*', label='i-band')

    plt.xlabel('HJD - '+str(dt))
    plt.ylabel('Magnitude')
    plt.legend()
    plt.grid('top-right')

    [xmin,xmax,ymin,ymax] = plt.axis()
    plt.axis([xmin,xmax,ymax,ymin])
    plt.title(args.title)

    plt.savefig(args.output_file, bbox_inches='tight')

    plt.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file_gp', help='Path to phased lightcurve file in g-band')
    parser.add_argument('input_file_rp', help='Path to phased lightcurve file in r-band')
    parser.add_argument('input_file_ip', help='Path to phased lightcurve file in i-band')
    parser.add_argument('title', help="String for title (no spaces)")
    parser.add_argument('output_file', help='Path to output plot file')
    args = parser.parse_args()

    plot_lightcurve(args)