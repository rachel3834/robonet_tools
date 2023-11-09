from os import path
import argparse
import matplotlib.pyplot as plt
import numpy as np

def plot_phased_lightcurve(args):

    data_gp = load_data(args.input_file_gp)
    data_rp = load_data(args.input_file_rp)
    data_ip = load_data(args.input_file_ip)

    fig = plt.figure(1, (10, 10))
    plt.rcParams.update({'font.size': 18})

    plt.errorbar(data_gp[:,0], data_gp[:,2], yerr=data_gp[:,3],
                 marker='v', fmt='mv', label='g-band')
    plt.errorbar(data_rp[:,0], data_rp[:,2], yerr=data_rp[:,3],
                 fmt='bs', label='r-band')
    plt.errorbar(data_ip[:,0], data_ip[:,2], yerr=data_ip[:,3],
                 fmt='r*', label='i-band')

    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    plt.legend()
    plt.grid('top-right')

    [xmin,xmax,ymin,ymax] = plt.axis()
    plt.axis([xmin,1.3,ymax,ymin])
    plt.title(args.title)

    plt.savefig(args.output_file, bbox_inches='tight')

    plt.close()

def load_data(file_path):
    if not path.isfile(file_path):
        raise IOError('Cannot find input file ' + file_path)

    line_list = open(file_path, 'r').readlines()
    data = []

    for line in line_list:
        if '#' not in line:
            entries = line.replace('\n','').split()
            if '0.0' in entries[-1]:
                data.append([float(x) for x in entries])

    data = np.array(data)

    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file_gp', help='Path to phased lightcurve file in g-band')
    parser.add_argument('input_file_rp', help='Path to phased lightcurve file in r-band')
    parser.add_argument('input_file_ip', help='Path to phased lightcurve file in i-band')
    parser.add_argument('output_file', help='Path to output plot file')
    parser.add_argument('title', help="Title")
    args = parser.parse_args()

    plot_phased_lightcurve(args)