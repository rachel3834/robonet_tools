from os import path
import argparse
import numpy as np

def phase_fold(args):

    data = load_data(args)

    phase = ((data[:,0] - float(args.epoch)) / float(args.period)) % 1

    output_phase_data(args, phase, data)

def output_phase_data(args, phase, data):

    output_file = args.input_file.replace('.dat', '_ph.dat')

    with open(output_file, 'w') as f:
        f.write('# Phase   JD [days]   Mag    Mag_error\n')
        for i in range(0,len(data),1):
            f.write(str(phase[i]) + ' ' + str(data[i,0]) + ' ' + str(data[i,1]) + ' ' + str(data[i,2]) + '\n')

def load_data(args):
    if not path.isfile(args.input_file):
        raise IOError('Cannot find input file ' + args.input_file)

    line_list = open(args.input_file, 'r').readlines()
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
    parser.add_argument('input_file', help='Path to pyDANDIA .dat lightcurve file')
    parser.add_argument('period', help='Period in days')
    parser.add_argument('epoch', help='Epoch in JD [days]')
    args = parser.parse_args()

    phase_fold(args)