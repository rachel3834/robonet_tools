from os import path
import phase_fold_lc
import argparse
import numpy as np

def filter_lc(args):

    data = phase_fold_lc.load_data(args)

    idx = np.where(data[:,3] == 0.0)[0]

    out_file = args.input_file.replace('.dat', '_cleaned.dat')
    with open(out_file, 'w') as f:
        f.write('# Photometry type: normalized\n')
        f.write('# HJD    mag       mag_error     QC_code\n')
        for i in range(0,len(idx),1):
            f.write(str(data[i,0]) + ' ' + str(data[i,1]) + ' ' + str(data[i,2]) + ' ' + str(data[i,3]) + '\n')
        f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to pyDANDIA .dat lightcurve file')
    args = parser.parse_args()

    filter_lc(args)