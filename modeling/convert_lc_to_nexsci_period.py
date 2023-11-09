from os import path
import argparse

def convert_lc_file(args):

    if not path.isfile(args.input_file):
        raise IOError('Cannot find input file ' + args.input_file)

    line_list = open(args.input_file, 'r').readlines()

    output_file = args.input_file.replace('.dat', '.txt')

    with open(output_file, 'w') as f:
        for line in line_list:
            if '#' not in line:
                entries = line.replace('\n','').split()
                if '0.0' in entries[-1]:
                    f.write(str(entries[0]) + ', ' + str(entries[1]) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to pyDANDIA .dat lightcurve file')
    args = parser.parse_args()

    convert_lc_file(args)