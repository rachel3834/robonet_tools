from os import path
import csv
import argparse
import numpy as np

def dat_to_csv(args):

    data = load_dat_lc(args)
    output_to_csv(args, data)

def output_to_csv(args, data):

    with open(args.output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['time','filter','magnitude','error'])
        for i in range(0,len(data),1):
            writer.writerow([str(data[i,0]),args.filter,str(data[i,1]),str(data[i,2])])

def load_dat_lc(args):
    if not path.isfile(args.input_file):
        raise IOError('Cannot find input file '+args.input_file)

    file_lines = open(args.input_file, 'r').readlines()
    data = []
    for line in file_lines[1:]:
        entries = line.replace('\n','').split()
        data.append( [float(entries[0]), float(entries[1]), float(entries[2])] )

    return np.array(data)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Path to the main JSON file of events", type=str)
    parser.add_argument("filter", help="Label for the filter used to take the data", type=str)
    parser.add_argument("output_file", help="Path to the table of Spitzer events", type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    dat_to_csv(args)
