from os import path
import json
import argparse

def count_variables_in_field(file_path):

    with open(file_path, 'r') as f:
        data = json.load(f)
        print(len(data))

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('events_file', help='Path to the events JSON file')
    parser.add_argument('variables_file', help='Path to the variables JSON file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    count_variables_in_field(args.events_file)
    count_variables_in_field(args.variables_file)