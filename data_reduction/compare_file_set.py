from os import path
import filecmp
import argparse
import glob

def compare_files(args):

    # Make a list of the files in dir1 which should be present in dir2
    files = glob.glob(path.join(args.dir1, '*'))
    common = [ path.basename(x) for x in files]

    # Run comparison of the files in both directories:
    (match, mismatch, errors) = filecmp.cmpfiles(args.dir1, args.dir2, common)
    print('Matching files: ', match)
    print('Mismatching files: ', mismatch)
    print('Errors: ', errors)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('dir1', help='Path to the first data directory')
    parser.add_argument('dir2', help='Path to the second data directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    compare_files(args)
