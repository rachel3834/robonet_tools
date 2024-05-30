import argparse
import glob
from os import path, rmdir
from shutil import rmtree

def review_data(args):
    """
    Function to review the data within a top-level directory, and remove older data where possible
    """

    # Identify reduction directories matching the search string, reviewing the list to eliminate
    # any stray files, since we are looking for reduction directories only
    dir_list = glob.glob(path.join(args.top_dir, args.search_key))
    data_list = []
    for dir in dir_list:
        if path.isdir(dir) and path.isdir(path.join(dir, 'data')):
            data_list.append(dir)

    # Walk over all of the directories in the data_list, removing the raw image data
    for dir in data_list:
        data_dir = path.join(dir, 'data')
        rmtree(data_dir)
        print('Removed the raw image data from ' + dir)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('top_dir', help='Path to top-level data directory')
    parser.add_argument('search_key', help="""
        Search string in double-quotes (wildcards alllowed) to identify reduction directories to review
        """)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    review_data(args)