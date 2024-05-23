from os import path
import argparse
import glob

def build_transfer_filelist(args):
    """
    Function to build a complete list of files to transfer to IPAC via the FDT file upload system.
    Bulk uploads can be performed with an ASCII input list of all files.
    """

    data_list_file = path.join(args.data_dir, args.file_name)
    with open(data_list_file, 'w') as output_file:

        # First look for a source catalog in the top-level data directory.  There should only be
        # one (!) so raise an error if this isn't the case
        flist = glob.glob(path.join(args.data_dir, '*_source_table.tbl'))
        if len(flist) != 1:
            raise IOError('Missing or multiple source catalog files found in ' + args.data_dir)
        output_file.write(flist[0] + '\n')

        # Now look for a set of sub-directories starting with 'FieldID...'
        # These contain the FITS-format lightcurve files, each of which needs to be indexed
        dir_list = glob.glob(path.join(args.data_dir, 'FieldID*'))
        for subdir in dir_list:
            lc_list = glob.glob(path.join(subdir, '*.fits'))
            for fpath in lc_list:
                output_file.write(fpath + '\n')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Path to top level of data directory structure')
    parser.add_argument('field_name', help='Name of the field')
    parser.add_argument('file_name', help='Name of file for data list')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    build_transfer_filelist(args)