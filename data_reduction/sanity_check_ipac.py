from os import path, walk
import argparse
import crossmatch_to_ipactable

def check_data_products(args):

    # Count the number of lightcurve files output by walking through the output directory
    lc_files = []
    for (root, dirs, files) in walk(args.output_dir):
        for file in files:
            if '_star_' in file:
                lc_files.append(path.join(root, file))
    print('Found a total of ' + str(len(lc_files)) + ' lightcurve data files')

    # Load the source catalog to get the expected number of stars
    source_catalog_file = path.join(args.output_dir, args.field_id + '_source_table.tbl')
    source_catalog = crossmatch_to_ipactable.read_ipactable(source_catalog_file)
    print('Loaded source catalog of ' + str(len(source_catalog)) + ' stars')

    # Load the original starcount files to verify the number of lightcurves expected
    starcounts = crossmatch_to_ipactable.load_starcounts(args)
    print(starcounts)

    # Verify that there is a lightcurve file for all catalog entries with a non-zero ndata
    # in at least one passband
    print('Verifying presence of lightcurve files...')
    for star in source_catalog:
        lc_path = path.join(args.output_dir, star[-1])
        if not path.isfile(lc_file):
            raise IOError('Cannot find expected lightcurve file ' + lc_path)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir', help='Path to top-level data directory')
    parser.add_argument('field_id', help='Name of field')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    check_data_products(args)
