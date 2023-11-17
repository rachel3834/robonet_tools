from pyDANDIA import vizier_tools, catalog_utils
from os import path
import argparse
from astropy.io import fits
from astropy.table import Table, Column
import numpy as np

def parse_gaia_vizier_catalog(args):

    # Load the Vizier catalog for the format of the survey selected:
    if args.cat_code == 'Gaia-DR2':
        catalog = load_vizier_gaiadr2(args)

    # Parse the Vizier-default format columns into those expected by
    # the pipeline:
    new_catalog = parse_catalog(args, catalog)

    # Output the revised catalog:
    catalog_utils.output_vizier_catalog(args.output_catalog, new_catalog,
                                        args.cat_code)

def parse_catalog(args, catalog):

    # Get the expected column definitions used by the pipeline
    supported_catalogs = vizier_tools.get_supported_catalogs()

    # Parse the table columns.  Since not all of the default columns are used
    # we can't just rename the column headers here.
    column_dict = supported_catalogs[args.cat_code][1]
    column_list = []
    column_order = column_dict.keys()
    for name1 in column_order:
        name2 = column_dict[name1]
        column_list.append(catalog[name1])

    new_catalog = Table(column_list)
    for name1, name2 in column_dict.items():
        new_catalog.rename_column(name1, name2)

    return new_catalog

def load_vizier_gaiadr2(args):

    if not path.isfile(args.input_catalog):
        raise IOError('Cannot find the input catalog downloaded from Vizier at '
                +args.input_catalog)

    with fits.open(args.input_catalog) as hdul:
        table_data = {}
        column_order = []
        for key, value in hdul[1].header.items():
            if 'TTYPE' in key:
                table_data[value] = []
                column_order.append(value)
        for row in hdul[1].data:
            for i,colname in enumerate(column_order):
                table_data[colname].append(row[i])

    column_list = []
    for colname in column_order:
        column_list.append(Column(name=colname, data=np.array(table_data[colname])))

    catalog = Table(column_list)

    return catalog

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_catalog', type=str,
                    help='Path to the input catalog downloaded from Vizier:')
    parser.add_argument("cat_code", type=str,
                    help='Which catalog format?  Currently supported: Gaia-DR2')
    parser.add_argument("output_catalog", type=str,
                    help='Path to the output catalog')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    parse_gaia_vizier_catalog(args)
