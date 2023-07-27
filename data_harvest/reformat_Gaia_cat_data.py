from os import path
from sys import argv
from astropy.io import votable
from astropy.io import fits
from astropy.table import Table, Column

def repackage_gaia_catalog():
    params = get_args()

    table_data = load_votable(params)
    output_fitstable(params, table_data)

def get_args():

    params = {}
    if len(argv) == 1:
        params['votable_path'] = input('Please enter the path to the VOTable file: ')
        params['fits_path'] = input('Please enter the path to the output FITS table: ')
    else:
        params['votable_path'] = argv[1]
        params['fits_path'] = argv[2]

    return params

def load_votable(params):

    if not path.isfile(params['votable_path']):
        raise IOError('Cannot find input VOTable at '+params['votable_path'])

    vo_object = votable.parse(params['votable_path'])
    table_data = vo_object.get_first_table().to_table(use_names_over_ids=True)

    return table_data

def output_fitstable(params, table_data):

    column_list = []
    for col_name in table_data.colnames:
        f = table_data[col_name].format
        if f == None:
            f = 'S200'
        col = Column(name=col_name, data=table_data[col_name], dtype=f)
        column_list.append(col)
    table = fits.BinTableHDU.from_columns(column_list)

    hdu = fits.PrimaryHDU(table)
    hdul = fits.HDUList([hdu])
    hdul.write(params['fits_path'], format='fits', overwrite=True)

if __name__ == '__main__':
    repackage_gaia_catalog()
