from pyDANDIA import crossmatch
import argparse
import numpy as np
from datetime import datetime
from os import path

VERSION = 0.1

def convert_to_ipactable(args):
    """Function to convert a pyDANDIA CrossMatchTable for a single field to
    an ASCII table in IPAC format, based on the documentation provided at
    https://exoplanetarchive.ipac.caltech.edu/docs/ddgen/ipac_tbl.html
    """

    # Load the field's crossmatch table
    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.crossmatch_file,log=None)

    # Data on each source is derived from a combination of two tables in the
    # crossmatch file, the field index and the stars table, so we extract that
    # first
    source_table = extract_source_data(xmatch)

    # Output the source table in IPAC table format
    output_to_ipactable(args, source_table)

def get_column_widths(table_columns):
    """Function to set the width of a column header, based on the length of the
    column header string"""

    for col in table_columns.keys():
        col_def = table_columns[col]
        if 'width' not in col_def.keys():
            # Width of the column is defined to be the maximum of
            # (length of the column name, length of the type name, length of the units name)
            # plus two spaces before and after, not including the leading and trailing |
            width = max(len(col), len(col_def['type']), len(col_def['unit'])) + 4
            col_def['width'] = width
        table_columns[col] = col_def

    return table_columns

def padd_header_entry(entry, width):

    while len(entry) < width:
        if len(entry)%2 == 0:
            entry = ' '+entry
        else:
            entry = entry+' '

    return entry

def padd_column_entry(entry, width):

    while len(entry) < width:
        entry = entry+' '

    return entry

def output_to_ipactable(args, source_table):
    """Function to output the source table to the IPAC table format"""

    # Formal definition of the columns in the source table:
    table_columns = {
                    'field_id': {'type': 'int', 'unit': 'null', 'nulls': 'null'},
                    'ra': {'type': 'double', 'unit': 'degrees', 'nulls': 'null'},
                    'dec': {'type': 'double', 'unit': 'degrees', 'nulls': 'null'},
                    'quadrant': {'type': 'int', 'unit': 'null', 'nulls': 'null'},
                    'quadrant_id': {'type': 'int', 'unit': 'null', 'nulls': 'null'},
                    'gaia_source_id': {'type': 'int', 'unit': 'null', 'nulls': 'null', 'width': 21},
                    'cal_mag_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'cal_mag_error_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'norm_mag_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'norm_mag_error_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'cal_mag_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'cal_mag_error_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'norm_mag_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'norm_mag_error_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'cal_mag_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'cal_mag_error_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'norm_mag_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'norm_mag_error_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0'},
                    'lc_file_path': {'type': 'chat', 'unit': 'null', 'nulls': 'null', 'width': 300},
                    }

    # Calculate the widths for each column:
    table_columns = get_column_widths(table_columns)

    # Construct file header
    tbl_file = open(args.ipactable_file, 'w')
    tbl_file.write('\catalog=romerea\n')
    tbl_file.write('\catalog_version='+VERSION+'\n')
    now = datetime.utcnow()
    tbl_file.write('\date='+now.strftime("%Y-%m-%dT%H:%M:%D")+'\n')
    header = '|'
    units = '|'
    null_values = '|'
    for col, col_def in table_columns.items():
        header += padd_header_entry(col, col_def['width']) + '|'
        units += padd_header_entry(col_def['unit'], col_def['width']) + '|'
        null_values += padd_header_entry(col_def['nulls'], col_def['width']) + '|'
    tbl_file.write(header+'\n')
    tbl_file.write(units+'\n')
    tbl_file.write(null_values+'\n')

    # Write file data
    for star in source_table:
        tbl_line = ' '
        for col, col_def in table_columns.items():
            if col == 'lc_file_path':
                tbl_line += get_lc_file_path(args, star['field_id'])
            elif col in ['ra', 'dec']:
                value = round(star[col],5)
                tbl_line += padd_column_entry(str(value), col_def['width']) + ' '
            elif 'mag' in col:
                value = round(star[col],3)
                tbl_line += padd_column_entry(str(value), col_def['width']) + ' '
            else:
                tbl_line += padd_column_entry(str(star[col]), col_def['width']) + ' '
        tbl_file.write(tbl_line+'\n')

    tbl_file.close()

def get_lc_file_path(args, field_id):
    """Function to assign a file path for the lightcurve of a star, based on
    the field name and the star's field ID within that field.

    Star lightcurves are organized into directories by field, and then into
    sub-directories of 1000 lightcurves each based on the star's field ID, i.e.:
    top_level_directory
            |- ROME-FIELD-01
            |- ROME-FIELD-02
            |- <etc>
            |- ROME-FIELD-20
                  |- FieldID000000_0001000
                  |- FieldID001000_0002000
                  |- FieldID002000_0003000
                  |- FieldID003000_0004000
                  |- <etc>

    The filenames assigned to each star are constructed as follows:
    <field_name>_star_<field_ID>.fits
    """
    def bin_num_string(num):
        num = str(num)
        while len(num) < 6:
            num = '0' + num
        return num

    bin_width = 1000
    len_num = len(str(field_id))

    if len_num < 4:                     # 1 - 999
        subdir = 'FieldID000000_0001000'
    if len_num >= 4:    # 1000 - 9999
        nmin = int(field_id / bin_width) * bin_width
        nmax = nmin + bin_width
        subdir = 'FieldID' + bin_num_string(nmin) + '_' + bin_num_string(nmax)

    filename = args.field_name+'_star_'+str(field_id)+'.fits'

    lc_path = path.join('./', args.field_name, subdir, filename)

    return lc_path

def extract_source_data(xmatch):
    """Function to extract the data on all stars within the field.  This
    combines information held in the field index and stars tables of the crossmatch
    table.

    The crossmatch table records the calibrated and normalized photometry for each
    star in all passbands and all datasets separately.  However, most stars are
    not measured in all datasets.  Where a star has valid data from LCO Chile,
    Dome A, those calibrated and normalized photometry have been used.
    If data from Chile are not available, the photometry from South Africa and
    then Australia were used for normalization in that order, so those
    magnitudes are provided instead.
    """
    # Star field indices, positional and quadrant data come from the field index
    col_list_field_idx = [
                          'field_id', 'ra', 'dec', 'quadrant', 'quadrant_id',
                          'gaia_source_id'
                         ]
    source_table = xmatch.field_index[col_list_field_idx]
    nstars = len(source_table)

    # Create columns to hold the photometry data
    phot_cols = [
                 'cal_mag_g', 'cal_mag_error_g', 'norm_mag_g', 'norm_mag_error_g',
                 'cal_mag_r', 'cal_mag_error_r', 'norm_mag_r', 'norm_mag_error_r',
                 'cal_mag_i', 'cal_mag_error_i', 'norm_mag_i', 'norm_mag_error_i',
                ]
    for col in phot_cols:
        source_table.add_column(np.zeros(nstars), name=col)

    # ROME/REA photometry were normalized by selecting the instruments used for
    # the survey data as primary references.
    filter_list = ['gp', 'rp', 'ip']
    reference_list = [
                    'lsc-doma-1m0-05-fa15',
                    'cpt-doma-1m0-10-fa16',
                    'coj-doma-1m0-11-fa12'
                    ]
    for bandpass in filter_list:
        f = bandpass.replace('p','')

        # For stars measured directly in the primary reference data, we use the
        # photometry from the star's table for the relevant columns.
        for pri_ref_code in reference_list:
            site_code = get_site_dome(facility_code=pri_ref_code)
            idx = select_stars_from_ref_dataset(xmatch, reference_list, pri_ref_code,
                                                bandpass, args.field_name)
            source_table['cal_mag_'+f][idx] = xmatch.stars['cal_'+f+'_mag_'+site_code][idx]
            source_table['cal_mag_error_'+f][idx] = xmatch.stars['cal_'+f+'_magerr_'+site_code][idx]

            # The normalized mag columns are only populated if this value is
            # different from the calibrated magnitude for that reference image,
            # i.e. if the reference was normalized to a primary reference:
            for j in idx:
                if xmatch.stars['norm_'+f+'_mag_'+site_code][j] > 0.0:
                    source_table['norm_mag_'+f][j] = xmatch.stars['norm_'+f+'_mag_'+site_code][j]
                    source_table['norm_mag_error_'+f][j] = xmatch.stars['norm_'+f+'_magerr_'+site_code][j]
                else:
                    source_table['norm_mag_'+f][j] = xmatch.stars['cal_'+f+'_mag_'+site_code][j]
                    source_table['norm_mag_error_'+f][j] = xmatch.stars['cal_'+f+'_magerr_'+site_code][j]

        # Some stars were only measured in one of the REA follow-up datasets
        # and none of the primary reference datasets:
        jdx = select_stars_no_ref_dataset(xmatch, reference_list, bandpass,
                                            args.field_name)

        # For REA-only stars, we work through the list of datasets in turn,
        # taking the photometry from the first dataset in the list to
        # provide valid measurements.

        # Fetch an ordered list of all of the datasets in this bandpass
        ddx = np.where(xmatch.datasets['dataset_filter'] == bandpass)[0]
        datasets = xmatch.datasets['dataset_code'][ddx]

        for j in jdx:
            for dset in datasets:
                if dset not in reference_list:
                    if xmatch.field_index[dset+'_index'][j] > 0:
                        site_code = get_site_dome(dataset_code=dset)
                        source_table['cal_mag_'+f][j] = xmatch.stars['cal_'+f+'_mag_'+site_code][j]
                        source_table['cal_mag_error_'+f][j] = xmatch.stars['cal_'+f+'_magerr_'+site_code][j]
                        source_table['norm_mag_'+f][j] = xmatch.stars['norm_'+f+'_mag_'+site_code][j]
                        source_table['norm_mag_error_'+f][j] = xmatch.stars['norm_'+f+'_magerr_'+site_code][j]

    return source_table

def get_site_dome(facility_code=None, dataset_code=None):

    if dataset_code:
        facility_code = dataset_code.split('_')[1]

    entries = facility_code.split('-')
    site_dome = entries[0]+'_'+entries[1]

    return site_dome

def select_stars_from_ref_dataset(xmatch, reference_list, pri_ref_code,
                                        bandpass, field_name):
    """Function selects all stars in the field index with measurements in the
    primary reference dataset in the given filter.
    The reference_list gives the order of preference to use the photometry from
    the three datasets that can be used as a primary reference: LSC, CPT, COJ (Dome A).
    """

    dset_index = reference_list.index(pri_ref_code)

    colname = field_name+'_'+pri_ref_code+'_'+bandpass+'_index'
    idx = set(np.where(xmatch.field_index[colname] > 0)[0])

    # If the primary reference in use is the first reference, select all stars
    # measured in that dataset.  If the primary reference dataset in use is
    # subsequent in the list, select all stars that are measured in the primary
    # reference, but then exclude those that have valid measurements in the
    # preceeding references:
    for i in range(0,dset_index,1):
        colname2 = field_name+'_'+reference_list[i]+'_'+bandpass+'_index'
        idx2 = np.where(xmatch.field_index[colname2] == 0)[0]
        idx = idx.intersection(set(idx2))

    return list(idx)

def select_stars_no_ref_dataset(xmatch, reference_list, bandpass, field_name):
    """Some stars were never measured in any of the Dome A datasets from different
    sites that were used as the primary reference.
    for these cases, the only available reference frame photometry has been
    calibrated and normalized from stars that do overlap the primary reference
    images.  This function identifies these stars from the field index"""

    # Start with a list of all of the stars in the field
    stars = range(0,len(xmatch.field_index),1)
    idx = set(stars)

    # Then exclude stars that have measurements in each of the primary
    # reference datasets:
    for dset in reference_list:
        colname = field_name+'_'+dset+'_'+bandpass+'_index'
        idx2 = np.where(xmatch.field_index[colname] == 0)[0]
        idx = idx.intersection(set(idx2))

    return list(idx)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('crossmatch_file', type=str,
                    help='Path to crossmatch file')
    parser.add_argument('ipactable_file', type=str,
                    help='Path to output IPAC table file')
    parser.add_argument('field_name', type=str,
                    help='Name of the field')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    convert_to_ipactable(args)
