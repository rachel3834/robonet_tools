from pyDANDIA import crossmatch
from pyDANDIA import logs
import argparse
import numpy as np
from datetime import datetime
from os import path
from astropy.table import Table, Column
import json

VERSION = '0.1'

def convert_to_ipactable(args):
    """Function to convert a pyDANDIA CrossMatchTable for a single field to
    an ASCII table in IPAC format, based on the documentation provided at
    https://exoplanetarchive.ipac.caltech.edu/docs/ddgen/ipac_tbl.html
    """

    log = logs.start_stage_log(args.output_dir, 'crossmatch_to_ipactable')

    # Load the field's crossmatch table
    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.crossmatch_file,log=log)

    # Load the field's event and variable catalogs
    variable_catalog = load_json_target_catalog(args.variable_lut_file, log)

    # Load the starcount files for each quadrant, produced by the fieldphot_to_ipaclc code
    starcounts = load_starcounts(args, log)

    # Data on each source is derived from a combination of two tables in the
    # crossmatch file, the field index and the stars table, so we extract that
    # first
    source_table = extract_source_data(args, xmatch, variable_catalog, starcounts, log)

    # Output the source table in IPAC table format
    output_to_ipactable(args, source_table, log)

    log.info('End of processing')
    logs.close_log(log)

def load_starcounts(args, log=None):
    """
    Function to load the JSON files that record the number of valid datapoints obtained for each
    star.  These files are output by fieldphot_to_ipaclc separately for each quadrant, so
    this function also concatenates the data into a single dictionary.
    """

    starcounts = {}
    for qid in range(1,5,1):
        input_file = path.join(args.output_dir, '..', args.field_name + '_starcounts_Q' + str(qid) + '.json')

        if not path.isfile(input_file):
            raise IOError('Missing star count data for field quadrants.  Looking for ' + input_file)

        with open(input_file, "r") as read_file:
            quad_counts = json.load(read_file)
            read_file.close()

        for j, counts in quad_counts.items():
            starcounts[int(j)] = counts

    if log:
        log.info('Loaded data on number of datapoints in all star lightcurves')

    return starcounts

def load_json_target_catalog(catalog_path, log=None):
    """Function to load a catalog of selected objects in JSON format"""

    with open(catalog_path, "r") as read_file:
        data = json.load(read_file)
        read_file.close()

    if log:
        log.info('Loaded look-up table of known variables and events')

    return data

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

def define_table_columns():
    """
    Formal definition of the columns in the source table.  Widths are set by either the
    maximum width of the likely parameter value or the width of the column name, whichever
    is larger
    """

    table_columns = {
        'name': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 20},
        'field': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 13},
        'field_id': {'type': 'int', 'unit': 'null', 'nulls': 'null', 'width': 8},
        'ra': {'type': 'double', 'unit': 'degrees', 'nulls': 'null', 'width': 10},
        'dec': {'type': 'double', 'unit': 'degrees', 'nulls': 'null', 'width': 9},
        'quadrant': {'type': 'int', 'unit': 'null', 'nulls': 'null', 'width': 8},
        'quadrant_id': {'type': 'int', 'unit': 'null', 'nulls': 'null', 'width': 11},
        'gaia_source_id': {'type': 'int', 'unit': 'null', 'nulls': 'null', 'width': 21},
        'ogle_event_id': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 18},
        'ogle_variable_id': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 20},
        'moa_event_id': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 16},
        'kmtnet_event_id': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 17},
        'spitzer_event': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 13},
        'vvv_variable_id': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 15},
        'cal_mag_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 9},
        'cal_mag_error_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 15},
        'norm_mag_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 10},
        'norm_mag_error_g': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 16},
        'cal_mag_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 9},
        'cal_mag_error_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 15},
        'norm_mag_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 10},
        'norm_mag_error_r': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 16},
        'cal_mag_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 9},
        'cal_mag_error_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 15},
        'norm_mag_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 10},
        'norm_mag_error_i': {'type': 'float', 'unit': 'mag', 'nulls': '0.0', 'width': 16},
        'ndata_g': {'type': 'int', 'unit': 'null', 'nulls': '0', 'width': 7},
        'ndata_r': {'type': 'int', 'unit': 'null', 'nulls': '0', 'width': 7},
        'ndata_i': {'type': 'int', 'unit': 'null', 'nulls': '0', 'width': 7},
        'lc_file_path': {'type': 'char', 'unit': 'null', 'nulls': 'null', 'width': 300}
    }

    return table_columns

def table_column_names():

    column_names = [
        'name',
        'field',
        'field_id',
        'ra',
        'dec',
        'quadrant',
        'quadrant_id',
        'gaia_source_id',
        'ogle_event_id',
        'ogle_variable_id',
        'moa_event_id',
        'kmtnet_event_id',
        'spitzer_event',
        'vvv_variable_id',
        'cal_mag_g',
        'cal_mag_error_g',
        'norm_mag_g',
        'norm_mag_error_g',
        'cal_mag_r',
        'cal_mag_error_r',
        'norm_mag_r',
        'norm_mag_error_r',
        'cal_mag_i',
        'cal_mag_error_i',
        'norm_mag_i',
        'norm_mag_error_i',
        'ndata_g',
        'ndata_r',
        'ndata_i',
        'lc_file_path'
    ]

    return column_names

def output_to_ipactable(args, source_table, log):
    """Function to output the source table to the IPAC table format"""

    # Formal definition of the columns in the source table.  Widths are set by either the
    # maximum width of the likely parameter value or the width of the column name, whichever
    # is larger

    log.info('Outputting source catalog')

    table_columns = define_table_columns()
    #table_columns = get_column_widths(table_columns)

    source_catalog_file = path.join(args.output_dir, args.field_name + '_source_table.tbl')

    # Construct file header
    tbl_file = open(source_catalog_file, 'w')
    tbl_file.write('\\catalog=romerea\n')
    tbl_file.write('\\catalog_version='+VERSION+'\n')
    now = datetime.utcnow()
    tbl_file.write('\\date='+now.strftime("%Y-%m-%dT%H:%M:%D")+'\n')
    header = '|'
    dtypes = '|'
    units = '|'
    null_values = '|'
    for col, col_def in table_columns.items():
        header += padd_header_entry(col, col_def['width']) + '|'
        dtypes += padd_header_entry(col_def['type'], col_def['width']) + '|'
        units += padd_header_entry(col_def['unit'], col_def['width']) + '|'
        null_values += padd_header_entry(col_def['nulls'], col_def['width']) + '|'
    tbl_file.write(header+'\n')
    tbl_file.write(dtypes+'\n')
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

def read_ipactable(file_path):
    """
    Function to read the source catalog in IPAC table format
    """

    if not path.isfile(file_path):
        raise IOError('Cannot find IPAC source catalog ' + file_path)

    file_lines = open(file_path, 'r').readlines()

    column_names = table_column_names()
    str_col_names = ['name', 'field', 'gaia_source_id', 'ogle_event_id', 'ogle_variable_id',
                     'moa_event_id', 'kmtnet_event_id', 'spitzer_event', 'vvv_variable_id', 'lc_file_path']
    int_col_names = ['field_id', 'quadrant', 'quadrant_id', 'ndata_g', 'ndata_r', 'ndata_i']
    source_catalog_data = {col: [] for col in column_names}
    for line in file_lines:
        if '\\' not in line[0:1] and '|' not in line[0:1]:
            entries = line.replace('\n','').split()
            for i,col in enumerate(column_names):
                print(col, entries[i])
                if col in str_col_names:
                    source_catalog_data[col].append(entries[i])
                elif col in int_col_names:
                    source_catalog_data[col].append(int(entries[i]))
                else:
                    source_catalog_data[col].append(float(entries[i]))

    columns = [Column(name=col, data=source_catalog_data[col]) for col in column_names]
    source_catalog = Table(columns)
    print(source_catalog.colnames)

    return source_catalog

def bin_num_string(num):
    num = str(num)
    while len(num) < 6:
        num = '0' + num
    return num

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

    bin_width = 1000
    len_num = len(str(field_id))

    if len_num < 4:                     # 1 - 999
        subdir = 'FieldID000000_001000'
    if len_num >= 4:    # 1000 - 9999
        nmin = int(field_id / bin_width) * bin_width
        nmax = nmin + bin_width
        subdir = 'FieldID' + bin_num_string(nmin) + '_' + bin_num_string(nmax)

    filename = args.field_name+'_star_'+str(bin_num_string(field_id))+'.fits'

    lc_path = path.join('./', args.field_name, subdir, filename)

    return lc_path

def extract_source_data(args, xmatch, variable_catalog, starcounts, log):
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

    log.info('Compiling data on all stars')

    # Sanity check that the starcounts dictionary has entries for all stars in the field index
    if len(starcounts) != len(xmatch.field_index['field_id']):
        raise IOError('Starcounts dictionary mismatch with field index')

    # Create an index of stars with data in at least one passband
    # There are legitimate reasons why a given star may be in the source catalog (=detected in
    # at least one reference image) but produce no timeseries data in one or more passbands
    # (=insufficient signal in subsequent images).  We want to output lightcurve files
    # for stars with timeseries data in one or more bands, but exclude stars with
    # no timeseries measurements in any band.
    invalid_star_idx = [field_id-1 for field_id, entry in starcounts.items()
                        if (entry['gp'] == 0 and entry['rp'] == 0 and entry['ip'] == 0)]
    log.info('Identifed ' + str(len(invalid_star_idx)) + ' stars with no timeseries data in any bandpass:')
    log.info(repr(invalid_star_idx))

    # TABLE CREATION
    # In order to keep the indexing consistent, the source table is created with the full
    # list of stars in the field.  Invalid stars are eliminated later.
    column_list = []
    table_columns = define_table_columns()

    # Create columns to hold the star and field identifiers
    nstars = len(xmatch.field_index['field_id'])
    star_ids = []
    for j in xmatch.field_index['field_id']:
        star_ids.append(args.field_name + '_' + zero_padd(j, 6))
    column_list.append(Column(name='name', data=star_ids))
    column_list.append(Column(name='field', data=[args.field_name]*nstars))

    # Star field indices, positional and quadrant data come from the field index
    col_list_field_idx = [
                          'field_id', 'ra', 'dec', 'quadrant', 'quadrant_id',
                          'gaia_source_id'
                         ]
    for col in col_list_field_idx:
        column_list.append(xmatch.field_index[col])

    # Include columns for crossmatches with other catalogs.  Default entry is the maximum width
    # of the column when empty
    cols = {
        'ogle_event_id': np.array([' '*table_columns['ogle_event_id']['width']]*nstars),
        'ogle_variable_id': np.array([' '*table_columns['ogle_variable_id']['width']]*nstars),
        'moa_event_id': np.array([' '*table_columns['moa_event_id']['width']]*nstars),
        'kmtnet_event_id': np.array([' '*table_columns['kmtnet_event_id']['width']]*nstars),
        'spitzer_event': np.array(['false']*nstars),
        'vvv_variable_id': np.array([' '*table_columns['vvv_variable_id']['width']]*nstars)
    }
    for field_id, entry in variable_catalog.items():
        field_idx = int(field_id) - 1
        for col_name, arr in cols.items():
            arr[field_idx] = entry[col_name]

    for col, arr in cols.items():
        column_list.append(Column(name=col, data=arr))

    # Create columns to hold the photometry data
    phot_cols = [
                 'cal_mag_g', 'cal_mag_error_g', 'norm_mag_g', 'norm_mag_error_g',
                 'cal_mag_r', 'cal_mag_error_r', 'norm_mag_r', 'norm_mag_error_r',
                 'cal_mag_i', 'cal_mag_error_i', 'norm_mag_i', 'norm_mag_error_i',
                ]
    for col in phot_cols:
        column_list.append(Column(name=col, data=np.zeros(nstars)))

    # Create columns to hold the number of datapoints
    ndata_g = np.zeros(nstars, dtype=int)
    ndata_r = np.zeros(nstars, dtype=int)
    ndata_i = np.zeros(nstars, dtype=int)
    for field_id, entry in starcounts.items():
        field_idx = int(field_id) - 1
        ndata_g[field_idx] = int(entry['gp'])
        ndata_r[field_idx] = int(entry['rp'])
        ndata_i[field_idx] = int(entry['ip'])
    column_list.append(Column(name='ndata_g', data=ndata_g))
    column_list.append(Column(name='ndata_r', data=ndata_r))
    column_list.append(Column(name='ndata_i', data=ndata_i))

    source_table = Table(column_list)
    log.info('Built source catalog table')

    # Filter the source_table for entries where no Gaia source ID is available,
    # and ensure that the entry is given as null not None
    jdx = source_table['gaia_source_id'] == 'None'
    source_table['gaia_source_id'][jdx] = 'null'

    # DATA POPULATION
    # ROME/REA photometry were normalized by selecting the instruments used for
    # the survey data as primary references.
    filter_list = ['gp', 'rp', 'ip']
    #reference_list = [
    #                'lsc-doma-1m0-05-fa15',
    #                'cpt-doma-1m0-10-fa16',
    #                'coj-doma-1m0-11-fa12'
    #                ]
    survey_reference = 'lsc-doma-1m0-05-fa15'
    for bandpass in filter_list:
        f = bandpass.replace('p','')

        # The cal_mag and norm_mag columns in the source catalog correspond to those
        # data from the survey's primary reference dataset, i.e. lsc-doma-fa15.
        # Where stars were not measured in these data, these columns may be zero;
        # this means they were detected only on other cameras.
        # Only the data from lsc-doma are included here so that all data provided have
        # a consistent calibration and are inter-comparable.
        source_table['cal_mag_' + f] = xmatch.stars['cal_' + f + '_mag_lsc_doma']
        source_table['cal_mag_error_' + f] = xmatch.stars['cal_' + f + '_magerr_lsc_doma']

        # If the stars were detected in the primary reference, they have no normalization mags
        # because they don't need to be normalized - they ARE the primary reference. So the norm_mag
        # columns are populated from the cal_mag columns.
        pri_ref_col = args.field_name + '_' + survey_reference + '_' + bandpass + '_index'
        jdx = np.where(xmatch.field_index[pri_ref_col] > 0)[0]
        source_table['norm_mag_' + f][jdx] = xmatch.stars['cal_' + f + '_mag_lsc_doma'][jdx]
        source_table['norm_mag_error_' + f][jdx] = xmatch.stars['cal_' + f + '_magerr_lsc_doma'][jdx]

        # For stars measured directly in the primary reference data, we use the
        # photometry from the star's table for the relevant columns.
        # for pri_ref_code in reference_list:
        #     site_code = get_site_dome(facility_code=pri_ref_code)
        #     idx = select_stars_from_ref_dataset(xmatch, reference_list, pri_ref_code,
        #                                         bandpass, args.field_name)
        #     source_table['cal_mag_'+f][idx] = xmatch.stars['cal_'+f+'_mag_'+site_code][idx]
        #     source_table['cal_mag_error_'+f][idx] = xmatch.stars['cal_'+f+'_magerr_'+site_code][idx]
        #
        #     # The normalized mag columns are only populated if this value is
        #     # different from the calibrated magnitude for that reference image,
        #     # i.e. if the reference was normalized to a primary reference:
        #     for j in idx:
        #         if xmatch.stars['norm_'+f+'_mag_'+site_code][j] > 0.0:
        #             source_table['norm_mag_'+f][j] = xmatch.stars['norm_'+f+'_mag_'+site_code][j]
        #             source_table['norm_mag_error_'+f][j] = xmatch.stars['norm_'+f+'_magerr_'+site_code][j]
        #         else:
        #             source_table['norm_mag_'+f][j] = xmatch.stars['cal_'+f+'_mag_'+site_code][j]
        #             source_table['norm_mag_error_'+f][j] = xmatch.stars['cal_'+f+'_magerr_'+site_code][j]
        #
        # # Some stars were only measured in one of the REA follow-up datasets
        # # and none of the primary reference datasets:
        # jdx = select_stars_no_ref_dataset(xmatch, reference_list, bandpass,
        #                                     args.field_name)
        #
        # # For REA-only stars, we work through the list of datasets in turn,
        # # taking the photometry from the first dataset in the list to
        # # provide valid measurements.
        #
        # # Fetch an ordered list of all of the datasets in this bandpass
        # ddx = np.where(xmatch.datasets['dataset_filter'] == bandpass)[0]
        # datasets = xmatch.datasets['dataset_code'][ddx]
        #
        # for j in jdx:
        #     for dset in datasets:
        #         if dset not in reference_list:
        #             if xmatch.field_index[dset+'_index'][j] > 0:
        #                 site_code = get_site_dome(dataset_code=dset)
        #                 source_table['cal_mag_'+f][j] = xmatch.stars['cal_'+f+'_mag_'+site_code][j]
        #                 source_table['cal_mag_error_'+f][j] = xmatch.stars['cal_'+f+'_magerr_'+site_code][j]
        #                 source_table['norm_mag_'+f][j] = xmatch.stars['norm_'+f+'_mag_'+site_code][j]
        #                 source_table['norm_mag_error_'+f][j] = xmatch.stars['norm_'+f+'_magerr_'+site_code][j]
    log.info('Populated source catalog table')

    # Now we filter the source catalog table, removing stars with no valid data:
    source_table.remove_rows(invalid_star_idx)

    return source_table

def zero_padd(ivalue, nsf):

    svalue = str(ivalue)
    nzeros = nsf - len(svalue)
    prefix = ''.join(['0']*nzeros)
    svalue = prefix + svalue

    return svalue

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
    parser.add_argument('variable_lut_file', type=str,
                        help='Path to look-up table of known events and variables')
    parser.add_argument('output_dir', type=str,
                    help='Path to output directory')
    parser.add_argument('field_name', type=str,
                    help='Name of the field')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    convert_to_ipactable(args)
