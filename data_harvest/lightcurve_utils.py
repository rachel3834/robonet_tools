from os import path
from sys import argv
import numpy as np
from astropy.table import Table, Column

def read_moa_lightcurve(input_file,ZP=27.4, flux_0=0.0):
    """Function to read a MOA lightcurve file"""

    filter = 'R'

    if path.isfile(input_file):

        file_lines = open(input_file, 'r').readlines()

        data = []
        for line in file_lines:
            if '#' not in line[0:1]:
                entries = line.replace('\n','').split()
                if len(entries) > 0:
                    (mag, merr) = moa_flux_to_mag(float(entries[1]), float(entries[2]), ZP=ZP, flux_0=flux_0)
                    if np.isnan(mag) == False and np.isnan(merr) == False:
                        data.append( [float(entries[0]), filter, mag, merr])
        data = np.array(data)

        lc = Table( [Column(name='time', data=data[:,0]),
                     Column(name='filter', data=data[:,1]),
                     Column(name='magnitude', data=data[:,2]),
                     Column(name='error', data=data[:,2])] )

    else:

        raise IOError('Cannot find lightcurve input file '+input_file)

    return lc

def moa_flux_to_mag(flux, flux_err, ZP=27.4, flux_0=0.0):

    if flux > 0.0 and flux_err > 0.0:
        mag = ZP - 2.5*np.log10(flux + flux_0)
        mag_err = (2.5 / np.log(10.0)) * flux_err / flux
    else:
        mag = np.nan
        mag_err = np.nan

    return mag, mag_err

def read_ASASSN_lightcurve(input_file):
    """Function to read a lightcurve in the standard ASAS-SN format"""

    filter = path.basename(input_file).split('.')[0].split('_')[-1]

    if path.isfile(input_file):

        file_lines = open(input_file, 'r').readlines()

        data = []
        for line in file_lines:
            if '#' not in line[0:1]:
                entries = line.replace('\n','').split()
                data.append( [float(entries[0]), filter, entries[1], entries[2]])

        data = np.array(data)

        lc = Table( [Column(name='time', data=data[:,0]),
                     Column(name='filter', data=data[:,1]),
                     Column(name='magnitude', data=data[:,2]),
                     Column(name='error', data=data[:,3])] )

    else:

        raise IOError('Cannot find lightcurve input file '+input_file)

    return lc

def output_tom_csv(lc, output_file):
    """Function to output a lightcurve in TOM-standard CSV format"""

    f = open(output_file, 'w')
    f.write('time,filter,magnitude,error\n')
    for i in range(0,len(lc),1):
        f.write(str(lc['time'][i])+','+str(lc['filter'][i])+','+\
                str(lc['magnitude'][i])+','+str(lc['error'][i])+'\n')
    f.close()

    print('Output lightcurve to '+output_file)

if __name__ == '__main__':

    supported_formats = {'MOA': read_moa_lightcurve,
                         'ASASSN': read_ASASSN_lightcurve}

    if len(argv) < 3:
        input_file = input('Please enter the path in an input file: ')
        format = input('Please indicate the format of the input file.  Options are: '\
                        +','.join(supported_formats.keys()))
    else:
        input_file = argv[1]
        format = argv[2]

    extn = path.basename(input_file).split('.')[-1]
    output_file = input_file.replace('.'+extn, '.csv')

    if format in supported_formats.keys():

        lc = supported_formats[format](input_file)
        output_tom_csv(lc, output_file)

    else:
        raise IOError('Unrecognized lightcurve format')

    print(lc)
