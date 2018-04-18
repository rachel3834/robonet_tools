import os
import sys
import glob

def verify_downloaded_data():
    """Function to verify that a WGET download of data from the LCO archive has properly obtained all frames.

    This function compares the TXT format table of frames produced from the archive with the list of images
    within the working directory.
    """

    data_dir = get_args()

    expected_images = read_archive_table_export(data_dir)

    downloaded_images = list_downloaded_images(data_dir)

    check_for_missing_images(expected_images, downloaded_images)


def get_args():
    """Function to request or gather the required arguments"""

    if len(sys.argv) == 1:

        data_dir = raw_input('Please enter the path to the working data directory: ')

    else:

        data_dir = sys.argv[1]

    return data_dir

def read_archive_table_export(data_dir):
    """Function to read the table of images produced by the archive.
    This table file should be produced using the TXT formatting option.
    """

    table_file = os.path.join(data_dir,'tableExport.txt')

    if not os.path.isfile(table_file):

        print('ERROR: Cannot find the archive-generated list of expected images')
        sys.exit()

    f = open(table_file,'r')
    flines = f.readlines()
    f.close()

    image_list = []

    for line in flines[1:]:

        entries = line.replace('\n','').split(',')

        image_list.append(entries[2])

    print(str(len(image_list))+' image(s) expected from the archive')

    return image_list

def list_downloaded_images(data_dir):
    """Function to build a list of image files that were actually downloaded from those available in the
    working directory, stripping the file extensions since these are not included in the archives data list."""

    file_list = glob.glob(os.path.join(data_dir,'*.fits.fz'))

    if len(file_list) == 0:

        file_list = glob.glob(os.path.join(data_dir,'*.fits'))

    if len(file_list) == 0:

        print('ERROR: No frames downloaded')
        sys.exit()

    image_list = []

    for f in file_list:

        image_list.append(os.path.basename(f).split('.')[0])

    print(str(len(image_list))+' image(s) downloaded')

    return image_list

def check_for_missing_images(expected_images, downloaded_images):
    """Function to compare the expected and downloaded images and flag a warning if some frames are missing."""

    n_missing = 0

    for eimage in expected_images:

        if eimage not in downloaded_images:

            print('MISSING IMAGE: '+eimage)

            n_missing += 1

    print(str(n_missing)+' image(s) missing from the downloaded dataset')


if __name__ == '__main__':

    verify_downloaded_data()

