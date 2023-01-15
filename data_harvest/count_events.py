from os import path
from sys import argv
import glob
import matplotlib.pyplot as plt
import numpy as np

def count_event_datasets(top_dir):
    """The purpose of this function is to review the datasets gathered on
    all microlensing observations and count the number of events observed"""

    # Number of images overwhich a dataset is considered to be substantial
    min_images = 10

    dir_list = glob.glob(path.join(top_dir, '*'))

    datasets = {}
    for item in dir_list:
        if '.txt' not in item:
            entries = str(path.basename(item)).split('_')
            target = entries[0]
            dataid = str(path.basename(item)).replace(target+'_','')
            data_dir = path.join(item, 'data')
            if path.isdir(data_dir):
                data_list = glob.glob(path.join(item, 'data', '*.fits'))
                if len(data_list) > min_images:
                    if target in datasets.keys():
                        data = datasets[target]
                        data[dataid] = len(data_list)
                    else:
                        data = {dataid: len(data_list)}
                    datasets[target] = data

    print('Number of targets: '+str(len(datasets)))

    nimages = []
    for target, data in datasets.items():
        for dataid, ndata in data.items():
            nimages.append(ndata)
    nimages = np.array(nimages)

    fig = plt.figure(1,(10,10))
    plt.hist(nimages)
    plt.xlabel('Number of images per dataset')
    plt.ylabel('Number of datasets')
    plt.title('OMEGA observations per dataset')

if __name__ == '__main__':
    if len(argv) > 1:
        top_dir = argv[1]
    else:
        top_dir = input('Please enter the path to the top level directory: ')
    count_event_datasets(top_dir)
