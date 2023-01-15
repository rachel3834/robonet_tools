from os import path
from sys import argv
import glob

def count_event_datasets(top_dir):
    """The purpose of this function is to review the datasets gathered on
    all microlensing observations and count the number of events observed"""

    dir_list = glob.glob(path.join(top_dir, '*'))

    datasets = {}
    for item in dir_list:
        if '.txt' not in item:
            entries = str(path.basename(item)).split('-')
            target = entries[0]
            dataid = str(path.basename(item)).replace(target+'_','')
            data_dir = path.join(item, 'data')
            if path.isdir(data_dir):
                data_list = glob.glob(path.join(item, 'data', '*.fits'))

                if target in datasets.keys():
                    data = datasets[target]
                    data[dataid] = len(data_list)
                else:
                    data = {dataid: len(data_list)}
                datasets[target] = data

    print('Number of targets: '+str(len(datasets)))
    print(datasets.keys())

if __name__ == '__main__':
    if len(argv) > 1:
        top_dir = argv[1]
    else:
        top_dir = input('Please enter the path to the top level directory: ')
    count_event_datasets(top_dir)
