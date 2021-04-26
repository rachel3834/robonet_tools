import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors

parser = argparse.ArgumentParser()
parser.add_argument("event_name")
args = parser.parse_args()

event = args.event_name

directory = '/data/messier/omega/data/photometry/'+event

dataset = [i for i in os.listdir(directory) if 'dat' in i]

for data in dataset:
    filter_name = str(data).replace('.dat','').split('_')[-1].replace('p','')
    if filter_name == 'i':
        col = 'g'
        mkr = '*'
    else:
        col = 'm'
        mkr = 'd'
    aa = np.loadtxt(directory+'/'+data)
    plt.errorbar(aa[:,0]-2450000,aa[:,1],aa[:,2],
                    fmt='.',mfc=col, mec=col,label=filter_name)

plt.gca().invert_yaxis()
plt.grid()
plt.xlabel('HJD-2450000.0')
plt.ylabel('Instrumental magnitude')
plt.legend()
plt.title(event)
plt.savefig(directory+'/'+event+'_omega_lc.png')
