# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 13:24:43 2018

@author: rstreet
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import datetime

def analyze_obs_execution():
    """Function to analyze the execution of observations"""
    
    log_path = get_args()
    
    (frames,data) = read_obs_execution_log(log_path)
    
    plot_percent_completed(frames,data)

def plot_percent_completed(frames,data):
    """Function to generate a plot of the percentage of observation 
    subrequests that are executed as a function of date and instrument"""
    
    cameras = np.unique(frames[:,3])
    
    dates = []
    
    for tstr in frames[:,2]:
        
        ts = datetime.datetime.strptime(tstr,"%Y-%m-%dT%H:%M:%S")
        
        dates.append(ts)
    
    dates = np.array(dates)
    
    start_date = datetime.datetime.strptime('2018-01-01T00:00:00',"%Y-%m-%dT%H:%M:%S")
    end_date = datetime.datetime.strptime('2018-12-31T00:00:00',"%Y-%m-%dT%H:%M:%S")
    dt = datetime.timedelta(days=7.0)
    dbin = datetime.timedelta(days=3.5)
    nbins = 52
    
    plt.figure(1,(10,20))
        
    for c in range(0,len(cameras),1):
        
        ax = plt.subplot(len(cameras),1,c+1)
        
        xdata = []
        ydata1 = []
        ydata2 = []
        ydata3 = []
        
        bin_start = start_date
        while (bin_start <= end_date):
            
            bin_end = bin_start + dt
                        
            nrequested = 0.0
            nsuccess = 0.0
            
            for i in range(0,len(dates),1):
            
                if dates[i] >= bin_start and dates[i] < bin_end and \
                        frames[i,3] == cameras[c]:
                    
                    nrequested += data[i,0]
                    nsuccess += data[i,1]
            
            if nrequested > 0:
                xdata.append( bin_start + dbin )
                ydata1.append( nsuccess )
                ydata2.append( nrequested )
                ydata3.append( (nsuccess/nrequested)*100.0 )
                
            bin_start = bin_end
        print xdata, ydata1, ydata2, ydata3
        
        ax.plot(xdata,ydata1,'g-',label='Obs successful')
        ax.fill_between(xdata, 0, ydata1, facecolor='g',alpha=0.2)
        ax.plot(xdata,ydata2,'m-',label='Obs requested')
        ax.fill_between(xdata, 0, ydata2, facecolor='m',alpha=0.2)
        
        ax.set_ylabel('N obs')
        ax.set_xlabel('Date')
        plt.legend(loc=2)
        
        ax2 = ax.twinx()
        ax2.plot(xdata,ydata3,'k-.',label='Fractional success')
        ax2.fill_between(xdata, 0, ydata3, facecolor='k',alpha=0.1)
        plt.legend(loc=1)
        
        ([xmin,xmax,ymin,ymax]) = plt.axis()
        #xmin = 736406.65
        #xmax= 736645.35
        plt.axis([xmin,xmax,ymin,ymax])
        ax2.set_ylabel('% obs successful')
        
        ax.set_title('Success rate of observation requests for '+cameras[c])
    
    
    plt.subplots_adjust(top=0.92, bottom=0.18, left=0.10, right=0.90, hspace=0.45,
                    wspace=0.15)
                    
    plt.savefig('obs_success_rate.png')
    
    plt.close()
    
        
def get_args():
    
    if len(sys.argv) > 1:
        
        log_path = sys.argv[1]
        
    else:
        
        log_path = raw_input('Please enter the path to the log file: ')
    
    return log_path

def read_obs_execution_log(log_path):
    """Function to parse the query log of the outcome of observations"""
    
    if not os.path.isfile(log_path):
        
        print('ERROR: Cannot find input log file '+log_path)
        
        sys.exit()
        
    f = open(log_path,'r')
    file_lines = f.readlines()
    f.close()
    
    frames = []
    data = []
    
    for line in file_lines:
        
        if line[0:1] != '#':
            
            entries = line.replace('\n','').split()
            
            frames.append( entries[0:4] )
            data.append( [int(entries[4]), int(entries[5])] )
    
    frames = np.array(frames)
    data = np.array(data)
    
    return frames, data
    
if __name__ == '__main__':
    
    analyze_obs_execution()