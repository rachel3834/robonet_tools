# Code to plot the distribution of confirmed exoplanets
# as a function of planet mass vs orbital semi-major axis
# using data from the NASA Exoplanet Archive
#
# By R. Street, with code from Y. Tsapras
import io
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.ticker
import re
import pandas as pd
from astropy.table import Table, Column

class KnownPlanets:
    def __init__(self):
        self.data = np.array([])
        self.methods = ['Radial Velocity', 'Transit', 'Imaging', 'Microlensing']
        self.symbols = ['D', 'o', 's', 'p']
        self.colours = ['#F1A94E', '#5D4C46', '#7B8D8E', '#C603FC']
        self.pointsize = [3,3,3,6]

    def fetch_data_from_archive(self):
        """Method to update the catalogue of data on known exoplanets by
        querying the NASA Exoplanet Archive"""

        url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+discoverymethod,pl_massj,pl_orbsmax,sy_dist,st_mass+from+ps&format=csv'
        querySet = pd.read_csv(url,sep=",") # use sep="," for coma separation.

        self.data = Table([Column(name='method', data=querySet.discoverymethod.values),
                           Column(name='planet_mass', data=querySet.pl_massj.values),
                           Column(name='planet_orb_sep', data=querySet.pl_orbsmax.values),
                           Column(name='distance', data=querySet.sy_dist.values),
                           Column(name='star_mass', data=querySet.st_mass.values)])

    def plot_m_vs_a(self):
        fig = plt.figure(1,(10,10))
        fontsize = 15

        for i,method in enumerate(self.methods):
            idx = np.where((self.data['method'] == method)
                    & (np.isnan(self.data['planet_mass']) == False))[0]
            label = method+': '+str(len(idx))
            plt.plot(self.data['planet_orb_sep'][idx],
                     self.data['planet_mass'][idx],
                     self.symbols[i], markerfacecolor=self.colours[i],
                     alpha=1, markersize=self.pointsize[i],
                     markeredgewidth=0, label=label)
        plt.loglog()
        ax = plt.gca()
        ax.add_patch(Polygon([[0.9,0.0001],[1.3,0.0001],
                              [1.3,100000.0],[0.9,100000.0]],
                              closed=True ,fill=True, alpha=0.7,linewidth=0,
                              label='Snowline'))
        #plt.annotate('snow line', xy=(0.95,0.11), size=12, rotation=90, color='white')

        plt.xlabel(r'Semi-major axis (AU)', size=fontsize)
        plt.ylabel(r'Planet mass ($M_{\oplus}$)', size=fontsize)
        date = str(datetime.datetime.now()).split()[0]
        nplanets = len(np.where(np.isnan(self.data['planet_mass']) == False)[0])
        plt.title(str(nplanets)+' confirmed exoplanets as of '+date,fontsize=13)
        plt.legend(loc=4, fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.grid()
        plt.axis([0.001,10000.0,1e-4,100.0])
        plt.savefig('exoplanets_mp_vs_a.png', bbox_inches='tight')
        plt.close(2)

def plot_mass_separations():
    planets =  KnownPlanets()
    planets.fetch_data_from_archive()
    planets.plot_m_vs_a()

if __name__ == '__main__':
    plot_mass_separations()
