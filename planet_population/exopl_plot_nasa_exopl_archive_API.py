#!/usr/bin/python

##########################################################
# Make a plot of the distrubition of exoplanets in Orbit
# size (AU) vs planet mass (M_earth).
# using the latest planet data from exoplanet.eu
#
# Requires: matplotlib, urllib2, HTMLParser, cStringIO, numpy
#           datetime
# version 1.0 by Y Tsapras 4 Feb 2011
# last modified 17 Mar 2016 by YT
#               11 Jan 2012 by YT - corrected bug with
#                     information parsing from website
#               15 Jan 2014 by YT - updated script to work with
#                   the new format of the exoplanet.eu website
#               7 Mar 2014 by YT - updated to record distances
#                   and stellar (host) masses
#               11 Apr 2014 by YT - Updated plots for snow-line
#               12 Jan 2015 by YT - Updated to work with CSV files
#               15 Feb 2017 by YT - Updated to work with single CSV file
#               15 May 2018 by YT - Updated to work with CSV file from NASA Exoplanet Archive
#               14 Sep 2020 by YT - Updated to use API
##########################################################

# Import required modules
import io
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.ticker
import re
import pandas as pd

# Get current date
now = datetime.datetime.now()
date_gen = str(now).split()[0]

# Define holding lists for parameters of all methods
# Define holding lists for parameters of all methods
mass_ml, orbit_ml, dist_ml, starmass_ml, snow_ml = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
mass_tr, orbit_tr, dist_tr, starmass_tr, snow_tr = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
mass_ti, orbit_ti, dist_ti, starmass_ti, snow_ti = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
mass_rv, orbit_rv, dist_rv, starmass_rv, snow_rv = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
mass_im, orbit_im, dist_im, starmass_im, snow_im = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

# Reset the number of planets found by each method
ml_plan, tr_plan, ti_plan, rv_plan, im_plan= 0, 0, 0, 0, 0

# Get the data by querying the NASA Exoplanet Archive API for the given column identifiers
# RAS: Updated following DB change of column names
#url ='https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=pl_discmethod,pl_bmassj,pl_orbsmax,st_dist,st_mass&format=csv'
url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+discoverymethod,pl_massj,pl_orbsmax,sy_dist,st_mass+from+ps&format=csv'
data = pd.read_csv(url,sep=",") # use sep="," for coma separation.


# You can use these columns:
# COLUMN pl_discmethod:  Discovery Method
# COLUMN pl_orbper:      Orbital Period [days]
# COLUMN pl_orbsmax:     Orbit Semi-Major Axis [AU])
# COLUMN pl_bmassj:      Planet Mass or M*sin(i) [Jupiter mass]
# COLUMN st_dist:        Distance [pc]
# COLUMN st_mass:        Stellar Mass [Solar mass]
# COLUMN pl_massj:       Planet Mass [Jupiter mass]


# Extract the information about the mass
# (converting to Earth masses), the semi-major axis in AU,
# The distance in pc and the stellar mass in solar masses.

# Update every list accordingly
method = data.discoverymethod.values
pl_mass = data.pl_massj.values
pl_orbit = data.pl_orbsmax.values
st_dist =  data.sy_dist.values
st_mass = data.st_mass.values

for i in np.arange(len(pl_mass)):
   # start going through each line, skipping lines that contain no useful info
   str_method = method[i]
   str_mass = pl_mass[i]
   str_orbit = pl_orbit[i]
   str_distance = st_dist[i]
   str_starmass = st_mass[i]
   if not str_mass: str_mass = '1000000.'
   if not str_orbit: str_orbit ='-1000.'
   if not str_distance: str_distance ='-1000.'
   if not str_starmass: str_starmass ='0.00000000000001'
   # Append the values to the lists
   if 'Microlensing' in str_method:
      mass_ml = np.append(mass_ml, 317.83 * float(str_mass)) # in Earth masses
      orbit_ml = np.append(orbit_ml, float(str_orbit)) # in AU
      dist_ml = np.append(dist_ml, float(str_distance))
      starmass_ml = np.append(starmass_ml, float(str_starmass)) # In solar masses
      snow_ml = np.append(snow_ml, float(str_orbit) / (2.848 * float(str_starmass)**(3./2)) )
      ml_plan = ml_plan + 1
   if 'Imaging' in str_method:
      mass_im = np.append(mass_im, 317.83 * float(str_mass)) # in Earth masses
      orbit_im = np.append(orbit_im, float(str_orbit)) # in AU
      if (float(str_distance) >= 1000): str_distance ='-1000.'
      if (float(str_starmass) <= 0.05): str_starmass ='0.00000000000001'
      dist_im = np.append(dist_im, float(str_distance))
      starmass_im = np.append(starmass_im, float(str_starmass)) # In solar masses
      snow_im = np.append(snow_im, float(str_orbit) / (2.848 * float(str_starmass)**(3./2)) )
      im_plan = im_plan + 1
   if 'Pulsar Timing' in str_method:
      mass_ti = np.append(mass_ti, 317.83 * float(str_mass)) # in Earth masses
      orbit_ti = np.append(orbit_ti, float(str_orbit)) # in AU
      dist_ti = np.append(dist_ti, float(str_distance))
      starmass_ti = np.append(starmass_ti, float(str_starmass)) # In solar masses
      snow_ti = np.append(snow_ti, float(str_orbit) / (2.848 * float(str_starmass)**(3./2)) )
      ti_plan = ti_plan + 1
   if 'Transit' in str_method:
      if (float(str_distance) >= 5000): str_distance ='-1000.'
      mass_tr = np.append(mass_tr, 317.83 * float(str_mass)) # in Earth masses
      orbit_tr = np.append(orbit_tr, float(str_orbit)) # in AU
      dist_tr = np.append(dist_tr, float(str_distance))
      starmass_tr = np.append(starmass_tr, float(str_starmass)) # In solar masses
      snow_tr = np.append(snow_tr, float(str_orbit) / (2.848 * float(str_starmass)**(3./2)) )
      tr_plan = tr_plan + 1
   if 'Radial Velocity' in str_method:
      mass_rv = np.append(mass_rv, 317.83 * float(str_mass)) # in Earth masses
      orbit_rv = np.append(orbit_rv, float(str_orbit)) # in AU
      dist_rv = np.append(dist_rv, float(str_distance))
      starmass_rv = np.append(starmass_rv, float(str_starmass)) # In solar masses
      snow_rv = np.append(snow_rv, float(str_orbit) / (2.848 * float(str_starmass)**(3./2)) )
      rv_plan = rv_plan + 1

#-------------------------------------------------------------------------------------#
# Generate the mass-vs-orbital_distance (AU) plot
plt.figure(2, edgecolor="k", figsize=[9,9])

plt.loglog()
#plt.grid(True,which="majorminor",ls="-")
#plt.xlim((0.01,200.0))
#plt.ylim((0.003,12000.0))
# set up the legends
rvleg = 'Radial velocity: '+str(rv_plan)
trleg = 'Transits: '+str(tr_plan)
imleg = 'Imaging: '+str(im_plan)
mlleg = 'Microlensing: '+str(ml_plan)
tileg = 'Timing: '+str(ti_plan)

# Overplot sensitivities
#ax = plt.gca()
# Transits from space
#ax.add_patch(Polygon([[0.001,0.1],[0.9,0.9],[0.9,100000.0],[0.001,100000.0]], closed=True ,fill=True, alpha=0.5,linewidth=0,color='#FF6699'))
#plt.annotate('transits (space)', xy=(0.003,0.5), size=14, rotation=0, color='#990000')
# Radial velocity
#ax.add_patch(Polygon([[0.001,0.9],[0.05,1.0],[0.08,1.4],[0.1,1.5],[0.45,3.5],[0.6,5],[0.8,7],[0.9,8.0],[1.2,10],[2,20],[3,70.0],[4,100],[7,700.0],[8,100000.0],[0.001,100000.0]], closed=True ,fill=True, alpha=0.9,linewidth=0,color='#808080'))
#plt.annotate('doppler', xy=(0.0014,2), size=14, rotation=0, color='#080808')
# Transits from ground
#ax.add_patch(Polygon([[0.001,90.0],[0.09,90.0],[0.09,100000.0],[0.001,100000.0]], closed=True ,fill=True, alpha=1,linewidth=0,color='#CC6699'))
#plt.annotate('transits', xy=(0.002,140), size=14, rotation=0, color='#990000')
#plt.annotate('(ground)', xy=(0.002,85), size=14, rotation=0, color='#990000')
# Microlensing from space
#ax.add_patch(Polygon([[0.1,12000.0],[0.6,0.8],[0.8,0.26],[0.9,0.24],[2.0,0.08],[4,0.06],[10,0.09],[100,0.1],[1100,0.1],[12000,12000],[0.1,12000]], closed=True ,fill=True, alpha=.3,linewidth=0,color='#00FF99'))
#plt.annotate('microlensing (space)', xy=(10,0.2), size=14, rotation=0, color='#006600')
#Imaging
#ax.add_patch(Polygon([[2,400],[10000.0,400],[10000.0,100000.0],[1.1,100000.0]], closed=True ,fill=True, alpha=.5,linewidth=0,color='#9900CC'))
#plt.annotate('imaging', xy=(100,500), size=14, rotation=0, color='#000099')
# Microlensing
#ax.add_patch(Polygon([[0.11,12000.0],[0.2,300],[0.3,95],[0.5,7],[0.8,1.9],[1,0.9],[4,0.9],[5,1.2],[6,1.8],[8,4],[10,6],[30,80],[50,300],[100,12000]], closed=True ,fill=True, alpha=.3,linewidth=0,color='#FF6600'))
#plt.annotate('microlensing', xy=(1.5,40), size=14, rotation=0, color='#006600')
#plt.annotate('(ground)', xy=(1.7,25), size=14, rotation=0, color='#006600')

# plot the survey findings (generic color #F1A94E)
plt.plot(orbit_rv, mass_rv,'D', markerfacecolor='#F1A94E', alpha=1, markersize=6, markeredgewidth=0, label=rvleg)
plt.plot(orbit_tr, mass_tr,'o', markerfacecolor='#5D4C46', alpha=1, markersize=6, markeredgewidth=0, label=trleg)
plt.plot(orbit_im, mass_im,'s', markerfacecolor='#7B8D8E', alpha=1, markersize=6, markeredgewidth=0, label=imleg)
plt.plot(orbit_ml, mass_ml,'p', markerfacecolor='#44B3C2', alpha=1, markersize=8, markeredgewidth=0, label=mlleg)
#plt.plot(orbit_rv, mass_rv,'o', markerfacecolor='#a64f08', alpha=1, markersize=6, markeredgewidth=0, label=rvleg)
#plt.plot(orbit_tr, mass_tr,'o', markerfacecolor='#d40606', alpha=1, markersize=6, markeredgewidth=0, label=trleg)
#plt.plot(orbit_im, mass_im,'o', markerfacecolor='#80e874', alpha=1, markersize=6, markeredgewidth=0, label=imleg)
#plt.plot(orbit_ml, mass_ml,'o', markerfacecolor='#fffd6e', alpha=1, markersize=6, markeredgewidth=0, label=mlleg)

# plot the positions of the solar planets
#plt.plot(0.39,0.055,'o', markerfacecolor='#8F8A78', alpha=1, markersize=13, markeredgewidth=0) # Mercury
#plt.text(0.366,0.05,'M',size=11)
plot_solsys = False
if plot_solsys:
    plt.plot(0.723,0.815,'o', markerfacecolor='#EACF60', alpha=1, markersize=13, markeredgewidth=0) # Venus
    plt.text(0.66,0.73,'V',size=11)
    plt.plot(1.,1.,'o', markerfacecolor='#48A8F2', alpha=1, markersize=13, markeredgewidth=0) # Earth
    plt.text(0.92,0.9,'E',size=11)
    plt.plot(1.52,0.107,'o', markerfacecolor='#C24024', alpha=1, markersize=13, markeredgewidth=0) # Mars
    plt.text(1.37,0.095,'M',size=11)
    plt.plot(5.203,318.0,'o', markerfacecolor='#FF5D00', alpha=1, markersize=18, markeredgewidth=0) # Jupiter
    plt.text(5.15,299.0,'J',size=11)
    plt.plot(9.536,93.6,'o', markerfacecolor='#A080B0', alpha=1, markersize=17, markeredgewidth=0) # Saturn
    plt.text(8.9,83.6,'S',size=11)
    plt.plot(19.18,15.0,'o', markerfacecolor='#2D59B0', alpha=1, markersize=15, markeredgewidth=0) # Uranus
    plt.text(17.5,13.6,'U',size=11)
    plt.plot(30.06,17.0,'o', markerfacecolor='#0400FF', alpha=1, markersize=15, markeredgewidth=0) # Neptune
    plt.text(27.6,15.2,'N',size=11)

# increase the size of the text that labels the tick marks on the axes
xticklabels = plt.getp(plt.gca(), 'xticklabels')
yticklabels = plt.getp(plt.gca(), 'yticklabels')
plt.setp(xticklabels, fontsize=16)
plt.setp(yticklabels, fontsize=16)
plt.xlabel(r'Semi-major axis (AU)', size=15)
plt.ylabel(r'Planet mass ($M_{\oplus}$)', size=15)
plt.title('Exoplanets as of '+date_gen,fontsize=13)
plt.setp(xticklabels, fontsize=16)
plt.setp(yticklabels, fontsize=16)

# set the legend font size to use
legend_font_par = {'legend.fontsize': 16}

plt.rcParams.update(legend_font_par)

#plt.legend(bbox_to_anchor=(0.0,1.02,1.0,0.102), ncol=4, loc=3, mode="expand", numpoints = 1,
#          labelspacing=0.02, borderpad=0.01, borderaxespad=0.01, handletextpad=0.3,
#          handlelength=0.2)

leg = plt.legend(loc=4)

#leg.get_frame().set_linewidth(2)

plt.axes().set_aspect('equal')
# Use scalar formatting for axes if needed (or leave commented out for power
#ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#plt.xticks(ticks=[0.01, 0.1,1,10,100,1000], labels=[0.01, 0.1,1,10,100,1000])
#plt.yticks(ticks=[0.1,1,10,100,1000,10000], labels=[0.1,1,10,100,1000,10000])
plt.xlim((0.01,1000.0))
plt.ylim((0.1,10000.0))    # set the legend font size to use
plt.grid(True)
#plt.show()

#x = np.ones(100)
#plt.plot(x,y,'b--',linewidth=3)
#plt.show()

print("Completed processing. Printing figure...")
# Save the figure as png
plt.savefig('exoplanets.png', dpi = (90) )
plt.close(2)
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Generate the mass-vs-orbital_distance (AU/snow_line) plot
plt.figure(3, edgecolor="k", figsize=[9,9])
#plt.figure(1, edgecolor="k", figsize=[20,10])

plt.loglog()
#plt.grid(True,which="majorminor",ls="-")
#plt.xlim((0.01,200.0))
#plt.ylim((0.003,12000.0))
# set up the legends
rvleg = 'Radial Velocity'
trleg = 'Transits'
imleg = 'Direct Imaging'
mlleg = 'Microlensing'
tileg = 'Timing'

# Overplot sensitivities
ax = plt.gca()

# Transits from space
#ax.add_patch(Polygon([[0.01,0.1],[0.4,0.9],[0.4,100000.0],[0.01,100000.0]], closed=True ,fill=True, alpha=0.5,linewidth=0,color='#FF6699'))
#plt.annotate('transits (space)', xy=(0.02,1.05), size=14, rotation=0, color='#990000')
#plt.annotate('transits (space)', xy=(0.02,1.05), size=14, rotation=0, color='#990000')
# Radial velocity
#ax.add_patch(Polygon([[0.001,0.9],[0.03,1.0],[0.1,1.5],[0.3,2.5],[0.5,6],[0.6,8],[0.8,15],[0.9,30.0],[1,60],[1.2,80],[2,200],[3,600],[4,100000.0],[0.001,100000.0]], closed=True ,fill=True, alpha=0.9,linewidth=0,color='#808080'))
#plt.annotate('Radial Velocity', xy=(0.05,210), size=16, rotation=0, color='#080808')
# Transits from ground
#ax.add_patch(Polygon([[0.001,90.0],[0.04,90.0],[0.04,100000.0],[0.001,100000.0]], closed=True ,fill=True, alpha=1,linewidth=0,color='#CC6699'))
#plt.annotate('transits', xy=(0.002,140), size=14, rotation=0, color='#990000')
#plt.annotate('Transits', xy=(0.005,35), size=16, rotation=0, color='#990000')
# Microlensing from space
#ax.add_patch(Polygon([[0.1,12000.0],[0.22,0.5],[0.8,0.08],[0.9,0.06],[2.0,0.05],[4,0.06],[10,0.08],[100,0.2],[1100,0.2],[12000,12000],[0.1,12000]], closed=True ,fill=True, alpha=.3,linewidth=0,color='#00FF99'))
#plt.annotate('microlensing (space)', xy=(10,0.45), size=14, rotation=0, color='#006600')
#Imaging
#ax.add_patch(Polygon([[2,400],[10000.0,400],[10000.0,100000.0],[1.1,100000.0]], closed=True ,fill=True, alpha=.5,linewidth=0,color='#9900CC'))
#plt.annotate('Direct Imaging', xy=(20,2000), size=16, rotation=0, color='#000099')
# Microlensing
#ax.add_patch(Polygon([[0.3,12000.0],[0.4,200],[1,10],[3,0.9],[8,0.9],[20,2.2],[110,1000],[200,12000.0]], closed=True ,fill=True, alpha=.3,linewidth=0,color='#FF6600'))
#plt.annotate('Microlensing', xy=(3,22), size=16, rotation=0, color='#006600')
#plt.annotate('Microlensing', xy=(1.9,25), size=14, rotation=0, color='#006600')

# Plot snow-line
ax = plt.gca()
ax.add_patch(Polygon([[0.9,0.0001],[1.3,0.0001],[1.3,100000.0],[0.9,100000.0]], closed=True ,fill=True, alpha=1,linewidth=0))
plt.annotate('snow line', xy=(0.95,0.11), size=12, rotation=90, color='white')

# plot the survey findings
plt.plot(snow_rv, mass_rv,'D', markerfacecolor='#F1A94E', alpha=1, markersize=6, markeredgewidth=0, label=rvleg)
plt.plot(snow_tr, mass_tr,'o', markerfacecolor='#5D4C46', alpha=1, markersize=6, markeredgewidth=0, label=trleg)
plt.plot(snow_im, mass_im,'s', markerfacecolor='#7B8D8E', alpha=1, markersize=6, markeredgewidth=0, label=imleg)
plt.plot(snow_ml, mass_ml,'p', markerfacecolor='#44B3C2', alpha=1, markersize=8, markeredgewidth=0, label=mlleg)

#plt.plot(snow_ti, mass_ti,'*', markerfacecolor='#FFCC00', alpha=1, markersize=8, markeredgewidth=0, label=tileg)
#plt.plot(snow_rv, mass_rv,'o', markerfacecolor='#a64f08', alpha=1, markersize=6, markeredgewidth=0, label=rvleg)
#plt.plot(snow_tr, mass_tr,'o', markerfacecolor='#d40606', alpha=1, markersize=6, markeredgewidth=0, label=trleg)
#plt.plot(snow_im, mass_im,'o', markerfacecolor='#80e874', alpha=1, markersize=6, markeredgewidth=0, label=imleg)
#plt.plot(snow_ml, mass_ml,'o', markerfacecolor='#fffd6e', alpha=1, markersize=6, markeredgewidth=0, label=mlleg)

# plot the positions of the solar planets
#plt.plot(0.39,0.055,'o', markerfacecolor='#8F8A78', alpha=1, markersize=13, markeredgewidth=0) # Mercury
#plt.text(0.366,0.05,'M',size=11)
plt.plot(float(0.723) / (2.848 * float(1.0)**(3./2)),0.815,'o', markerfacecolor='#EACF60', alpha=1, markersize=13, markeredgewidth=0) # Venus
plt.text(float(0.660) / (2.848 * float(1.0)**(3./2)),0.73,'V',size=11)
plt.plot(float(1.000) / (2.848 * float(1.0)**(3./2)),1.,'o', markerfacecolor='#48A8F2', alpha=1, markersize=13, markeredgewidth=0) # Earth
plt.text(float(0.920) / (2.848 * float(1.0)**(3./2)),0.9,'E',size=11)
plt.plot(float(1.520) / (2.848 * float(1.0)**(3./2)),0.107,'o', markerfacecolor='#C24024', alpha=1, markersize=13, markeredgewidth=0) # Mars
plt.text(float(1.370) / (2.848 * float(1.0)**(3./2)),0.095,'M',size=11)
plt.plot(float(5.203) / (2.848 * float(1.0)**(3./2)),318.0,'o', markerfacecolor='#FF5D00', alpha=1, markersize=18, markeredgewidth=0) # Jupiter
plt.text(float(5.150) / (2.848 * float(1.0)**(3./2)),299.0,'J',size=11)
plt.plot(float(9.536) / (2.848 * float(1.0)**(3./2)),93.6,'o', markerfacecolor='#A080B0', alpha=1, markersize=17, markeredgewidth=0) # Saturn
plt.text(float(8.900) / (2.848 * float(1.0)**(3./2)),83.6,'S',size=11)
plt.plot(float(19.18) / (2.848 * float(1.0)**(3./2)),15.0,'o', markerfacecolor='#2D59B0', alpha=1, markersize=15, markeredgewidth=0) # Uranus
plt.text(float(17.50) / (2.848 * float(1.0)**(3./2)),13.6,'U',size=11)
plt.plot(float(30.06) / (2.848 * float(1.0)**(3./2)),17.0,'o', markerfacecolor='#0400FF', alpha=1, markersize=15, markeredgewidth=0) # Neptune
plt.text(float(27.60) / (2.848 * float(1.0)**(3./2)),15.2,'N',size=11)

# increase the size of the text that labels the tick marks on the axes
# increase the size of the text that labels the tick marks on the axes
xticklabels = plt.getp(plt.gca(), 'xticklabels')
yticklabels = plt.getp(plt.gca(), 'yticklabels')
plt.setp(xticklabels, fontsize=16)
plt.setp(yticklabels, fontsize=16)
plt.xlabel(r'Semi-major axis/Snow line', size=18)
plt.ylabel(r'Planet mass ($M_{\oplus}$)', size=18)
plt.xlim((0.003,800.0))
plt.ylim((0.08,11000.0))

# set the legend font size to use
#legend_font_par = {'legend.fontsize': 13}
#plt.rcParams.update(legend_font_par)
#plt.legend(bbox_to_anchor=(0.0,1.07,1.0,0.102), ncol=5, loc=3, mode="expand", numpoints = 1,
#           labelspacing=0.01, borderpad=0.1, borderaxespad=0.1, handletextpad=0.2,
#	   handlelength=0.4)
#plt.grid(True)
#plt.annotate('transits (space)', xy=(0.02,1.05), size=14, rotation=0, color='#990000')
#plt.annotate('microlensing (space)', xy=(10,0.45), size=14, rotation=0, color='#006600')
#plt.annotate('microlensing', xy=(1.8,40), size=14, rotation=0, color='#006600')
#plt.annotate('(ground)', xy=(1.9,25), size=14, rotation=0, color='#006600')

plt.axes().set_aspect('equal')

# set the legend font size to use
legend_font_par = {'legend.fontsize': 16}

plt.rcParams.update(legend_font_par)

#plt.legend(bbox_to_anchor=(0.0,1.02,1.0,0.102), ncol=4, loc=3, mode="expand", numpoints = 1,
#          labelspacing=0.02, borderpad=0.01, borderaxespad=0.01, handletextpad=0.3,
#          handlelength=0.2)

leg = plt.legend(loc=4)

leg.get_frame().set_linewidth(2)

#plt.show()

print("Completed processing. Printing figure...")
# Save the figure as png
plt.savefig('exoplanets_snowline.png', dpi = (90) )
#-------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------#
# Plot distances
plt.figure(4, edgecolor="k", figsize=[14,6])

plt.loglog()
#plt.grid(True,which="majorminor",ls="-")
#plt.xlim((0.01,200.0))
#plt.ylim((0.003,12000.0))
# set up the legends
rvleg = 'Radial Velocity'
trleg = 'Transits'
imleg = 'Direct Imaging'
mlleg = 'Microlensing'
tileg = 'Timing'

# plot the survey findings
plt.plot(dist_rv, mass_rv,'D', markerfacecolor='#F1A94E', alpha=1, markersize=6, markeredgewidth=0, label=rvleg)
plt.plot(dist_tr, mass_tr,'o', markerfacecolor='#5D4C46', alpha=1, markersize=6, markeredgewidth=0, label=trleg)
plt.plot(dist_im, mass_im,'s', markerfacecolor='#7B8D8E', alpha=1, markersize=6, markeredgewidth=0, label=imleg)
plt.plot(dist_ml, mass_ml,'p', markerfacecolor='#44B3C2', alpha=1, markersize=8, markeredgewidth=0, label=mlleg)
#plt.plot(dist_ti, mass_ti,'*', markerfacecolor='#FFCC00', alpha=.8, markersize=8, markeredgewidth=0, label=tileg)

# increase the size of the text that labels the tick marks on the axes
xticklabels = plt.getp(plt.gca(), 'xticklabels')
yticklabels = plt.getp(plt.gca(), 'yticklabels')
plt.xlabel(r'Distance (pc)', size=16)
plt.ylabel(r'Planet Mass ($M_{\oplus}$)', size=16)
#plt.title('Exoplanets as of '+date_gen,fontsize=13)
plt.setp(xticklabels, fontsize=16)
plt.setp(yticklabels, fontsize=16)
#xticksize=plt.getp(plt.gca(), 'xticklines')
#yticksize=plt.getp(plt.gca(), 'yticklines')
#plt.setp(xticksize, markersize=12)
#plt.setp(yticksize, markersize=12)
plt.xlim((0.01,1200.0))
plt.ylim((0.09,12000.0))
# set the legend font size to use
plt.legend(loc=2)
plt.xlim((1.0,12000.0))
plt.ylim((0.09,12000.0))

#x = np.ones(100)
#y = np.arange(0.001,10000,100)
#plt.plot(x,y,'b--',linewidth=3)

print("Completed processing. Printing figure...")
# Save the figure as png
plt.savefig('exoplanets_distance.png', dpi = (60) )
#-------------------------------------------------------------------------------------#
