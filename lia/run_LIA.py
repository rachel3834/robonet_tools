import numpy as np
import matplotlib.pyplot as plt
from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA import microlmodels
from LIA import microlensing_classifier
from LIA import models
from LIA import training_set
from sklearn.manifold import TSNE as tsne
from sklearn.metrics import confusion_matrix
from pyDANDIA import lightcurves
from sys import argv
import sqlite3
from os import getcwd, path, remove, environ
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units, table
import matplotlib.pyplot as plt
#from pyDANDIA import  phot_db
#from pyDANDIA import  hd5_utils
#from pyDANDIA import  pipeline_setup
#from pyDANDIA import  metadata
#from pyDANDIA import  logs
import matplotlib
matplotlib.use('TkAgg')
import h5py
import pandas as pd
from pyDANDIA import  metadata
from pyLIMA import microlsimulator
from astropy.coordinates import EarthLocation
import astropy
from astropy.coordinates import AltAz

from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA import microlmodels


def get_Moon(times,ra=270,dec=-30):

    moon = []
    illumination = []

    for time in times:

        ti = astropy.time.Time(time,format='jd')

        aa = astropy.coordinates.get_moon(ti)
        bb = astropy.coordinates.get_sun(ti)
       
        moon.append([aa.ra.value,aa.dec.value])
        illumination.append(microlsimulator.moon_illumination(bb,aa).value)

    Moon = np.array(moon)
    distance = np.arccos(np.sin(Moon[:,1]*3.14/180)*np.sin(dec*3.14/180)+np.cos(Moon[:,1]*3.14/180)*np.cos(dec*3.14/180)*np.cos((Moon[:,0]-ra)*3.14/180))*180/3.14
    
    return Moon,distance,np.array(illumination)
    
### Setup paths etc    
phot_path = '/home/ebachelet/Work/ROMEREA/LIA/ROME1/LSCA_ip/'
hdf_files = h5py.File(phot_path+'photometry.hdf5', 'r')


reduction_metadata = metadata.MetaData()

reduction_metadata.load_all_metadata(metadata_directory=phot_path,
                                     metadata_name='pyDANDIA_metadata.fits')


target = SkyCoord(ra=267.83589537*units.deg,dec=-30.060817819*units.deg)


### Extract images quality metrics


pscale =     reduction_metadata.images_stats[1]['PSCALE']
sky =     reduction_metadata.images_stats[1]['SKY']
fwhm =     reduction_metadata.images_stats[1]['FWHM']
nstars =     reduction_metadata.images_stats[1]['NSTARS']
shift_x =     reduction_metadata.images_stats[1]['SHIFT_X']
shift_y =     reduction_metadata.images_stats[1]['SHIFT_Y']
sig_x =     reduction_metadata.images_stats[1]['SIGMA_X']
sig_y =     reduction_metadata.images_stats[1]['SIGMA_Y']
distance = reduction_metadata.headers_summary[1]['MOONDKEY'].data.astype(float)
illumination = reduction_metadata.headers_summary[1]['MOONFKEY'].data.astype(float)
exptime = reduction_metadata.headers_summary[1]['EXPKEY'].data.astype(float)
TIME = astropy.time.Time(reduction_metadata.headers_summary[1]['DATEKEY'].data).jd
ind_ref = np.where(reduction_metadata.images_stats[1]['IM_NAME']==reduction_metadata.data_architecture[1]['REF_IMAGE'])[0][0]


### Compare time with what in the metadata (old pieces)
#LSC = EarthLocation.of_site('Cerro Tololo')
#dico_time = []
#time = hdf_files['dataset_photometry'][10][:,9]
#altaz = AltAz(location=COJ, obstime=astropy.time.Time(time,format='jd'))
#airmass = target.transform_to(altaz).secz
#for index in np.arange(0,len(pscale)):

#    if time[int(index)]>0:
#        
#        dist = np.abs(TIME-time[int(index)])
#        
#        ind = dist.argmin()
#        dico_time.append([int(index), ind])
#    else:   
#        dico_time.append([0, 0])

#dico_time = np.array(dico_time)

times = []
for lc in np.random.randint(0,2000,100):

    time = hdf_files['dataset_photometry'][lc][:,9]
    mag = hdf_files['dataset_photometry'][lc][:,11]
    emag = hdf_files['dataset_photometry'][lc][:,12]
    try:
        mask = (time>1) &  (exptime[dico_time[:,1]]>200) #& (mag>0) & (pscale>0)  & (np.median(mag)<19) & (Illumination<1.0) #& (fwhm[dico_time[:,1]]<6) #& (distance>40) #& (pscale>0.8) & (pscale<1.2)
    
    except:
        #Moon,Distance,Illumination = get_Moon(time,ra=267.83589537,dec=-30.060817819) 
        mask = (time>0) #& (mag>0) & (pscale>0) &  (fwhm[dico_time[:,1]]<6) 
    if len(time[mask])>5:
        times.append(np.sort(time[mask]))



### Comment this when done once
#n_class = 1000
#training_set.create(times, min_mag=10, max_mag=20, noise=None, cv_n1=2,n_class=n_class,test=False)


model = models.create_models('all_features.txt', model='rf')

### ploting TSNE if needed

#pcapca = np.loadtxt('pca_features.txt',dtype=str)
#vis_data = tsne(n_components=2).fit_transform(pcapca[:,2:].astype(float))
#vis_x = vis_data[:,0]
#vis_y = vis_data[:,1]

#colors = ['red', 'green', 'blue', 'yellow', 'orange', 'purple', 'black']
#fig, ax = plt.subplots()
#for i,object_class in enumerate(['Variable', 'Constant', 'CV', 'ML']):
#    x=vis_x[i*n_class:(i+1)*n_class]
#    y=vis_y[i*n_class:(i+1)*n_class]
#    ax.scatter(x, y, c=colors[i], label=object_class)

#plt.legend()
#plt.show()


### Testing on ROME data
These lc should have be some ML events, including OB190011
for lc in [207194,132119,177748,177461,78283,121424,174610,215315]:



#randi = np.arange(30000,200000,1)
for lc in randi:

    lc = int(lc)
    print(lc)
    try:
        time = hdf_files['dataset_photometry'][lc][:,9]
        mag = hdf_files['dataset_photometry'][lc][:,11]
        emag = hdf_files['dataset_photometry'][lc][:,12]
        back = hdf_files['dataset_photometry'][lc][:,-2]
        ppscale = hdf_files['dataset_photometry'][lc][:,-4]
        
        mask = (time>1) & (mag>1) & (np.abs(ppscale-exptime[ind_ref]/exptime)<0.2) & (np.abs(back)<250)
        order = time[mask].argsort()
        if np.median(mag[mask])<30:
            classification= microlensing_classifier.predict(time[mask][order],mag[mask][order],emag[mask][order], model)
            if float(classification[3][1])>0.6:
                plt.scatter(time[mask],mag[mask])
                plt.gca().invert_yaxis()
                plt.show()
                import pdb; pdb.set_trace()
    except:
        pass
import pdb; pdb.set_trace()

