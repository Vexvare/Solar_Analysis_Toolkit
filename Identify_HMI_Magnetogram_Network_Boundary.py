# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 12:23:40 2022

@author : Landon Wells

Purpose : To make HMI magnetogram network boundaries from any given input data.
          
Inputs : The user will have an option to input their own information regarding 
         the time of observation, OR point to where
         the path where their data is and load up the first map from there.
         
Outputs : Maps are saved locally within ./Solar_Analysis_Toolkit/Data/Classification_Methods/CHIMERA/ . . . 
          These maps will be called upon in other programs (if they exist). 
"""

# General Module Imports
import sunpy.map as mp
from sunpy.net import Fido
from sunpy.net import attrs as a
from astropy import units as u
import matplotlib
matplotlib.use('Qt5Agg') 
from aiapy.calibrate import register
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
homedir = str(Path.home())
#import warnings
#warnings.filterwarnings('ignore')

# Local Module Imports
from QOL_Programs.Choose_GUI_togetdata import getting_info
from SolarEventClass.Magnetogram_Boundaries import HMIMagnetogramBoundary

gaussianlevel = 30

def make_network_boundary(time, gaussianlevel, *args):
    
    # Figure settings (to make them look nice)
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['figure.figsize'] = [12.0,12.0]
    plt.rcParams['figure.dpi'] = 80
    plt.rcParams['savefig.dpi'] = 80
    
    # First query the contemporanious AIA image. This image will be used to rescale
    # and reproject the HMI data to make sure that the solar disk is the same size
    # as the AIA solar disk.
    result_aia = Fido.search(time,
                         a.Instrument('AIA'), # Atmospheric Imaging Assembly.
                         a.Wavelength(171*u.angstrom),  # Physical observables
                         a.Sample(1800 * u.s) # Only take a shot every $var seconds.
                         ) #More observables at http://sdac.virtualsolar.org/cgi/show_details?keyword=PHYSOBS
    # Save the results to fits files, keeping in mind the path of the data.
    aia171_name = Fido.fetch(result_aia, site='ROB')
    # Use the path to make a sunpy map.
    aia171_ = mp.Map(aia171_name, sequence = True)
    
    # Now we query the contemporanious HMI data.
    # HMI
    result_hmi = Fido.search(time,
                          a.Instrument('HMI'), # Helioseismic and Magnetic Imager.
                          a.Physobs('LOS_magnetic_field'),
                          a.Sample(1800 * u.s)) 
    hmi_name = Fido.fetch(result_hmi[0])
    hmi_ = mp.Map(hmi_name, sequence = True)
    
    # Align and scale the HMI to the aia171 level 1.5 image.
    # Enumerate through the list of maps . . .
    hmi_reproj_list = []
    print('Registering and scaling the HMI maps to match the AIA maps. This will take some time.')
    for i, hmi_map in enumerate(hmi_):
        
        aia171 = register(aia171_[i])
        hmi = hmi_map.reproject_to(aia171.wcs)
        hmi_reproj_list.append(hmi)
        
    # Now we create the HMI Magnetic Network Boundary using a simple gaussian 
    # thresholding techinque.
    hmi = mp.Map(hmi_reproj_list, sequence = True)
    mag_boundary =  HMIMagnetogramBoundary(hmi, classification_method = 'GaussianContours' , GaussianLevel = gaussianlevel)
    mag_boundary_map = mag_boundary.magnetogram_boundary_maps
    
    # Now we show the sequence of maps which make up the magnetogram boundaries.
    fig = plt.figure(1)
    ani = mag_boundary_map.plot()
    
    # Middle image hmi full-disk with overlay of boundary.
    fig = plt.figure(2)
    ax = plt.subplot(projection = hmi[int(len(hmi)/2)])
    image = hmi[int(len(hmi)/2)].plot()
    hmi[int(len(hmi)/2)].draw_grid()
    plt.title('')
    boundary_labels = ['HMI Network Boundary']
    boundary_mag = plt.contour(mag_boundary_map[int(len(hmi)/2)]._data, levels = [0.5,1.5], alpha = 0.85,
                                 color = 'black', cmap = cm.cool, linewidths = 0.5)

    boundary_mag.collections[0].set_label(boundary_labels[0])
    
    # Show all of the plots
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            plt.close(1)
            plt.close(2)
            break
    

if __name__ == '__main__':
    
    # First we ask the user to input/search for data . . . 
    info = getting_info()
    instrument = info.telescope
    wavelength = info.wavelength
    starttime = info.starttime
    endtime = info.endtime
    time = a.Time(starttime, endtime)
    
    make_network_boundary(time, gaussianlevel)
    