# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 14:19:36 2022

@author : Landon Wells

Purpose : Create Coronal Hole Boundaries using the CHARM classification
          method.
          
Inputs : The user will have an option to input their own information regarding 
         the time of observation, OR point to where
         the path where their data is and load up the first map from there.
         
Outputs : Maps are saved locally within ./Solar_Analysis_Toolkit/Data/Classification_Methods/CHARM/ . . . 
          These maps will be called upon in other programs (if they exist). 
"""

# General imports
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

# Module Imports
from QOL_Programs.Choose_GUI_togetdata import getting_info
from SolarEventClass.CoronalHoleBoundaries import CoronalHoleBoundary



def make_charm_boundary(time, *args):
    
    # Figure settings (to make them look nice)
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['figure.figsize'] = [12.0,12.0]
    plt.rcParams['figure.dpi'] = 80
    plt.rcParams['savefig.dpi'] = 80
    
    # if args is given, it is assumed to be the maps that were queried previously.
    
    # But for all intensive purposes, CHARM only the 193 AIA data. So regardless 
    # of the previous query, we use the time of interest to download and search 
    # for data in this specific wavelength.
    
    # Build the query via Fido, this can return one item, or a list of them to DL matching
    # the given filters. Full documentation on Fido and it's attributes can be found on the Sunpy website.
    
    # Obtain the corresponding AIA data showing the EUV images at the time of interest.
    
    # Searching for the 193 data via Sunpy Fido.

    result_aia = Fido.search(time,
                         a.Instrument('AIA'),
                         a.Wavelength(193*u.angstrom),
                         a.Sample(1800 * u.s)) 
    aia193_name = Fido.fetch(result_aia, site='ROB')
    # We use the header information of the first map, but the data
    # information of the middle map. For my purposes,
    # this will be used specifically to overlay on
    # time series, and at any point can be changed by the user. 
    aia193_maps = mp.Map(aia193_name, sequence = True)
    
    # Align and scale the aia171 level 1.5 image.
    # Enumerate through the list of maps . . .
    aia_list = []
    for i, aia_map in enumerate(aia193_maps):
        aia_list.append(register(aia193_maps[i]))
    aia193_ = mp.Map(aia_list, sequence = True)

    # CHARM Boundary 
    # We are going to call the Solar_Events_Class and define a CH Boundary 
    # event that uses the CHARM boundary. Calling this class will save the 
    # CHARM boundary locally, so that other programs may make use of it.
    
    boundary = CoronalHoleBoundary(aia193_, 
                                   classification_method = 'CHARM', 
                                   intensitylevel = 63)
    charm_boundary = boundary.charm_map
    
    
    # Plotting middle image full-disk, as well as the map sequence.
    
    # Map sequence of the boundaries
    fig = plt.figure(1)
    ani = charm_boundary.plot()
    
    # Middle image 193 full-disk with overlay of boundary.
    fig = plt.figure(2)
    ax = plt.subplot(projection = aia193_[int(len(aia193_)/2)])
    image = aia193_[int(len(aia193_)/2)].plot()
    aia193_[int(len(aia193_)/2)].draw_grid()
    plt.title('')
    boundary_labels = ['CHARM Boundary']
    boundary_charm = plt.contour(charm_boundary[int(len(aia193_)/2)]._data, 
                                 levels = [0.5,1.5], alpha = 0.85,
                                 color = 'black', cmap = cm.cool, linewidths = 0.5)

    boundary_charm.collections[0].set_label(boundary_labels[0])
    
    
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
    
    make_charm_boundary(time)
    