# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 14:02:21 2022

@author : Landon Wells

Purpose : Create Coronal Hole Boundaries using the CHIMERA classification
          method.
          
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
from SolarEventClass.CoronalHoleBoundaries import CoronalHoleBoundary


def make_chimera_boundary(time, *args):
    
    # Figure settings (to make them look nice)
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['figure.figsize'] = [12.0,12.0]
    plt.rcParams['figure.dpi'] = 80
    plt.rcParams['savefig.dpi'] = 80
    
    # if args is given, it is assumed to be the maps that were queried previously.
    
    # But for all intensive purposes, CHIMERA uses 171, 193, and 211 full disk 
    # AIA data. So regardless of the previous query, we use the time of interest
    # to download and search for data in those specific wavelengths. ALSO, CHIMERA 
    # in the original IDL program uses the HMI full disk image. However, due to 
    # time constraints, the PYTHON version does not use the HMI image, but we still
    # query it regardless because it might one day be incorporated.
    
    
    # Build the query via Fido, this can return one item, or a list of them to DL matching
    # the given filters. Full documentation on Fido and it's attributes can be found on the Sunpy website.
    
    # Obtain the corresponding AIA data showing the EUV images at the time of interest.
    
    # Searching for the 171, 193, and 211 data via Sunpy Fido. (These are the Coronal lines).
    result_aia = Fido.search(time,
                         a.Instrument('AIA'), # Atmospheric Imaging Assembly.
                         a.Wavelength(171*u.angstrom),  # Physical observables
                         a.Sample(1800 * u.s) # Only take a shot every $var seconds.
                         ) #More observables at http://sdac.virtualsolar.org/cgi/show_details?keyword=PHYSOBS
    # Save the results to fits files, keeping in mind the path of the data.
    aia171_name = Fido.fetch(result_aia, site='ROB')
    # Use the path to make a sunpy map.
    #aia171_ = mp.Map(aia171_name[0])
    aia171_ = mp.Map(aia171_name, sequence = True)
    # Repeat for the 193
    result_aia = Fido.search(time,
                         a.Instrument('AIA'),
                         a.Wavelength(193*u.angstrom),
                         a.Sample(1800 * u.s)) 
    aia193_name = Fido.fetch(result_aia, site='ROB')
    #aia193_ = mp.Map(aia193_name[0])
    aia193_ = mp.Map(aia193_name, sequence = True)
    # Repeat for the 211
    result_aia = Fido.search(time,
                         a.Instrument('AIA'),
                         a.Wavelength(211*u.angstrom),
                         a.Sample(1800 * u.s))
    aia211_name = Fido.fetch(result_aia, site='ROB')
    #aia211_ = mp.Map(aia211_name[0])
    aia211_ = mp.Map(aia211_name, sequence = True)
    
    # HMI
    result_hmi = Fido.search(time,
                          a.Instrument('HMI'), # Helioseismic and Magnetic Imager.
                          a.Physobs('LOS_magnetic_field'),
                          a.Sample(1800 * u.s)) 
    hmi_name = Fido.fetch(result_hmi[0])
    hmi_ = mp.Map(hmi_name, sequence = True)
    
    # Converting the level 1 AIA data into level 1.5 AIA data
    #aia171 = register(aia171_)
    #aia193 = register(aia193_)
    #aia211 = register(aia211_)
    
    # Align and scale the HMI to the aia171 level 1.5 image.
    
    # After some testing, SunPY claims that you can use the below function : 
    # https://aiapy.readthedocs.io/en/stable/api/aiapy.calibrate.register.html#aiapy.calibrate.register
    #hmi_reproj = register(hmi)
    
    # BUT . . . it seems like the below method is . . . better. 
    # Check for yourself on the full-disk HMI image and compare with the AIA maps.
    # https://docs.sunpy.org/en/stable/generated/gallery/map_transformations/reprojection_align_aia_hmi.html
    #hmi = hmi_.reproject_to(aia171_[0].wcs)

    # CHIMERA Boundary 
    # We are going to call the Solar_Events_Class and define a CH Boundary 
    # event that uses the CHIMERA boundary. Calling this class will save the 
    # chimera boundary locally, so that other programs may make use of it.
    boundary = CoronalHoleBoundary(aia171_, aia193_, aia211_, hmi_, classification_method = 'CHIMERA')
    
    # This object contains (most) of the possible maps that chimera can make. 
    # For the most part however, we're really only concerned with the first map.
    chimera_boundaries = boundary.chimera_map

    # Now we show the sequence of maps which make up the magnetogram boundaries.
    fig = plt.figure(1)
    ani = chimera_boundaries.plot()

    # Plotting full-disk LEVEL 1.5 AIA maps
    
    #    171 FULL DISK
    fig = plt.figure(2)
    ax = plt.subplot(projection = aia171_[int(len(aia171_)/2)])
    aia171_[int(len(aia171_)/2)].plot()
    aia171_[int(len(aia171_)/2)].draw_grid()
    plt.title('')
    boundary_labels = ['CHIMERA Boundary']
    boundary_chimera = plt.contour(chimera_boundaries[int(len(chimera_boundaries)/2)]._data, 
                                   levels = [0.5,1.5], alpha = 0.85,
                                   color = 'white', cmap = cm.winter, linewidths = 0.5)
    boundary_chimera.collections[0].set_label(boundary_labels[0])
    
    #    193 FULL DISK
    fig = plt.figure(3)
    ax = plt.subplot(projection = aia193_[int(len(aia193_)/2)])
    aia193_[int(len(aia193_)/2)].plot()
    aia193_[int(len(aia193_)/2)].draw_grid()
    plt.title('')
    boundary_labels = ['CHIMERA Boundary']
    boundary_chimera = plt.contour(chimera_boundaries[int(len(chimera_boundaries)/2)]._data, 
                                   levels = [0.5,1.5], alpha = 0.85,
                                   color = 'white', cmap = cm.winter, linewidths = 0.5)
    boundary_chimera.collections[0].set_label(boundary_labels[0])
    
    #    211 FULL DISK
    fig = plt.figure(4)
    ax = plt.subplot(projection = aia211_[int(len(aia211_)/2)])
    aia211_[int(len(aia211_)/2)].plot()
    aia211_[int(len(aia211_)/2)].draw_grid()
    plt.title('')
    boundary_labels = ['CHIMERA Boundary']
    boundary_chimera = plt.contour(chimera_boundaries[int(len(chimera_boundaries)/2)]._data, 
                                   levels = [0.5,1.5], alpha = 0.85,
                                   color = 'white', cmap = cm.winter, linewidths = 0.5)
    boundary_chimera.collections[0].set_label(boundary_labels[0])
    
    

    # Show all of the plots
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            plt.close(1)
            plt.close(2)
            plt.close(3)
            plt.close(4)
            break



if __name__ == '__main__':
    
    # First we ask the user to input/search for data . . . 
    info = getting_info()
    instrument = info.telescope
    wavelength = info.wavelength
    starttime = info.starttime
    endtime = info.endtime
    time = a.Time(starttime, endtime)
    
    make_chimera_boundary(time)
    
    