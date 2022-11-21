# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 21:39:08 2022

@author : Landon Wells

Purpose : To identify coronal loops using the TRACE method.

This program is based on the following link: 
    
https://docs.sunpy.org/projects/sunkit-image/en/latest/generated/gallery/tracing_loops.html
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
import sunkit_image.trace as trace
import numpy as np
#import warnings
#warnings.filterwarnings('ignore')

# Module Imports
from QOL_Programs.Choose_GUI_togetdata import getting_info
from SolarEventClass.ActiveRegionCoronalLoops import ActiveRegionCoronalLoops

def coronal_loop_tracing(time, *args):
    """
    Implements the Oriented Coronal CUrved Loop Tracing (OCCULT-2) algorithm
    for loop tracing in images.

    Parameters
    ----------
    time : astropy Time unit.

    args are passed onto ActiveRegionCoronalLoops

    Returns
    -------
    No returns, data is saved using ActiveRegionCoronalLoops.py

    References (I just took this from the trace.py program)
    ----------
    * Markus J. Aschwanden, Bart De Pontieu, Eugene A. Katrukha.
      Optimization of Curvi-Linear Tracing Applied to Solar Physics and Biophysics.
      Entropy, vol. 15, issue 8, pp. 3007-3030
      https://doi.org/10.3390/e15083007
    """
    
    # Figure settings (to make them look nice)
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['figure.figsize'] = [12.0,12.0]
    plt.rcParams['figure.dpi'] = 80
    plt.rcParams['savefig.dpi'] = 80
    
    
    # There is currently not an instrument specific class for the TRACE instrument.
    result_aia = Fido.search(time,
                         a.Instrument('AIA'), # Atmospheric Imaging Assembly.
                         a.Wavelength(193*u.angstrom),  # Physical observables
                         a.Sample(1800 * u.s) # Only take a shot every $var seconds.
                         ) #More observables at http://sdac.virtualsolar.org/cgi/show_details?keyword=PHYSOBS
    # Save the results to fits files, keeping in mind the path of the data.
    aia_name = Fido.fetch(result_aia, site='ROB')
    # Use the path to make a sunpy map.
    aia_ = mp.Map(aia_name, sequence = True)
    
    # Align and scale the aia171 level 1.5 image.
    # Enumerate through the list of maps . . .
    aia_list = []
    for i, aia_map in enumerate(aia_):
        aia_list.append(register(aia_map))
    aia = mp.Map(aia_list, sequence = True)
    
    # Now we call ActiveRegionCoronalLoops to draw the loops. The map containing the
    # pixel locations of the loops will be saved accordingly in Data/Classification_Data/ . . . 
    # We can also change other parameters of the TRACE loop map by altering the **kwargs given.
    # E.g. 
    # trace_loops = ActiveRegionCoronalLoops(aia, classification_method='CoronalLoops', 
    #                                        nsm1 = 3, rmin = 30, lmin = 25, nstruc = 1000, 
    #                                        ngap = 0, qthresh1 = 0.0, qthresh2 = 3.0, *args)
    trace_loops = ActiveRegionCoronalLoops(aia, classification_method = 'TRACECoronalLoops', *args)
    loop_map = trace_loops.trace_map
    
    # Map sequence of the loops
    fig = plt.figure(1)
    ani = loop_map.plot()
    
    # Middle image full-disk with overlay of loops.
    fig = plt.figure(2)
    ax = plt.subplot(projection = aia[int(len(aia)/2)])
    image = aia[int(len(aia)/2)].plot()
    aia[int(len(aia)/2)].draw_grid()
    plt.title('')
    boundary_labels = ['TRACE Loops']
    trace_loop_contour = plt.contour(loop_map[int(len(aia)/2)]._data, 
                                   levels = [0.5,1.5], alpha = 0.85,
                                   color = 'white', cmap = cm.winter, linewidths = 0.5)
    trace_loop_contour.collections[0].set_label(boundary_labels[0])
    
    
    
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
    
    # Query the data using sunpy Fido. This query is purely for the sake of
    # convience. The program only needs to know the time for the data of interest
    # (for making the full-disk CHIMERA maps).
    #result = Fido.search(time, a.Instrument(instrument), 
    #                     a.Wavelength(wavelength*u.angstrom), 
    #                     a.Sample(3600 * u.s))
    
    # Save the results to fits files.
    #query = Fido.fetch(result, site = 'ROB')
    #query = mp.Map(query[0])
    
    coronal_loop_tracing(time)
    