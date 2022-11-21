# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 13:30:10 2022

@author : Landon Wells

Purpose : To analyze solar time sequence data using MFDFA.

Input : None. However, using a QOL_Program, you will need
        to point to the file where the time series 
        data is stored. If you do not have any time-series
        yet, please run the Create_AIA_Time_Series.py program
        first to build a time series. Once you have done that,
        then the time series will be located in the following
        folder : home/user/Solar_Analysis_Toolkit/Data/ . . . / [data].
        
Outputs : A 4-D MFDFA numpy array that is saved locally to be used in the secondary
          program Run_Solar_MFDFA_Analysis.py .

"""

# General Module imports
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map as mp
import matplotlib.patches as patches
import matplotlib
matplotlib.use('Qt5Agg')

# Local Module imports
from QOL_Programs.GUI_FilePath import ObtainMapFiles
from SolarEventClass import Unspecified_Event_Classification
from SolarEventClass import Event_Matcher
from SolarEventClass.Solar_Events_Class import SolarEvent

def run_solar_mfdfa(timeseries: np.ndarray, 
                    qrange: np.ndarray = 2, 
                    orderfit: int = 1, 
                    scales: np.ndarray = None, **kwargs):

    # Test to see if any event matches the time series of interest.
    eventmatchs = Event_Matcher.solar_event_matcher(sunpymap = map_sequence) 
    
    print(eventmatchs)
    
    # Building an animation of the time-series.
    fig = plt.figure()
    ax = fig.add_subplot(projection = map_sequence[0])
    ani = map_sequence.plot()
    
    # If events match the timeseries, we can use those events and overlay them 
    # onto the figures.
    # We can call over a plotting method for specific events. We just need to pass
    # over the eventsmatchs variable containing all of the matched events.
    if eventmatchs != []:
        SolarEvent(map_sequence).plot_matched_events(eventlist = eventmatchs, fig = fig, ax = ax)

    # Allowing the user to specify an area over which the MFDFA analysis is to be ran over.
    print('\nSelect an area of which the MF-DFA analysis is to be ran over.')
    global num
    num = 0
    plot_function = interactive_plot_regionselect(fig, ax, sunpymap = ani)
    
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            # In order to limit large regions of interest, the area of the chosen region is considered.
            if abs(xf-xb)*abs(yf-yb) >= 10000 :
                areyousure = input('The chosen area is LARGE.\n The MFDFA calculation is not friendly to slow computers.\n Are you sure you want to preform the MFDFA analysis over your region? \n (Y,N)')
                if areyousure == 'Y' or areyousure == 'y' : 
                    plt.close()
                    break
                if areyousure != 'Y' or areyousure != 'y':
                    plt.waitforbuttonpress()
            else:
                plt.close()
                break
    
    bot_left = map_sequence[0].bottom_left_coord
    # We use the bottom left coord of the last map and add the top right of the 
    # chosen box to this value. This gives the chosen coordinate for the final map.
    top_right = map_sequence[-1].bottom_left_coord
    
    xb_as = bot_left.Tx
    yb_as = bot_left.Ty
    xf_as = top_right.Tx
    yf_as = top_right.Ty
    
    # These are the chosen coords in arc-seconds with respect to the (0,0) of the map
    as_xb = map_sequence[0].meta.get('cdelt1')*(xb*u.arcsec)
    as_yb = map_sequence[0].meta.get('cdelt2')*(yb*u.arcsec)
    as_xf = map_sequence[-1].meta.get('cdelt1')*(xf*u.arcsec)
    as_yf = map_sequence[-1].meta.get('cdelt2')*(yf*u.arcsec)
    
    chosen_xb = xb_as + as_xb 
    chosen_yb = yb_as + as_yb
    chosen_xf = xf_as + as_xf
    chosen_yf = yf_as + as_yf
    
    print('\nThe area of interest has been chosen.\nThe coordinates in Arc-Seconds are ',
          '(bottom left x, bottom left y) (top right x, top right y):',
          '\n(' + str(chosen_xb/u.arcsec) + ', ' + str(chosen_yb/u.arcsec) + ') (',
          str(chosen_xf/u.arcsec) + ', ' + str(chosen_yf/u.arcsec) + ').\n')
    
    # Although we just obtained the coordinates of the cut, we convert the timeseries
    # into an array to actually do the cutting (purely for the data, the headers
    # will be cut later).
    map_sequence_array = map_sequence.as_array()
    map_sequence_array_cut = map_sequence_array[int(yb):int(yf),int(xb):int(xf),:]
    
    bottom_left = SkyCoord(chosen_xb, chosen_yb, frame = map_sequence[0].coordinate_frame)
    top_right = SkyCoord(chosen_xf, chosen_yf, frame = map_sequence[-1].coordinate_frame)
    
    submap_list = []
    for i, map_ in enumerate(map_sequence):
        # Information that cuts the meta/header information of the submap
        bot_left = map_sequence[i].bottom_left_coord
        xb_as = bot_left.Tx
        yb_as = bot_left.Ty
        cut_xb = xb_as + as_xb
        cut_yb = yb_as + as_yb
        cut_xf = xb_as + as_xf
        cut_yf = yb_as + as_yf
        bottom_left = SkyCoord(cut_xb, cut_yb, frame = map_sequence[i].coordinate_frame)
        top_right = SkyCoord(cut_xf, cut_yf, frame = map_sequence[i].coordinate_frame)
        # Making a sunpymap out of the newly cut data and header information 
        # and appending it into the submap list
        submap_list.append(mp.Map([map_sequence_array_cut[:,:,i], map_.submap(bottom_left = bottom_left, top_right = top_right).meta]))
    # Making a map sequence out of the submaps
    map_sequence_submaps = mp.Map(submap_list, sequence = True)
    
    
    fig = plt.figure()
    ax = plt.subplot(projection = map_sequence_submaps[0])
    image = map_sequence_submaps.plot()
    
    # Show the chosen region.
    print('\nPress any button to continue.\n')
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            plt.close()
            break

    # Creating a second Unspecified Event, now we use the SolarMFDFA classification
    # on the newly cut map_sequence_submap and obtain the mfdfamap object. We then
    # choose to save the data, so that it may be used in future analysis.
    mfdfa_map_object = Unspecified_Event_Classification.UnspecifiedClassification(map_sequence_submaps, 
                                                                                  classification_method = 'SolarMFDFA',
                                                                                  qrange = qrange, 
                                                                                  orderfit = orderfit, 
                                                                                  scales = scales)
    
    # Now we save the data! We will use the data that we just saved in another 
    # program, 'Run_Solar_MFDFA_Analysis.py'. We do this because there are 
    # many different way's to interpret the data. We will discuss them in the
    # afformentioned program.
    mfdfa_map_object.save_MFDFA_map()
    

    
def interactive_plot_regionselect(fig, ax, sunpymap):
    
    def onclick(event):
        global go, x, y, num, num_tracker
        go = True
        x = event.xdata
        y = event.ydata
        num_tracker = False
        if num == 0 and num_tracker == False:
            global xb, yb
            xb = x
            yb = y
            num += 1
            num_tracker = True
        if num == 1 and num_tracker == False:
            global xf, yf, rect
            xf = x
            yf = y
            num += 1
            num_tracker = True
            rect = patches.Rectangle((int(xb),int(yb)), int(xf-xb), int(yf-yb), fill = True, color = 'k', alpha = 0.5)
            ax.add_patch(rect)
        if num == 2 and num_tracker == False:
            num = 0
            num_tracker = False
            rect.remove()
    fig.canvas.mpl_connect('button_press_event', onclick)
    

if __name__ == '__main__':
    
    # Obtaining the time-series using QOL programs . . . 
    info = ObtainMapFiles()
    map_sequence = info.map_sequence
    # TODO MAKE A GUI FUNCTION TO SPECIFY Q'S AND SCALES FOR MFDFA ANALYSIS
    qrange = np.linspace(-10, 10, num = 21, endpoint = True)
    
    run_solar_mfdfa(map_sequence, qrange)