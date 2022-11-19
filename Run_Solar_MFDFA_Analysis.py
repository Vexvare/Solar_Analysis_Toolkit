# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 12:10:23 2022

@author : Landon Wells

Purpose : To observe multifractal 2-d image data and conduct a proper analysis 
          over the previously mentioned data.
          
          Steps of the analysis : 
              1] Observe the un-scaled alpha = 2 map, this is used to gain an idea
                 of a scaling range which might fit the map according to the image data.
              2] Using the scaling range made in the previous step, we scale the 
                 multifractal map to those scales and gain an idea of the evolution
                 of the maps degree of multifractality across the given scales. We
                 then use the multifractal evolution map to pick an area of interest.
                 This area will average the fluctuations across scales to be used in 
                 the next step.
              3] The averaged fluctuations of the previously chosen pixels across all 
                 scales will be plotted in loglog plots of f(q) vs scales. The user 
                 is then prompted to choose a scaling range that best fits the data.
              4] The previously chosen scaling range is used to calculate the h(q) vs. q
                 plots. A polynomial of order - 1 is fit across the fluctuation functions
                 for each value of q, giving the value of h(q) for each moment q. 
              
Inputs : The user is prompted using a simple GUI fucntion to point to the folder 
         where the MFDFA classification has been ran. If there is no such folder
         (found within home/user/SolarToolkitData/Classification_Data/SolarMFDFA/wavelength/starttime_---_endtime_---/lowerleft_---_upperright_---/ . . . )
         then it is required to run the analysis using the Run_Solar_MFDFA.py program.
         
         You will also need to point in a second GUI to the time-series 
         that was used to create the MFDFA data is stored.
         This is used to gather background that is relevant to the data of interest,
         such as cadence and solar events in the region of interest.
         
Outputs : No specific output is given, but the loglog plots of f(q) vs scales and
          h(q) vs q is saved in the following folder:
          home/user/SolarToolkitData/Analysis_Data/Run_Solar_MFDFA_Analysis/wavelength/starttime_---_endtime_---/fullmap_lowerleft_---_upperright_---/submap_lowerleft_---_upperright_---/

"""

# General Module Imports
import numpy as np
from numpy.polynomial.polynomial import polyfit
import sunpy.map as mp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.animation as anim
import matplotlib.cm as cm
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.patches as patches
import matplotlib.animation as animation
from astropy import units as u
from astropy.coordinates import SkyCoord
from pathlib import Path
import warnings
from matplotlib.offsetbox import AnchoredText
warnings.filterwarnings("ignore")

# Local Module imports
from SolarEventClass import Unspecified_Event_Classification
from QOL_Programs.Build_Data_SavePath import buildsavepath_submapanalysisdata
from SolarEventClass import Event_Matcher
from SolarEventClass.Solar_Events_Class import SolarEvent

plt.rcParams['font.family'] = 'serif'

def Solar_MFDFA_MapScaling_Analysis(submap_sequence, mfdfa_map, qrange, scales):



    # First we need to use the alpha = 2 map to identify a good scaling range
    # for the MFDFA map.
    # (For those unfamiliar with why this particular value is of interest, 
    # it relates to the generalized Hurst exponent and from it much information
    # about the time series can be obtained (TIME TO READ UP ON MULTIFRACTALITY 
    # OF TIME SERIES)).
    alpha_2_index = np.where(( qrange == 2 ))[0][0]
    alpha_map_ = mfdfa_map[:,:, :, alpha_2_index]

    # We also want to identify the temporal difference between the images (cadence)
    image_cadence = abs(int(submap_sequence[0].meta.get('date-obs')[17:19]) - int(submap_sequence[1].meta.get('date-obs')[17:19]))
    
    # We first initialize the alpha map by choosing the fitting range across
    # all scales to gain an idea of where we should choose a scaling range.
    print('\nThe DFA map will have coordinates with respect to the middle image of the time series.\n')
    alpha_map = scaling_dfa_analysis(alpha_map = alpha_map_, scales = scales, order = 1)
    sunpy_alpha_map = mp.Map([alpha_map], submap_sequence[int(len(submap_sequence)/2)+1].meta)
    
    # Plotting a time series of the cut map of interest.
    global ax
    fig = plt.figure(1)
    ax = fig.add_subplot(projection = submap_sequence[0])
    anix = submap_sequence.plot()
    
    # Test to see if any event matches the time series of interest.
    eventmatchs = Event_Matcher.solar_event_matcher(sunpymap = submap_sequence)
    
    # If events match the timeseries, we can use those events and overlay them 
    # onto the figures.
    # We can call over a plotting method for specific events. We just need to pass
    # over the eventsmatchs variable containing all of the matched events.
    if eventmatchs != []:
        SolarEvent(submap_sequence).plot_matched_events(eventlist = eventmatchs, 
                                                        fig = fig, ax = ax)
    
    # Plotting the initial figure with the devices we will use to adjust the figure
    fig = plt.figure(2)
    ax = fig.add_subplot(projection = sunpy_alpha_map)
    if eventmatchs != []:
        SolarEvent(submap_sequence).plot_matched_events(eventlist = eventmatchs, 
                                                        fig = fig, ax = ax)
    ani = sunpy_alpha_map.plot(norm = matplotlib.colors.Normalize(vmin = 0.38, 
                                                                  vmax = np.max(alpha_map)), 
                               cmap = cm.viridis)
    plt.set_cmap("viridis")
    
    fig.subplots_adjust(bottom = 0.25)
    axscales_start = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    axscales_end = fig.add_axes([0.25, 0.05, 0.65, 0.03])
    

    # Calling the interactive_plot
    startscale_slider  = Slider(ax = axscales_start, 
                                label = 'Range Start (mins)', 
                                valmin = scales[0]*image_cadence/60, 
                                valmax = scales[-1]*image_cadence/60,
                                valinit = scales[0]*image_cadence/60,
                                valstep = scales*image_cadence/60,
                                facecolor = 'blue',
                                track_color = 'black')
    
    endscale_slider  = Slider(ax = axscales_end, 
                              label = 'Range End (mins)', 
                              valmin = scales[0]*image_cadence/60, 
                              valmax = scales[-1]*image_cadence/60,
                              valinit = scales[-1]*image_cadence/60,
                              valstep = scales*image_cadence/60,
                              facecolor = 'blue',
                              track_color = 'black')
        
    def update_start(val):
        scalestart = int(np.where((scales == int(startscale_slider.val*60/image_cadence)))[0][0])
        alpha_map = scaling_dfa_analysis(alpha_map = alpha_map_, scales = scales, 
                                         order = 1, startscale = scalestart)
        sunpy_alpha_map = mp.Map([alpha_map], submap_sequence[int(len(submap_sequence)/2)+1].meta)
        global ax
        ax.clear()
        ax = fig.add_subplot(projection = sunpy_alpha_map)
        if eventmatchs != []:
            SolarEvent(submap_sequence).plot_matched_events(eventlist = eventmatchs, 
                                                            fig = fig, ax = ax)
        global ani
        ani = sunpy_alpha_map.plot(norm = matplotlib.colors.Normalize(vmin = 0.38, 
                                                                      vmax = np.max(alpha_map)), 
                                   cmap = cm.viridis)
        fig.canvas.draw_idle()
        
    def update_end(val):
        scaleend = int(np.where((scales == int(endscale_slider.val*60/image_cadence)))[0][0])
        alpha_map = scaling_dfa_analysis(alpha_map = alpha_map_, scales = scales, order = 1, endscale = scaleend)
        sunpy_alpha_map = mp.Map([alpha_map], submap_sequence[int(len(submap_sequence)/2)+1].meta)
        global ax
        ax.clear()
        ax = fig.add_subplot(projection = sunpy_alpha_map)
        if eventmatchs != []:
            SolarEvent(submap_sequence).plot_matched_events(eventlist = eventmatchs, 
                                                            fig = fig, ax = ax)
        global ani
        ani = sunpy_alpha_map.plot(norm = matplotlib.colors.Normalize(vmin = 0.38, 
                                                                      vmax = np.max(alpha_map)), 
                                   cmap = cm.viridis)
        fig.canvas.draw_idle()
        
    print('\nPress any button to continue once you have identified your approximate scaling range.\n')
    plt.waitforbuttonpress()
    while True:
        anim.FuncAnimation(fig, startscale_slider.on_changed(update_start))
        anim.FuncAnimation(fig, endscale_slider.on_changed(update_end))
        plt.show()
        if plt.waitforbuttonpress():
            plt.close()
            break

    # Now we plot the new scaling range across the MFDFA evolution map. We will 
    # use this map to identify an area of interest where the multifractality
    # evolves the most (or whatever other features might catch our interest).
    print('\nWe will use the chosen scaling range of : [' + str(startscale_slider.val) + '-' + str(endscale_slider.val) + '] minutes to scale the next image.\n')
    scalestart = int(np.where((scales == int(startscale_slider.val*60/image_cadence)))[0][0])
    scaleend = int(np.where((scales == int(endscale_slider.val*60/image_cadence)))[0][0])
    scaled_maps = scaling_mfdfa_analysis(mfdfa_map, 
                                         scales, 
                                         qrange, 
                                         order = 1, 
                                         startscale = scalestart,
                                         endscale = scaleend)
    
    print('\nThe MFDFA map will have coordinates with respect to the middle image of the time series.\n')
    map_list = []
    for qcount in range(len(scaled_maps[0,0,:])):
        map_list.append(mp.Map([scaled_maps[:,:,qcount], submap_sequence[int(len(submap_sequence)/2)+1].meta]))
    mfdfa_map_series = mp.Map(map_list, sequence = True)
    
    
    fig = plt.figure(3)
    ax = fig.add_subplot(projection = mfdfa_map_series[0])
    legend_elements = '[' + str(startscale_slider.val) + '-' + str(endscale_slider.val) + '] minutes'
    at = AnchoredText(legend_elements, prop = dict(size=10), frameon = True, loc = 'upper right')
    at.patch.set_boxstyle('round,pad=0,rounding_size=0.2')
    ax.add_artist(at)
    ani = mfdfa_map_series.plot(norm = matplotlib.colors.Normalize(vmin = 0.38, 
                         vmax = np.max(scaled_maps)), 
                         cmap = cm.viridis)
    plt.set_cmap("viridis")
    plt.colorbar()
    global num
    num = 0
    plot_function = interactive_plot_regionselect(fig, ax, sunpymap = mfdfa_map_series)
    print('\nPress any button to continue once you have identified an area of interest.\n')
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            break

    # Now that we have identified an area of interest, we will observe the true 
    # scaling and fit the best scaling range across the pixels of interest.
    
    # First we build the coordinates of the chosen region.
    global xb, xf, yb, yf
    bot_left = submap_sequence[int(len(submap_sequence)/2)+1].bottom_left_coord
    top_right = submap_sequence[int(len(submap_sequence)/2)+1].bottom_left_coord
    
    xb_as = bot_left.Tx
    yb_as = bot_left.Ty
    xf_as = top_right.Tx
    yf_as = top_right.Ty
    
    # These are the chosen coords in arc-seconds with respect to the (0,0) of the map
    as_xb = submap_sequence[0].meta.get('cdelt1')*(xb*u.arcsec)
    as_yb = submap_sequence[0].meta.get('cdelt2')*(yb*u.arcsec)
    as_xf = submap_sequence[-1].meta.get('cdelt1')*(xf*u.arcsec)
    as_yf = submap_sequence[-1].meta.get('cdelt2')*(yf*u.arcsec)
    
    chosen_xb = xb_as + as_xb 
    chosen_yb = yb_as + as_yb
    chosen_xf = xf_as + as_xf
    chosen_yf = yf_as + as_yf
    
    print('\nThe area of interest has been chosen.\nThe coordinates in Arc-Seconds are ',
          '(bottom left x, bottom left y) (top right x, top right y):',
          '\n(' + str(chosen_xb/u.arcsec) + ', ' + str(chosen_yb/u.arcsec) + ') (',
          str(chosen_xf/u.arcsec) + ', ' + str(chosen_yf/u.arcsec) + ').\n',
          'These coordinates are with respect to the middle image`s coordinates.\n')
    
    bottom_left = SkyCoord(xb_as, yb_as, frame = submap_sequence[int(len(submap_sequence)/2)+1].coordinate_frame)
    top_right = SkyCoord(xf_as, yf_as, frame = submap_sequence[int(len(submap_sequence)/2)+1].coordinate_frame)
    
    # # This next map is purely to help build the save path.
    subsubmap_sequence = mp.Map([submap_sequence[int(len(submap_sequence)/2)+1].submap(bottom_left = bottom_left,
                                                                                       top_right = top_right)], sequence = True)
    
    # Building the save path.
    true_path = buildsavepath_submapanalysisdata(sunpymap_sequence = submap_sequence, 
                                                 sunpysubmap_sequence = subsubmap_sequence, 
                                                 analysis_method = 'Run_Solar_MFDFA_Analysis')
    
    
    # Close and save the previous plot (We choose to close it here because we need the code above to save it correctly).
    writervideo = animation.FFMpegWriter(fps=2)
    ani.save(true_path + 'SubReg.mp4', writer = writervideo)
    plt.close()
    
    # Cutting the MFDFA_MAP across the pixels that were chosen in the previous plot.    
    # + 1 because in the next line, mfdfa_map_cut doesn't include yf or xf rows/columns
    # Also, we can't use .as_array() because it`s a 4-D map data
    mfdfa_map_cut = mfdfa_map[int(yb):int(yf)+1,int(xb):int(xf)+1,:,:]
    
    # Quick save to the fluctuation functions, should one need to call back on them.
    np.save(true_path + 'mfdfa_map_cut.npy', mfdfa_map_cut)
    # While we're at it, lets save the qvalues and all of the scales, should we need to call back on them.
    np.save(true_path + 'qrange.npy', qrange)
    np.save(true_path + 'scales.npy', scales)
    
    # Averaging the fluctuation functions across scales for all pixels
    mfdfa_cut_averaged = np.mean((np.mean(mfdfa_map_cut, axis = 0)), axis = 0)
    # Quick save to the averaged fluctuation functions, should one need to call back on them.
    np.save(true_path + 'mfdfa_cut_averaged.npy', mfdfa_cut_averaged)

    
    # Creating log log plots of the averaged fluctuation 
    # function vs scales for the area of interest
    print('\nPick the scaling range that best fits the averaged fluctuation functions.')
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot()
    plt.loglog(scales*image_cadence/60, mfdfa_cut_averaged, '.b' , markersize = 0.3, picker = point_picker)
    #fluct_2q = plt.errorbar(scales*image_cadence/60, mfdfa_cut_averaged[:,alpha_2_index], 
    #                      yerr = mfdfa_cut_err[:,alpha_2_index], errorevery = 20, ecolor = 'black')
    #fluct_minq = plt.errorbar(scales*image_cadence/60, mfdfa_cut_averaged[:,0], 
    #                       yerr = mfdfa_cut_err[:,0], errorevery = 20, ecolor = 'black')
    #fluct_maxq = plt.errorbar(scales*image_cadence/60, mfdfa_cut_averaged[:,-1], 
    #                       yerr = mfdfa_cut_err[:,-1], errorevery = 20, ecolor = 'black')
    fluct_2q = plt.loglog(scales*image_cadence/60, mfdfa_cut_averaged[:,alpha_2_index])
    fluct_minq = plt.loglog(scales*image_cadence/60, mfdfa_cut_averaged[:,0])
    fluct_maxq = plt.loglog(scales*image_cadence/60, mfdfa_cut_averaged[:,-1])
    
    ax.set_xticks(np.arange(10,scales[-1]*image_cadence/60, step = 20).astype(int))
    ax.set_xticklabels(np.arange(10,scales[-1]*image_cadence/60, step = 20).astype(int), fontsize = 15, weight = 'semibold')
    plt.xlabel('$log_{10}[time]$ (minutes)', fontsize = 15, labelpad = 5)
    plt.yticks(fontsize = 15, weight = 'semibold')
    plt.ylabel('$log_{10}[F_{q}]$', fontsize = 15, labelpad = 5)
    plt.legend(handles = [fluct_minq, fluct_2q, fluct_maxq], 
               labels = ['q = ' + str(qrange[0]), 'q = ' + str(qrange[alpha_2_index]), 'q = ' + str(qrange[-1])], 
               fontsize = 15)
    fig.canvas.mpl_connect('pick_event', pick_scalingrange)
    plt.waitforbuttonpress()
    
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            plt.savefig(fname = true_path + 'Fluct.png' , format = 'png')
            plt.savefig(fname = true_path + 'Fluct.jpg' , format = 'jpg')
            plt.close()
            break
    
    global scaling_list
    # Assumes that the first and last picked events. (Could be optimized by
    # looking for the largest and smallest scales in the list, but im lazy).
    print('\nThe scaling range for the next plot is from : [' + str(scaling_list[0]) + '-' + str(scaling_list[-1]) + '] minutes.')
    
    beginning_scale = int(scaling_list[0]*60/image_cadence)
    ending_scale = int(scaling_list[-1]*60/image_cadence)
    
    # Quick save to the new scaling range.
    np.save(true_path + 'chosen_scales.npy', scales[beginning_scale:ending_scale])
    
    # Now we observe the degree of multifractality by looking at h(q) vs q.
    mfdfa_hq = polyfit(np.log(scales[beginning_scale:ending_scale]),
                       np.log(mfdfa_cut_averaged[beginning_scale:ending_scale,:]), 1, full = True)[0][1]
    mfdfa_hq_rounded = np.round_(mfdfa_hq, decimals = 2)
    
    mfdfa_hq_error = np.sqrt((polyfit(np.log(scales[beginning_scale:ending_scale]),
                       np.log(mfdfa_cut_averaged[beginning_scale:ending_scale,:]), 
                       1, full = True)[1][0]))/np.sqrt(len(scales[beginning_scale:ending_scale]))
    
    rounded_mfdfa_hq_error = np.round_(mfdfa_hq_error, decimals = 2)
    
    
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot()
    for i, q in enumerate(qrange):
        if q == qrange[0] or q == qrange[alpha_2_index] or q == qrange[-1]:
            plt.errorbar(q, mfdfa_hq[i], 
                         yerr = rounded_mfdfa_hq_error[i], 
                         ecolor = 'black', marker = '.', 
                         markersize = 10, color = 'red')
        else:
            plt.errorbar(q, mfdfa_hq[i], 
                         yerr = rounded_mfdfa_hq_error[i], 
                         ecolor = 'black', marker = '.',
                         markersize = 10, color = 'blue')
    
    ax.set_xticks(np.arange(qrange[0],qrange[-1]+1, step = 2).astype(int))
    ax.set_xticklabels(np.arange(qrange[0],qrange[-1]+1, step = 2).astype(int), fontsize = 15, weight = 'semibold')
    plt.xlabel('q', fontsize = 15, weight = 'semibold', labelpad = 5)
    ax.set_ylim([0,2])
    plt.yticks(fontsize = 15, weight = 'semibold') 
    plt.ylabel('h(q)', fontsize = 15, weight = 'semibold', labelpad = 5)
    legend_elements = 'h('+ str(int(qrange[0])) +') = ' + str(mfdfa_hq_rounded[0]) + ' +/- ' + str(rounded_mfdfa_hq_error[0]) + '\nh('+ str(int(qrange[alpha_2_index])) +') = ' + str(mfdfa_hq_rounded[alpha_2_index]) + ' +/- ' + str(rounded_mfdfa_hq_error[alpha_2_index]) + '\nh('+ str(int(qrange[-1])) +') = ' + str(mfdfa_hq_rounded[-1]) + ' +/- ' + str(rounded_mfdfa_hq_error[-1]) + '\nScaling Range : ' + str(scaling_list[0]) + '-' + str(scaling_list[-1]) + ' [mins]'
    print(legend_elements)
    at = AnchoredText(legend_elements, prop = dict(size=15), frameon = True, loc = 'upper right')
    at.patch.set_boxstyle('round,pad=0,rounding_size=0.2')
    ax.add_artist(at)
    
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            plt.savefig(fname = true_path + 'h_q_vs_q.png' , format = 'png')
            plt.savefig(fname = true_path + 'h_q_vs_q.jpg' , format = 'jpg')
            plt.close()
            plt.close(1)
            break

    print('All done! Refer back to the previously mentioned save path to revisit the analysis.')





global scaling_list, scales
scaling_list = []
    

def pick_scalingrange(event):
    global scaling_list
    scaling_list.append(event.pickx)
    
    
def point_picker(point, mouseevent):
    '''
    Finds the x data points within a certain distance from the mouseclick in 
    data coords and attaches some extra attributes, pickx which are the data
    points that were picked.
    '''

    if mouseevent.xdata is None:
        return False, dict()
    xdata = point.get_xdata()
    max_dist = 0.01
    d = xdata-mouseevent.xdata
    
    ind = np.nonzero(d <= max_dist)

    if len(ind):
        pickx = xdata[ind][-2]
        print(pickx)
        props = dict(ind = ind, pickx = pickx)
        return True, props
    else:
        return False, dict()
    
    
    
    
    
def scaling_dfa_analysis(alpha_map, scales, order = 1, startscale = 0, endscale = -1):
    '''
    Function to scale the dfa map (q = 2) generalized hurst exponent map.
    '''
    alpha_map_scaled = np.zeros([len(alpha_map[:,0,0]), len(alpha_map[0,:,0])])
    for xc in range (len(alpha_map[0,:,0])):
        for yc in range (len(alpha_map[:,0,0])):
            alpha_map_scaled[yc,xc] = polyfit(np.log(scales[startscale:endscale]), 
                                       np.log(alpha_map[yc,xc,startscale:endscale]), 
                                       order)[order]
    return alpha_map_scaled

def scaling_mfdfa_analysis(mfdfa_map, scales, qrange, order = 1, startscale = 0, endscale = -1):
    '''
    Function to scale all of the fluctuation function maps across all moments (values of q).
    '''
    mfdfa_map_scaled = np.zeros([len(mfdfa_map[:,0,0,0]), len(mfdfa_map[0,:,0,0]), len(qrange)])
    
    for xc in range (len(mfdfa_map[0,:,0,0])):
        for yc in range (len(mfdfa_map[:,0,0,0])):
            for qc in range (len(mfdfa_map[0,0,0,:])):
                mfdfa_map_scaled[yc,xc,qc] = polyfit(np.log(scales[startscale:endscale]), 
                                                     np.log(mfdfa_map[yc,xc,startscale:endscale,qc]), 
                                                     order)[order]
    return mfdfa_map_scaled

def interactive_plot_regionselect(fig, ax, sunpymap):
    '''
    Function to interact with the plot of interest and choose a rectangular region.
    '''
    def onclick(event):
        global x, y, num, num_tracker
        x = event.xdata
        y = event.ydata
        num_tracker = False
        if num == 0 and num_tracker == False:
            global xb, yb
            xb = int(x + 0.5) - 0.5
            yb = int(y + 0.5) - 0.5
            num += 1
            num_tracker = True
        if num == 1 and num_tracker == False:
            global xf, yf, rect
            xf = int(x + 0.5) + 0.5
            yf = int(y + 0.5) + 0.5
            num += 1
            num_tracker = True
            if xf < xb :
                xtemp = xf
                xf = xb
                xb = xtemp
            if yf < yb :
                ytemp = yf
                yf = yb
                yb = ytemp
            rect = patches.Rectangle((xb,yb), xf-xb, yf-yb, fill = True, color = 'k', alpha = 0.5)
            ax.add_patch(rect)
        if num == 2 and num_tracker == False:
            num = 0
            num_tracker = False
            rect.remove()
    fig.canvas.mpl_connect('button_press_event', onclick)



if __name__ == '__main__' :
    # Obtaining the data . . . 
    submap_sequence, mfdfa_map, qrange, scales, datafilepath = Unspecified_Event_Classification.UnspecifiedClassification.load_MFDFA_map()
    
    # Conducting the analysis . . . 
    Solar_MFDFA_MapScaling_Analysis(submap_sequence, mfdfa_map, qrange, scales)