# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:03:12 2022

@author : Landon Wells

Purpose : To call on a magnetogram potential field extrapolation program named
          `Solarbextrapolation`.
"""
# General Module Imports
from astropy import units as u
from sunpy.net import attrs as a
from sunpy.net import Fido
from aiapy.calibrate import register
import sunpy.map as mp
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import matplotlib
matplotlib.use('Qt5Agg')

# Local Module Imports
from QOL_Programs.Choose_GUI_togetdata import getting_info
from SolarEventClass.HMI_Potential_Extrapolation import HMIMagnetogramExtrapolation

def plot_potential_extrapolaiton_of_hmi_data(hmi_magnetogram, map_base, zshape, zrange, resample_divide):

    map_base_0 = map_base[0]
    hmi_magnetogram_0 = hmi_magnetogram[0]

    # Make the user choose a region over which the analysis will be ran. This
    # helps to limit the size of the data to be extrapolated.
    
    print('\nSelect an area of interest to build the extrapolation.\n')
    fig = plt.figure()
    global ax
    ax = fig.add_subplot(projection = hmi_magnetogram[0])
    ani = hmi_magnetogram.plot() # DONT CHANGE ###### SUPER IMPORTANT TO STORE THE ANIMATION IN A VARIABLE OTHERWISE IT WILL GET GARBAGE COLLECTED
    global num
    num = 0
    plot_function = interactive_plot_regionselect(fig, ax, sunpymap = hmi_magnetogram)
    
    print('\nPress any button to continue.\n')
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            # In order to limit large regions of interest, the area of the chosen region is considered.
            if abs(xf-xb)*abs(yf-yb) >= 1000000 :
                areyousure = input('\nThe chosen area is LARGE.\nExtrapolations over regions this large is going to save a LOT of data.\nAre you sure you want to build a time series this large?\n(Y,N)')
                if areyousure == 'Y' or areyousure == 'y' : 
                    plt.close()
                    break
                if areyousure != 'Y' or areyousure != 'y':
                    plt.waitforbuttonpress()
            else:
                plt.close()
                break

    ref_pix = (hmi_magnetogram_0.reference_pixel)*u.pix

    xb_as, xf_as, yb_as, yf_as = (xb-ref_pix[0]/u.pix)*0.6*u.arcsec, (xf-ref_pix[0]/u.pix)*0.6*u.arcsec, (yb-ref_pix[1]/u.pix)*0.6*u.arcsec, (yf-ref_pix[1]/u.pix)*0.6*u.arcsec
    bottom_left = SkyCoord(xb_as, yb_as, frame = hmi_magnetogram_0.coordinate_frame)
    top_right = SkyCoord(xf_as, yf_as, frame = hmi_magnetogram_0.coordinate_frame)
    
    print('\nThe area of interest has been chosen.\nThe coordinates in Arc-Seconds are ',
          '(bottom left x, bottom left y) (top right x, top right y):',
          '\n(' + str(xb_as/u.arcsec) + ', ' + str(yb_as/u.arcsec) + ') (',
          str(xf_as/u.arcsec) + ', ' + str(yf_as/u.arcsec) + ').\n')
    
    # Crop the HMI and map_base maps around the area of interest.
    hmi_submap = hmi_magnetogram_0.submap(bottom_left=bottom_left, top_right=top_right)
    map_base_submap = map_base_0.submap(bottom_left=bottom_left, top_right=top_right)
    
    xrange = u.Quantity([xb_as, xf_as])
    yrange = u.Quantity([yb_as, yf_as])
    
    xrange_resample = ((abs(xb_as - xf_as))/resample_divide)*(u.pix/u.arcsec)
    yrange_resample = ((abs(yb_as - yf_as))/resample_divide)*(u.pix/u.arcsec)
    print('\nResampling Ranges: (x,y)','x= ', xrange_resample, 'area pixels', 'y= ', yrange_resample, 'area pixels.\n')
    hmi_submap_resampled = hmi_submap.resample(u.Quantity([xrange_resample, yrange_resample]), method='linear')
    hmi_submap_resampled = mp.Map(hmi_submap_resampled, sequence = True)
    
    # Now we have to call over HMI_Potential Extrapolation. This will complete the program.
    mag_extrap = HMIMagnetogramExtrapolation(hmi_magnetogram, 
                                             classification_method = 'Solarbextrapolation',
                                             hmi_submap = hmi_submap_resampled,
                                             map_base = map_base_submap, zshape = zshape,
                                             zrange = zrange, xrange = xrange, yrange = yrange,
                                             resample_divide = resample_divide)




def interactive_plot_regionselect(fig, ax, sunpymap):
    
    def onclick(event):
        global x, y, num, num_tracker
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
            rect = patches.Rectangle((xb,yb), xf-xb, yf-yb, fill = True, color = 'k', alpha = 0.5)
            ax.add_patch(rect)
        if num == 2 and num_tracker == False:
            num = 0
            num_tracker = False
            rect.remove()
    fig.canvas.mpl_connect('button_press_event', onclick)



if __name__ == '__main__':
    
    print('\nPlease input information. This information \nwill be used to build the base of the \nmap once the extrapolation is finished.\n')
    # First we ask the user to input/search for data . . . 
    info = getting_info()
    # Obtain information needed to query data
    instrument = info.telescope
    wavelength = info.wavelength
    starttime = info.starttime
    endtime = info.endtime
    time = a.Time(starttime, endtime)
    
    # For AIA maps
    if isinstance(wavelength, int):
        result = Fido.search(time,
                                 a.Instrument('AIA'),
                                 a.Wavelength(wavelength*u.angstrom),
                                 a.Sample(1800 * u.s))
        # Save the results to fits files.
        result_map = Fido.fetch(result, site='ROB')
        # Use the path to make a sunpy map.
        result_map = mp.Map(result_map[int(len(result_map)/2)], sequence = True)
        # Make the maps into level 1.5 images if they are AIA type
        result_base_map = register(result_map[0])
        result_base_map = mp.Map(result_base_map, sequence = True)
        
    # For HMI maps
    if isinstance(wavelength, str):
        result = Fido.search(time, 
                             a.Instrument('HMI'), 
                             a.Physobs(wavelength), 
                             a.Sample(1800 * u.s))
        result_map = Fido.fetch(result, site='ROB')
        # Use the path to make a sunpy map.
        result_base_map = mp.Map(result_map[int(len(result_map)/2)], sequence = True)


    # Query the contemperanious HMI data. We will need to do so because this is
    # the map that we will use to run the extrapolation method.
    result_hmi = Fido.search(time, a.Instrument('HMI'), 
                             a.Physobs('LOS_magnetic_field'), 
                             a.Sample(1800 * u.s))
    
    # Save the results to fits files.
    data_hmi = Fido.fetch(result_hmi[int(len(result_map)/2)])
    hmi_ = mp.Map(data_hmi[int(len(result_map)/2)], sequence = True)
    # Reproject the HMI data to match that of the base map.
    
    hmi_reproj_list = []
    print('Registering and scaling the HMI maps. This might take some time.')
    for i, hmi_ in enumerate(hmi_):
        
        hmi = hmi_.reproject_to(result_base_map[i].wcs)
        hmi_reproj_list.append(hmi)
    
    hmi = mp.Map(hmi_reproj_list, sequence = True)
    
    
    # This variable controls the size of the Z interpolation in units of pixels
    z_shape_size = 20
    # This variable controls the height of the Z axis used in the interpolation.
    zrange = u.Quantity([0*u.arcsec,150*u.arcsec])
    # This variable controls the size of resample of the HMI map . . . 
    # the higher the number the larger the resample size, and thus the less reliable
    # the extrapolation becomes. Vise versa, a small number will take much longer
    # to preform the extrapolation.
    resample_divide = 5
    
    plot_potential_extrapolaiton_of_hmi_data(hmi_magnetogram = hmi, 
                                             map_base = result_base_map, 
                                             zshape = z_shape_size,
                                             zrange = zrange, 
                                             resample_divide = resample_divide)