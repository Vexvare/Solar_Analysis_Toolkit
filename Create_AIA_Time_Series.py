# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 12:45:32 2022

@author : Landon Wells


Purpose : This program is to create a Time Series as a stack of images, of particular 
          interest is the AIA EUV wavelenghts, 171, 193, 211, 304 A ect. We will be 
          obtaining the time series from the JSOC repository. We will be using SunPY 
          modules to search for the time series from JSOC directly, as Fido does not 
          'seperate the staging and downloading steps'. This is a 3 stage process, with 
          which we should end up with a MapSequence that has been co-aligned wrt the 
          reference image (the zero'th image in the stack), and can be called back 
          for futher analysis after the data has been saved.

Inputs :  The user is required to input the following parameters : 
          start_time - The time that the data of interest starts. I recommend 
                       visiting sites like https://www.helioviewer.org/ to gain 
                       an idea of a good start time.
          end_time - The time that the data of interest ends.
          wavelength - Because this is AIA data that we are making a time series of,
                       we will need to specifiy a wavelength of interest.
          jsoc_email - A registered email address from JSOC.

Outputs : Multiple Sunpy Maps. These maps are saved in the following directory : 
          .\SolarToolkitData\Solar_Data\wavelength\starttime_---_endtime_---\lowerleft_---_upperright_---\ . . . '
          These maps can be used to call programs that use time series as data inputs, 
          such as Run_Solar_MFDFA.py.


NOTE : You will need an email registered with JSOC in order to query and 
download data.

This most recent version is based off of the following webpage
https://docs.sunpy.org/en/stable/generated/gallery/acquiring_data/downloading_cutouts.html
and the following link will be used to coalign the map sequence once it is
queried and downloaded from the JSOC repository. 
https://docs.sunpy.org/projects/sunkit-image/en/latest/generated/gallery/mapsequence_coalignment.html#sphx-glr-generated-gallery-mapsequence-coalignment-py

TIP : RUN THE PROGRAM TO COMPLETION. SOMETIMES THERE MAY BE ONE OR TWO MAPS WHICH CAN`T BE READ IN CORRECTLY.
    IN THESE CASES, DELETE THE MAP THAT CANT BE READ IN CORRECTLY IN THE SUNPY DATA DIRECTORY, THEN
    RUN THE PROGRAM OVER THE SAME REGION. THE AREA MIGHT BE SLIGHTLY DIFFERENT, BUT THIS CAN BE CUT LATER.
    

"""

# General Module Imports
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.map as mp
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunkit_image import coalignment 
import matplotlib.patches as patches
import matplotlib
matplotlib.use('Qt5Agg')
from pathlib import Path

# Local Module Imports
from QOL_Programs.Build_Data_SavePath import buildsavepath_solardata
from QOL_Programs.GUI_Enter_Information import gather_info


info = gather_info()
instrument = info.instrument
wavelength = info.wavelength
starttime = info.starttime
endtime = info.endtime
# # JSOC Email address is needed . . . 
# # We will use our submap to create a cutout request from JSOC (WE WILL NEED A 
# # EMAIL REGISTERED WITH JSOC). IF YOU DON'T HAVE ONE, PLEASE MAKE ONE AT THE 
# # FOLLOWING LINK : http://jsoc.stanford.edu/ajax/register_email.html
jsoc_email = info.jsoc_email

timeaia = a.Time(starttime, endtime)

if instrument == 'HMI':
    raise ValueError('The instument should be AIA, not HMI.')

def create_time_series_AIA(timeaia, wavelength, jsoc_email):
    
    # First we query a single full-disk AIA image
    query = Fido.search(timeaia,
                        a.Instrument('AIA'),
                        a.Wavelength(wavelength*u.angstrom),
                        a.Sample(20000 * u.s)
                        )
    files = Fido.fetch(query)
    # We trick the full disk image into being a sequence to allow for
    # matplotlib plot functions to be used (normally they can only be used on 
    # animated plots).
    amap = mp.Map(files, sequence = True)
    
    # The solar data where the user chooses the 
    # rea/object of interest to build 
    # the time series over.
    print('\nSelect an area of interest to build the time-series.\n')
    fig = plt.figure()
    global ax
    ax = fig.add_subplot(projection = amap[0])
    plt.title('Choose the area of interest.')
    ani = amap.plot()
    global num
    num = 0
    plot_function = interactive_plot_regionselect(fig, ax, sunpymap = amap)
    plt.waitforbuttonpress()
    while True:
        plt.show()
        if plt.waitforbuttonpress():
            # In order to limit large regions of interest, the area of the chosen region is considered.
            if abs(xf-xb)*abs(yf-yb) >= 1000000 :
                areyousure = input('\nThe chosen area is LARGE.\nBuilding a time series this large is going to save a LOT of data.\nAre you sure you want to build a time series this large?\n(Y,N)')
                if areyousure == 'Y' or areyousure == 'y' : 
                    plt.close()
                    break
                if areyousure != 'Y' or areyousure != 'y':
                    plt.waitforbuttonpress()
            else:
                plt.close()
                break

    ref_pix = (amap[0].reference_pixel)*u.pix

    xb_as, xf_as, yb_as, yf_as = (xb-ref_pix[0]/u.pix)*0.6*u.arcsec, (xf-ref_pix[0]/u.pix)*0.6*u.arcsec, (yb-ref_pix[1]/u.pix)*0.6*u.arcsec, (yf-ref_pix[1]/u.pix)*0.6*u.arcsec
    aia_bottom_left = SkyCoord(xb_as, yb_as, frame = amap[0].coordinate_frame)
    aia_top_right = SkyCoord(xf_as, yf_as, frame = amap[0].coordinate_frame)

    print('\nThe area of interest has been chosen.\nThe coordinates (x1,y1) (x2,y2) in Arc-Seconds are :' + ' (' + str(xb_as/u.arcsec) + ', ' + str(yb_as/u.arcsec) + ') (' + str(xf_as/u.arcsec) + ', ' + str(yf_as/u.arcsec) + ').')

    # Next create a submap from this single image. We want to crop the FOV of the 
    # area of interest by using a SunPY submap. This WILL be used to query the
    # JSOC cutouts.
    #aia_bottom_left = SkyCoord(x1, y1, frame = amap.coordinate_frame)
    #aia_top_right = SkyCoord(x2, y2, frame = amap.coordinate_frame)
    
    
    cutouts = amap[0].submap(bottom_left = aia_bottom_left, 
                             top_right = aia_top_right)
    
    # Now we wish to obtain a time series of the cutouts and view the evolution of 
    # the region in time. We also want to AVOID downloading the full disk image. . . 
    # as this is a TON of data and most computers can't handle this.
    cutout = a.jsoc.Cutout(bottom_left = cutouts.bottom_left_coord, 
                           top_right = cutouts.top_right_coord,
                           tracking = True,
                           register = True
                           )
    
    timelength = float(input('How many hours of data do you want?\n'))
    cadence = int(input('Cadence of the data set (in seconds)?\n'))
    query = Fido.search(a.Time(cutouts.date - (timelength/2)*u.h, cutouts.date + (timelength/2)*u.h),
                        a.Wavelength(cutouts.wavelength),
                        a.Sample(cadence*u.s),
                        a.jsoc.Series.aia_lev1_euv_12s,
                        a.jsoc.Notify(jsoc_email),
                        cutout
                        )

    # Submit the export request and download the data.
    files = Fido.fetch(query)
    files.sort()
    
    # Now we create a MapSequence from the downloaded files.
    sequence = mp.Map(files, sequence = True)
    
    # Finally, construct an animation in time from out stack of cutouts and 
    # interactively flip through each image in the sequence. Uncomment this 
    # to view the sequence of images. They are currently ONLY aligned using the 
    # JSOC procedure. We need to align them after (Uncomment and compare to see
    # why I believe the next function aligns the images better). 
    #plt.figure()
    #ani = sequence.plot()

    # In order to make sure that the data is coaligned with respect to the solar
    # rotation  . . . even though the 'tracking' option is considered with the JSOC
    # download of the data, this should be used as a secondary measure to ensure that
    # the data is perfectly aligned. In my comparison, the data appears to be
    # better aligned using this function than without it.
    print('Coaligning the data . . . may take a second.')
    derotated_sequence = coalignment.mapsequence_coalign_by_rotation(sequence, layer_index=int(len(sequence)/2-1))


    # Viewing the newly aligned data. This is purely for visual purpose.
    #plt.figure()
    #ani = derotated_sequence.plot()
    #plt.waitforbuttonpress()
    #while True:
    #    plt.show()
    #    if plt.waitforbuttonpress():
    #        plt.close()
    #        break
    
    # Saving the maps into a local filepath (they also save into your sunpy 
    # data user directory). But I find it more convienent to save it locally.
    pathrot = buildsavepath_solardata(sunpymap_sequence = derotated_sequence) + '\\Coalignbyrotation\\'
    Path(pathrot).mkdir(parents=True, exist_ok=True)
    
    derotated_sequence.save(pathrot + 'Submap_{index}.fits')
    
    # Also, I want to save the un-coaligned sequence so that I can see the difference between
    # the two in future references. This is not needed, but if you wish to do the
    # same feel free to uncomment the following line.
    pathorig = buildsavepath_solardata(sunpymap_sequence = derotated_sequence) + '\\NoCoalignbyrotation\\'
    Path(pathorig).mkdir(parents=True, exist_ok=True)
    
    sequence.save(pathorig + 'Submap_{index}.fits')

    
    
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




if __name__ == '__main__' :
    
    create_time_series_AIA(timeaia, wavelength, jsoc_email)