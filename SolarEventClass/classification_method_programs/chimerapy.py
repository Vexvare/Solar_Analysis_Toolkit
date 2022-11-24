# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 11:29:20 2021

@author: Landon Wells

Completed on 1/30/2022 and cleaned up on 2/3/2022 and 3/16/2022 and 8/17/2022

Purpose: Replicate an (somewhat) exact copy of the IDL CHIMERA program created by the authors below in python.

"""

#+ Original IDL comments 
# Project     : Coronal Hole Identification via Multi-thermal Emission Reconition Algorithm (CHIMERA)
#
# Name        : chimera
#
# Purpose     : Generate a coronal hole segmented tricolour image and corresponding property .txt file
#
# Syntax      : chimera
#
# Inputs      : 171A .fits file
#		193A .fits file
#		211A .fits file
#		hmi .fits file
#
# Outputs     : chimera.png
#	      : chimera.txt
#
# Keywords    : temp= input directory
#		outpath= output directory
#		track= directory containing tracking information
#
# Example    : IDL> chimera,temp='location/string/', outpath='location/string/', track='location/string/'
#                  --------- we in python territory now . . . you use MY CODE or none at all
#
# History     : Written 01-jun-2016, Tadhg Garton, TCD
#
# Contact     : gartont@tcd.ie
#		info@solarmonitor.org
#
#-

import numpy as np
import sunpy.map as mp
from astropy import units as u
import warnings
warnings.filterwarnings('ignore')
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import copy
from matplotlib.path import Path
from astropy.wcs import WCS
from reproject import reproject_interp
# Note! Shapely only works specifically in the cartesian plane. 
# To convert from one plane to another this module will not work.
from shapely import geometry

def Chimera(aia171, aia193, aia211, *args):
    print('\nUsing CHIMERA to find the locations of the Coronal Hole boundaries . . .\n')

    if(aia171 == '' or aia193 == '' or aia211 == ''):
    	raise ValueError('Not all files (AIA 171, 193, 211) are present.')
    ############### Reading in the data #################
    
    #hin = hmi.meta
    #hd = hmi._data
    
    # The AIA data can NOT be combined into a single object due to the fact 
    # that the 'maybase.py' does not currently have a '__getitem__' function...
    # Instead, I choose to keep the three AIA map objects seperate
    
    # Deep copy's of the data are needed to be made to protect the integrity of 
    # the original maps that were downloaded using fido.Search()
    aia171_copy = copy.deepcopy(aia171)
    aia193_copy = copy.deepcopy(aia193)
    aia211_copy = copy.deepcopy(aia211)
    
    ind1 = aia171_copy.meta
    ind2 = aia193_copy.meta
    ind3 = aia211_copy.meta
    
    aia171data = aia171_copy._data
    aia193data = aia193_copy._data
    aia211data = aia211_copy._data
    
    
    ############ Rotates magnetograms if necessary ############
    # The magnetogram should already be rotated. If you have yet to rotate it 
    # use this Sunpy command on your hmi data : hmi_rot = hmi.rotate(order = 3)
    # Alongside that, if your data is not already level 1.5 you can use this 
    # aiapy library command to make it that 
    # level : aia_xxx = register(#your_aia_data)
    
    
    ############ Resize and smooth image(s) #############
    
    xrange_resample = (1024)*(u.pix)
    yrange_resample = (1024)*(u.pix)
    shape = u.Quantity([yrange_resample, xrange_resample])
    data1_r1 = aia171_copy.resample(shape, method='linear')
    data2_r1 = aia193_copy.resample(shape, method='linear')
    data3_r1 = aia211_copy.resample(shape, method='linear')
    
    
    ############# Alternative Coordinate system ###########
    
    #wcs = aia193_copy.reference_coordinate
    ind1['naxis1'] = 4096
    ind1['naxis2'] = 4096
    ind2['naxis1'] = 4096
    ind2['naxis2'] = 4096
    ind3['naxis1'] = 4096
    ind3['naxis2'] = 4096
    xrange_resample = (4096)*(u.pix)
    yrange_resample = (4096)*(u.pix)
    shape = u.Quantity([yrange_resample, xrange_resample])
    data1_r2 = data1_r1.resample(shape, method='linear')
    data2_r2 = data2_r1.resample(shape, method='linear')
    data3_r2 = data3_r1.resample(shape, method='linear')
    #transf = wcs.transform_to('heliographic_stonyhurst')    
    data1 = data1_r2._data
    data2 = data2_r2._data
    data3 = data3_r2._data
    s = data2.shape
    
    frame_out171 = SkyCoord(0,0,unit = u.deg,
                         frame = "heliographic_stonyhurst",
                         obstime = aia171_copy.date,
                         rsun = aia171_copy.coordinate_frame.rsun
                         )
    header_car171 = mp.make_fitswcs_header(s, frame_out171,                    
                                    scale = (180 / s[0], 
                                             180 / s[1]) * u.deg/u.pix,
                                    projection_code = "CAR"
                                    )
    out_wcs171 = WCS(header_car171)
    aia171_copy2 = copy.deepcopy(aia171_copy)
    
    # reproject_interp is the fastest; yet the least accurate method of reprojecting data between coordinate frames.
    # for more accurate, yet slower reprojection it is recommended to use reproject_adaptive
    coord171, footprint171 = reproject_interp(aia171_copy2,out_wcs171, shape_out=s)
    
    frame_out193 = SkyCoord(0,0,unit = u.deg,
                         frame = "heliographic_stonyhurst",
                         obstime = aia193_copy.date,
                         rsun = aia193_copy.coordinate_frame.rsun
                         )
    header_car193 = mp.make_fitswcs_header(s, frame_out193,                    
                                    scale = (180 / s[0], 
                                             180 / s[1]) * u.deg/u.pix,
                                    projection_code = "CAR"
                                    )
    out_wcs193 = WCS(header_car193)
    aia193_copy2 = copy.deepcopy(aia193_copy)
    
    coord193, footprint193 = reproject_interp(aia193_copy2,out_wcs193, shape_out=s)
    
    frame_out211 = SkyCoord(0,0,unit = u.deg,
                         frame = "heliographic_stonyhurst",
                         obstime = aia211_copy.date,
                         rsun = aia211_copy.coordinate_frame.rsun
                         )
    header_car211 = mp.make_fitswcs_header(s, frame_out211,                    
                                    scale = (180 / s[0], 
                                             180 / s[1]) * u.deg/u.pix,
                                    projection_code = "CAR"
                                    )
    out_wcs211 = WCS(header_car211)
    aia211_copy2 = copy.deepcopy(aia211_copy)
    
    coord211, footprint211 = reproject_interp(aia211_copy2,out_wcs211, shape_out=s)
    

    # # Just functions to make sure the reprojection is correct.. 
    # feel free to verify for yourself that the reprojection is correctly done
    # outmap.plot_settings = aia193_copy2.plot_settings
    # fig = plt.figure()
    # ax = plt.subplot(projection = outmap)
    # outmap.plot(ax)
    # ax.set_xlim(0,s[1])
    # ax.set_ylim(0,s[0])

    
    # This method is to get the 'coord' array of the coordinate pixel values in
    # arcseconds . . . If you are importing level 1.5 data the number's wont match
    # EXACTLY with those found in the IDL program ; but they should be fairly close
    arrcol = np.zeros((s))
    arrrow = np.zeros((s))
    x_scl = ind2['cdelt1']
    y_scl = ind2['cdelt2']
    x=0
    
    x_ref = ind2['CRVAL1']
    y_ref = ind2['CRVAL2']
    numx = -(len(data2[0,:])*.6)/2 + x_ref
    numy = -(len(data2[:,0])*.6)/2 + y_ref
    while(x < len(data2)):
        
        arrcol[x,:] = numx
        arrrow[:,x] = numy
        numx = numx + x_scl
        numy = numy + y_scl
        x = x+1

    xcol = arrcol
    ycol = arrrow
    coord = np.zeros((2,s[0],s[1]))
    coord[0,:,:] = xcol
    coord[1,:,:] = ycol
    
    
    ############# Setting up arrays to be used ############
    ident = 1 # Used to track the numerical-identity of the coronal holes
    iarr = np.zeros((s[0],s[1]))
    offarr = np.zeros((s[0],s[1]))
    onarr = np.zeros((s[0],s[1]))
    bound_by_index = np.zeros((s[0],s[1]))
    #ch_of_interest = np.zeros((s[0],s[1]))
    bound_arr_FD = np.zeros((s[0],s[1]))
    bound_arr_FD_w_holes = np.zeros((s[0],s[1]))
    areas_list = []
    mas = np.zeros((s[0],s[1]))
    mak = np.zeros((s[0],s[1]))
    msk = np.zeros((s[0],s[1]))
    #tem = np.zeros((4000,4000))
    #tmp = np.zeros((4000,4000))
    #tep = np.zeros((4000,4000))
    deff = np.zeros((s[0],s[1]))
    circ = np.zeros((s[0],s[1]), dtype=int)
    #n = np.zeros(1, dtype=int)
    #x = np.zeros(1)
    #y = np.zeros(1)
    #ch = np.zeros(1, dtype=int)

    
    ############### Creation of a 2d gaussian for magnetic cut offs ##################
    r = (s[0]/2.0)-450
    onesx = np.ones((s[1]))
    onesy = np.ones((s[0]))
    xranges = np.arange(0,s[1])
    yranges = np.arange(0,s[0])
    xgrid, ydummy = np.meshgrid(xranges,onesx)
    xdummy, ygrid = np.meshgrid(onesy,yranges)
    center = ([int(s[0]/2),int(s[1]/2)])
    wy = np.where(((xgrid-center[1])**2 + (ygrid-center[0])**2 > r**2 ))[0]
    wx = np.where(((xgrid-center[1])**2 + (ygrid-center[0])**2 > r**2 ))[1]
    circ[wy,wx] = 1
    # Astropy comes with a function titled Gaussian2D that we will use to make
    # the point spread function gaussians with a specified FWHM
    # we want the FWHM to be of size [2000,2000] for x and y respectivly and
    # the entire gaussian 2d needs to be the size of the x and y pixel values
    # s[1] and s[0] respectivly
    #x_stddev = (2000.0/(2*np.sqrt(2*np.log(2))))
    #y_stddev = (2000.0/(2*np.sqrt(2*np.log(2))))
    #garr0 = models.Gaussian2D(y_stddev = y_stddev, x_stddev = x_stddev)
    #garr0.amplitude = 1
    #garr01 = garr0(ygrid,xgrid)
    #garr01[wy,wx]=1
    #garr02 = garr01
    #garr02[wy,wx]=1
    
    
    ############### Creation of array for CH properties ##################
    #... will be implemented at a later date. for now i'm going to focus on the rest of the program 
    
    
    ############### Sort data by wavelength ###############
    # The data is already in seperated array variables; so the data doesn't need to be sorted by the wavelength.
    
    
    ############## Normalize the data with respect to exposure time #########
    exptime171 = ind1['exptime']
    new_meta1 = copy.deepcopy(ind1)
    #new_meta1['exptime'] = 1.0
    m_normalized_171 = data1/exptime171
    
    exptime193 = ind2['exptime']
    new_meta2 = copy.deepcopy(ind2)
    #new_meta2['exptime'] = 1.0
    m_normalized_193 = data2/exptime193
    
    exptime211 = ind3['exptime']
    new_meta3 = copy.deepcopy(ind3)
    #new_meta3['exptime'] = 1.0
    m_normalized_211 = data3/exptime211
    
    
    ######### removes negative data values ##########
    data171 = m_normalized_171
    data193 = m_normalized_193
    data211 = m_normalized_211
    
    data171 = np.nan_to_num(data171)
    data193 = np.nan_to_num(data193)
    data211 = np.nan_to_num(data211)
    
    data171[np.where(data171<0)] = 0 
    data193[np.where(data193<0)] = 0
    data211[np.where(data211<0)] = 0 
    
    
    ################## Readies maps, specifies solar radius and calculates conversion value of pixel to arcsec #############
    # Not needed. They should already be SunPY maps.
    #map171 = mp.Map(data171, new_meta1)
    #map193 = mp.Map(data193, new_meta2)
    #map211 = mp.Map(data211, new_meta3)

    rs = new_meta2['rsun_obs']

    if( new_meta2['cdelt1'] > 1 ):
        new_meta1['cdelt1'] = new_meta1['cdelt1']/4
        new_meta2['cdelt1'] = new_meta2['cdelt1']/4
        new_meta3['cdelt1'] = new_meta3['cdelt1']/4
        
        new_meta1['cdelt2'] = new_meta1['cdelt2']/4
        new_meta2['cdelt2'] = new_meta2['cdelt2']/4
        new_meta3['cdelt2'] = new_meta3['cdelt2']/4
        
        new_meta1['crpix1'] = new_meta1['crpix1']*4
        new_meta2['crpix1'] = new_meta2['crpix1']*4
        new_meta3['crpix1'] = new_meta3['crpix1']*4
        
        new_meta1['crpix2'] = new_meta1['crpix2']*4
        new_meta2['crpix2'] = new_meta2['crpix2']*4
        new_meta3['crpix2'] = new_meta3['crpix2']*4
 
    dattoarc = new_meta2['cdelt1']


    ####################### Seperates each image to an individual array ##############
    # No need to do this because the images should be seperated
    #dat0=data171
    #dat1=data193
    #dat2=data211
        
    
    ###################### Get pixels with useful intensities and on disk #################
    r = new_meta1['r_sun'] 
    x2y2 = (xgrid-center[1])**2 + (ygrid-center[0])**2
    wy = np.where(((data171 < 4000) & (data193 < 4000) & (data211 < 4000) & (x2y2 < r**2)))[0]
    wx = np.where(((data171 < 4000) & (data193 < 4000) & (data211 < 4000) & (x2y2 < r**2)))[1]
    
    
    ########## Create intensity ratio arrays
    # Commented out for now as this is the part of the program tha takes the 
    # longest to run and doesn't seem to be required for the program . . .
    # i=0
    #while( (i < len(wx)-1) ):
    #    tem[int(round(data171[wy[i],wx[i]])),
    #        int(round(data193[wy[i],wx[i]]))] = tem[int(round(data171[wy[i],wx[i]])),
    #                                                int(round(data193[wy[i],wx[i]]))]+1
    #    tmp[int(round(data171[wy[i],wx[i]])),
    #        int(round(data211[wy[i],wx[i]]))] = tmp[int(round(data171[wy[i],wx[i]])),
    #                                                int(round(data211[wy[i],wx[i]]))]+1
    #    tep[int(round(data171[wy[i],wx[i]])),
    #        int(round(data211[wy[i],wx[i]]))] = tep[int(round(data193[wy[i],wx[i]])),
    #                                                int(round(data211[wy[i],wx[i]]))]+1
    #    i=i+1
    
    ########## Make a multi-wavelength image for contours ################
    log171 = np.log10(data171)
    log193 = np.log10(data193)
    log211 = np.log10(data211)
    
    truecolorimage171logscaled = bytescale(log171, cmin = 1.2, cmax = 3.9)
    truecolorimage193logscaled = bytescale(log193, cmin = 1.4, cmax = 3.0)
    truecolorimage211logscaled = bytescale(log211, cmin = 0.8, cmax = 2.7)
    
    t0 = truecolorimage211logscaled
    t1 = truecolorimage193logscaled
    t2 = truecolorimage171logscaled
    
    
    ######### Create a 3 segmented bitmasks bitmasks #############
    condit1 = (np.mean(data171)*0.6375)/(np.mean(data211))
    condit2 = 0.7*(np.mean(data193)+np.mean(data211))
    condit3 = (np.mean(data171)*1.5102)/(np.mean(data193))

    test1 = np.divide(t2.astype(np.float),t0.astype(np.float),
                      out = t2.astype(np.float),
                      where=(t0.astype(np.float)>0)
                      )
    test2 = (t0 + t1).astype(np.float)
    test3 = np.divide(t2.astype(np.float),t1.astype(np.float), 
                      out = t2.astype(np.float), 
                      where=(t1.astype(np.float)>0)
                      )
    
    msk[np.where(test1 >= condit1)] = 1
    mak[np.where(test2 < condit2)] = 1
    mas[np.where(test3 >= condit3)] = 1
    
    
    ######### Plot tricolor image with lon/lat contours ##########
    # I could not get this to plot a proper tricolor image . . . so I will not be using this image 
    #tricolormap0 = mp.Map(t0, new_meta3)
    #tricolormap1 = mp.Map(t1, new_meta2)
    #tricolormap2 = mp.Map(t2, new_meta1)
    #comp_map = mp.Map(tricolormap2, tricolormap1, tricolormap0, composite=True)
    #comp_map.set_alpha(0, .33)
    #comp_map.set_alpha(1, .33)
    #comp_map.set_alpha(2, .33)
    #comp_map.set_zorder(0,0)
    #comp_map.set_zorder(1,0)
    #comp_map.set_zorder(2,0)
    
    
    ########## Removes off detector mis-identifications and seperates on-disk and off-limb CH's ###########
    circ[:] = 1
    rm = (s[0]/2.0) - 100
    r = rs/dattoarc
    xgrid, ydummy = np.meshgrid(xranges,onesx)
    xdummy, ygrid = np.meshgrid(onesy,yranges)
    center = ([int(s[0]/2),int(s[1]/2)])
    xxyy = (xgrid-center[1])**2 + (ygrid-center[0])**2
    wy = np.where( (xxyy >= rm**2) | ((xxyy >= (r-10)**2) & (xxyy <= (r+40)**2)) )[0]
    wx = np.where( (xxyy >= rm**2) | ((xxyy >= (r-10)**2) & (xxyy <= (r+40)**2)) )[1]
    county = len(wy)
    if( county > 0 ):
        circ[wy,wx] = 0
    # This array is used to identify the possible CHs via thresholding intensities
    deff = mas*msk*mak*circ 
    
    
    ######### Contours the identified datapoints ##########
    levels = [0.5,1.5] 
    cs = plt.contourf(deff, levels = levels, alpha = 0)
    # Removing the contour lines
    for c in cs.collections:
        c.set_edgecolor("face") 
    
    
    ######### plotting the lon/lat contours ##########
    #comp_map.draw_grid(axes=ax)
    #aia193.draw_grid(axes=ax)
    
    
    ######### Cycles through contours ############
    # The lengths of the contour arrays identify possible areas containing CHs
    contour = cs.collections[0]
    info = contour.get_paths()
    plt.show()
    plt.close()
    
    main_ch_ = False
    
    for j in range(len(info)):
        
        
        ########### Only takes values of minimum surface length and calculates area ###########
        if( len(info[j]) >= 100 ):
            area_of_individual_polygons = []
            inside_boundaries_xs = []
            offs = info[j].vertices
            
            
            #########################################################
            # Calculating the area of the coronal hole with considerations of the islands
            ##############          This next bit of code is rather complicated; to understand 
            ###########       what it does it is recommended to read the following stackoverflow discussion : 
            ########## https://stackoverflow.com/questions/48634934/contour-area-calculation-using-matplotlib-path
            sign = 1
            codes = info[j].codes
            idx = np.where(codes == Path.MOVETO)[0]
            if(len(idx) != 0):
                vert_segs = np.split(offs,idx)[1:]
                code_segs = np.split(codes,idx)[1:]
                for code, vert in zip(code_segs,vert_segs):
                    
                    
                    # This will show exactly where the CH boundary and islands are located... 
                    # Uncomment the function below to see where they exist.
                    #new_path = Path(vert,code)
                    #patch = PathPatch(
                    #    new_path,
                    #    edgecolor = 'black' if sign == 1 else 'midnightblue',
                    #    facecolor = 'none',
                    #    lw = 1)
                    #ax.add_patch(patch)
                    
                    # Computing the exact coronal hole area (WITHOUT the area of the islands)
                    area_of_individual_polygons.append(sign*geometry.Polygon(vert[:-1]).area)
                    if(sign == 1):
                        xs = offs[:,1]
                        ys = offs[:,0]
                        outside_vert = vert[:-1]
                    if(sign == -1):
                        inside_boundaries_xs.append(vert[:-1])
                    sign = -1
            total_area_no_holes = np.sum(area_of_individual_polygons)
            arcar = total_area_no_holes *dattoarc *dattoarc
            ########################## End of complex code ########################################
            
            
            if(arcar > 10000):
                ################ Find's Centroid ##################
                chptsy = np.arange(len(ys), dtype = int)
                chptsx = np.arange(len(xs), dtype = int)
                cent = [np.mean(ys[chptsy]).astype('int'),np.mean(xs[chptsx]).astype('int')]
                
                
                ################ Removes quiet sun regions encompassed by coronal holes ##############
                yy1 = (max( ys[chptsy] )+1).astype('int')
                xx1 = (xs[min( np.where(ys[chptsy] == max(ys[chptsy])) )]).astype('int')
                iarry1 = (max(ys[chptsy])+1).astype('int')
                iarrx1 = (xs[min(np.where(ys[chptsy] == max(ys[chptsy])))]).astype('int')
                
                if( np.any(deff[yy1,xx1]>0) == True and np.any(iarr[iarry1,iarrx1]>0) == True ):
                    # Creating the exterior of a polygon using the boundary of the coronal hole
                    # This converts the contour points into a shapely object to be used in finding the boundaries . . . 
                    subscripts = geometry.Polygon(outside_vert[:-1], inside_boundaries_xs) 
                    # Max and min values of the boundary (IN PIXEL-COORDS, NOT INTENSITY) defined by the polygon.
                    minx, miny, maxx, maxy = subscripts.bounds 
                    # Calculating all of the lattice points that lie both inside and on the boundary of the given polygon . . . 
                    x_c = np.arange(np.floor(minx),np.ceil(maxx)+1)
                    y_c = np.arange(np.floor(miny),np.ceil(maxy)+1)
                    points = geometry.MultiPoint(np.transpose([np.tile(x_c,len(y_c)), np.repeat(y_c,len(x_c))]))
                    result = (points.intersection(subscripts))
                    list_arr = [np.array((geom.xy[0][0],geom.xy[1][0])) for geom in result]
                    # Adding the boundary points to list_arr just incase the points were not considered in the algorithm above
                    list_boundary_x, list_boundary_y = subscripts.exterior.coords.xy
                    
                    # Appending the x and y pixel coordinates 
                    boundary_arr = np.zeros((len(list_boundary_x),2))
                    ii=0
                    while(ii < len(list_boundary_x)):
                        boundary_arr[ii,0] = int(list_boundary_x[ii])
                        boundary_arr[ii,1] = int(list_boundary_y[ii])
                        ii += 1
                    new_arr = np.concatenate((list_arr,boundary_arr), axis = 0)
                    
                    ###################################################################################################################
                    ################ Just to recap what happened above; we now currently have THREE arrays of interest ################
                    ################ boundary_arr is the array of PURELY THE BOUNDARY POINTS ##########################################
                    ################ list_arr is the array of the INSIDE CORONAL HOLE POINTS ##########################################
                    ################ ch_arr is the array of the ENTIRE CORONAL HOLE INCLUDING THE BOUNDARY ############################
                    ################ ALL ARRAYS ARE IN TERMS OF PIXEL COORDS ##########################################################
                    ###################################################################################################################
                    
                    # Adding the coronal holes points into iarr
                    iy = 0
                    while(iy < len(list_arr)):
                        iarr[int(list_arr[iy][1]),int(list_arr[iy][0])] = 0
                        iy += 1
                        
                        
                    ############## Create a simple center point ##############   
                    
                else:
                    
                    
                    arccent0 = coord[0,cent[0],0] 
                    arccent1 = coord[1,0,cent[1]]
                    
                    
                    ############## Classifies off limb CH regions ############
                    
                    test_ys = coord[0,ys[chptsy].astype('int'),0]**2
                    test_xs = coord[1,0,xs[chptsx].astype('int')]**2
                    
                    if( ((arccent0**2 + arccent1**2) > rs**2) | (np.any((test_ys + test_xs) > rs**2) == True) ):
                        # Creating the exterior of a polygon using the boundary of the coronal hole
                        #keep in mind; this converts the contour points into a shapely object . . . 
                        subscripts = geometry.Polygon(outside_vert[:-1], inside_boundaries_xs) 
                        #Max and min values of the boundary defined by the polygon
                        minx, miny, maxx, maxy = subscripts.bounds 
                        # Calculating all of the lattice points that lie both inside and on the boundary of the given polygon . . . 
                        x_c = np.arange(np.floor(minx),np.ceil(maxx)+1)
                        y_c = np.arange(np.floor(miny),np.ceil(maxy)+1)
                        points = geometry.MultiPoint(np.transpose([np.tile(x_c,len(y_c)), np.repeat(y_c,len(x_c))]))
                        result = (points.intersection(subscripts))
                        list_arr = [np.array((geom.xy[0][0],geom.xy[1][0])) for geom in result]
                        # Adding the boundary points to list_arr just incase the points were not considered in the algorithm above
                        list_boundary_x, list_boundary_y = subscripts.exterior.coords.xy
                        boundary_arr = np.zeros((len(list_boundary_x),2))
                        ii=0
                        while(ii < len(list_boundary_x)):
                            boundary_arr[ii,0] = int(list_boundary_x[ii])
                            boundary_arr[ii,1] = int(list_boundary_y[ii])
                            ii += 1
                        new_arr = np.concatenate((list_arr,boundary_arr), axis = 0)
                        
                        # Adding the coronal holes points into offarr
                        iy = 0
                        while(iy < len(new_arr)):
                            offarr[int(new_arr[iy][1]),int(new_arr[iy][0])] = 1
                            iy += 1
                    
                        
                    
                    ############ Classifies on disk coronal holes ############
                    
                    else:
                        
                        print('The area [in arcseconds] of the on-disk coronal hole is : ' , arcar, '[arcseconds]')
                        
                        areas_list.append(arcar)
                        
                        # Creating the exterior of a polygon using the boundary of the coronal hole
                        subscripts = geometry.Polygon(outside_vert[:-1], inside_boundaries_xs) 
                        minx, miny, maxx, maxy = subscripts.bounds

                        # Calculating all of the lattice points that lie both 
                        # inside and on the boundary of the given polygon . . .
                        x_c = np.arange(np.floor(minx),np.ceil(maxx)+1)
                        y_c = np.arange(np.floor(miny),np.ceil(maxy)+1)
                        points = geometry.MultiPoint(np.transpose([np.tile(x_c,len(y_c)), np.repeat(y_c,len(x_c))]))
                        result = (points.intersection(subscripts))

                        # list_arr is all of the interior coordinates; excluding
                        # the interior boundary of islands and exterior boundaries
                        # of the CH
                        list_arr = [np.array((geom.xy[0][0],geom.xy[1][0])) for geom in result] 

                        
                        # Adding the boundary points to list_arr just incase 
                        # the points were not considered in the algorithm above
                        list_boundary_x, list_boundary_y = subscripts.exterior.coords.xy
                        boundary_arr = np.zeros((int(len(list_boundary_x)),2))

                        ii=0
                        while(ii < len(list_boundary_x)):
                            boundary_arr[ii,0] = int(list_boundary_x[ii])
                            boundary_arr[ii,1] = int(list_boundary_y[ii])
                            ii += 1
                        # This array contains all of the coronal hole points
                        # including boundary and interior points
                        ch_arr = np.concatenate((list_arr,boundary_arr), axis = 0) 
                        
                        interior_coords = []
                        for interior in subscripts.interiors:
                            interior_coords += interior.coords[:]
                        
                        interior_boundary_arr = np.zeros((int(len(interior_coords)),2))
                        i_int = 0
                        while(i_int < len(interior_coords)):
                            interior_boundary_arr[i_int,0] = interior_coords[i_int][0]
                            interior_boundary_arr[i_int,1] = interior_coords[i_int][1]
                            i_int += 1
                        
                        ################ Just to recap what happened above; we now currently have three arrays of interest ################
                        ################ boundary_arr is the array of purely the BOUNDARY POINTS ##########################################
                        ################ list_arr is the array of the INSIDE CORONAL HOLE POINTS ##########################################
                        ### ch_arr is the array of the INDEX NUMBERS/PIXEL COORDINATE OF THE ENTIRE CORONAL HOLE INCLUDING THE BOUNDARY ###
                        
                        # Now we have several LIST ARRAYS; 
                        # interior_boundary_arr / / / which might also be the same as interior_coords (i.e. island boundaries), 
                        # ch_arr which is the BOUNDARY AND INTERIOR points, 
                        # boundary_arr which is PURELY THE BOUNDARY points,
                        # and list_arr which is all of the INTERIOR with OUT boundary points

                        ix = 0
                        while( ix < len(boundary_arr) ):
                            bound_arr_FD[int(boundary_arr[ix][1]),int(boundary_arr[ix][0])] = 1
                            bound_arr_FD_w_holes[int(boundary_arr[ix][1]),int(boundary_arr[ix][0])] = 1
                            bound_by_index[int(boundary_arr[ix][1]),int(boundary_arr[ix][0])] = ident
                            ix += 1
                        ixx = 0
                        while( ixx < len(interior_boundary_arr) ):
                            bound_arr_FD_w_holes[int(interior_boundary_arr[ixx][1]),int(interior_boundary_arr[ixx][0])] = 1
                            bound_by_index[int(interior_boundary_arr[ixx][1]),int(interior_boundary_arr[ixx][0])] = ident
                            ixx += 1
                            
                        iy = 0 
                        while(iy < len(ch_arr)):
                            onarr[int(ch_arr[iy][1]),int(ch_arr[iy][0])] = 1
                            iarr[int(ch_arr[iy][1]),int(ch_arr[iy][0])] = ident
                            iy += 1

                        ident += 1

                        ############ Finding the largest CH, usually this is the one of interest #################
                        
                        if arcar > 20000:
                            
                            main_ch_ = True
                        
                            ch_boundary_of_interest = np.zeros((s[0],s[1]))
                            
                            iy = 0
                            while(iy < len(boundary_arr)):
                                ch_boundary_of_interest[int(boundary_arr[iy][1]),int(boundary_arr[iy][0])] = 1
                                iy += 1    
                            iy = 0
                            while(iy < len(interior_boundary_arr)) :
                                ch_boundary_of_interest[int(interior_boundary_arr[iy][1]),int(interior_boundary_arr[iy][0])] = 1
                                iy += 1
                                
                        
                        
                        ############ Create an array for magnetic polarity ################
                        # Not implemented due to time constraints
                        
                        ################### Create an accurate center point ##################
                        # Not implemented due to time constraints
                                    
                        ################### Calculate average angle coronal hole is subjected histogram too ###############
                        # Not implemented due to time constraints
                                    
                        ################### Calculate area of CH with minimal projection effects ################
                        # Not implemented due to time constraints
                                    
                        ################### Find CH extent in latitude and longitude #########################
                        # Not implemented due to time constraints
                                    
                        ################## CH centroid in lat/lon #########################
                        # Not implemented due to time constraints
                                    
                        ################## Calculate the mean magnetic field ##################
                        # Not implemented due to time constraints
                                    
                        ################## Finds coordinates of CH boundaries #################
                        # Not implemented due to time constraints
                                    
                        ################## Insertions of CH properties into property array #############
                        # Not implemented due to time constraints
  
    ############ Sets ident back to max value of iarr ##############
    # ident = ident - 1
    
    ############ Looks for a previous segmentation array #############
    # Not implemented due to time constraints
    
    ############ Finds time difference for tracking ################
    # Not implemented due to time constraints
    
    ############ Only track if previous segmentation given ################
    # Not implemented due to time constraints
    
    ############ Calculate centroids of old segmentation ################
    # Not implemented due to time constraints
    
    ############ Rotate old segmented array for comparison ################
    # Not implemented due to time constraints

    if main_ch_ == True:
        return iarr, onarr, bound_arr_FD, bound_arr_FD_w_holes, bound_by_index, ch_boundary_of_interest
    else:
        return iarr, onarr, bound_arr_FD, bound_arr_FD_w_holes, bound_by_index
    


def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    """
    Byte scales an array (image).

    Byte scaling means converting the input image to uint8 dtype and scaling
    the range to ``(low, high)`` (default 0-255).
    If the input image already has dtype uint8, no scaling is done.

    Parameters
    ----------
    data : ndarray
        PIL image data array.
    cmin : scalar, optional
        Bias scaling of small values. Default is ``data.min()``.
    cmax : scalar, optional
        Bias scaling of large values. Default is ``data.max()``.
    high : scalar, optional
        Scale max value to `high`.  Default is 255.
    low : scalar, optional
        Scale min value to `low`.  Default is 0.

    Returns
    -------
    img_array : uint8 ndarray
        The byte-scaled array.

    Examples
    --------
    >>> img = array([[ 91.06794177,   3.39058326,  84.4221549 ],
                     [ 73.88003259,  80.91433048,   4.88878881],
                     [ 51.53875334,  34.45808177,  27.5873488 ]])
    >>> bytescale(img)
    array([[255,   0, 236],
           [205, 225,   4],
           [140,  90,  70]], dtype=uint8)
    >>> bytescale(img, high=200, low=100)
    array([[200, 100, 192],
           [180, 188, 102],
           [155, 135, 128]], dtype=uint8)
    >>> bytescale(img, cmin=0, cmax=255)
    array([[91,  3, 84],
           [74, 81,  5],
           [52, 34, 28]], dtype=uint8)

    """
    if data.dtype == np.uint8:
        return data

    if high < low:
        raise ValueError("`high` should be larger than `low`.")

    if cmin is None:
        cmin = data.min()
    if cmax is None:
        cmax = data.max()

    cscale = cmax - cmin
    if cscale < 0:
        raise ValueError("`cmax` should be larger than `cmin`.")
    elif cscale == 0:
        cscale = 1

    scale = float(high - low) / cscale
    bytedata = (data * 1.0 - cmin) * scale + 0.4999
    bytedata[bytedata > high] = high
    bytedata[bytedata < 0] = 0
    return np.cast[np.uint8](bytedata) + np.cast[np.uint8](low)