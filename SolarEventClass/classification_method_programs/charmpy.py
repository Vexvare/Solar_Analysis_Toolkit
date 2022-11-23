# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:57:56 2022

@author: Landon Wells

Purpose : To make a VERY simplified approach to the CHARM boundary threshold method.


Output : Multiple CHARM 2d numpy map containing all the points where this approach identifies
         the boundaries of the coronal hole.
"""



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

def Charm(sunpymap, intensitylevel, *args):
    print('\nUsing CHARM to find the locations of the Coronal Hole boundaries . . .\n')

    if(sunpymap == ''):
    	raise ValueError('A sunpy map is not given.')
    ############### Reading in the data #################


    sunpymap_copy = copy.deepcopy(sunpymap)
    ind = sunpymap_copy.meta
    sunpymapdata = sunpymap_copy._data

    
    ############ Resize and smooth image(s) #############
    
    xrange_resample = (1024)*(u.pix)
    yrange_resample = (1024)*(u.pix)
    shape = u.Quantity([yrange_resample, xrange_resample])
    data2_r1 = sunpymap_copy.resample(shape, method='linear')
    
    
    ############# Alternative Coordinate system ###########
    
    #wcs = aia193_copy.reference_coordinate
    ind['naxis1'] = 4096
    ind['naxis2'] = 4096
    xrange_resample = (4096)*(u.pix)
    yrange_resample = (4096)*(u.pix)
    shape = u.Quantity([yrange_resample, xrange_resample])
    data2_r2 = data2_r1.resample(shape, method='linear')
    data2 = data2_r2._data
    s = data2.shape
    
    frame_out = SkyCoord(0,0,unit = u.deg,
                         frame = "heliographic_stonyhurst",
                         obstime = sunpymap_copy.date,
                         rsun = sunpymap_copy.coordinate_frame.rsun
                         )
    header_car = mp.make_fitswcs_header(s, frame_out,                    
                                    scale = (180 / s[0], 
                                             180 / s[1]) * u.deg/u.pix,
                                    projection_code = "CAR"
                                    )
    out_wcs = WCS(header_car)
    sunpymap_copy2 = copy.deepcopy(sunpymap_copy)
    
    # reproject_interp is the fastest; yet the least accurate method of reprojecting data between coordinate frames.
    # for more accurate, yet slower reprojection it is recommended to use reproject_adaptive
    coord171, footprint171 = reproject_interp(sunpymap_copy2,out_wcs, shape_out=s)
    
    
    ####################################################################################################################### 
    

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
    x_scl = ind['cdelt1']
    y_scl = ind['cdelt2']
    x=0
    
    x_ref = ind['CRVAL1']
    y_ref = ind['CRVAL2']
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
    bound_arr_FD = np.zeros((s[0],s[1]))
    bound_arr_FD_w_holes = np.zeros((s[0],s[1]))
    deff = np.zeros((s[0],s[1]))
    circ = np.zeros((s[0],s[1]), dtype=int)
    mak = np.zeros((s[0],s[1]))

    
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
    
    
    ############### Creation of array for CH properties ##################
    #... will be implemented at a later date. for now i'm going to focus on the rest of the program 
    
    
    ############### Sort data by wavelength ###############
    # The data is already in seperated array variables; so the data doesn't need to be sorted by the wavelength.
    
    
    ############## Normalize the data with respect to exposure time #########

    exptime = ind['exptime']
    new_meta = copy.deepcopy(ind)
    #new_meta2['exptime'] = 1.0
    m_normalized = data2/exptime

    
    
    ######### removes negative data values ##########
    data = m_normalized
    data = np.nan_to_num(data)
    data[np.where(data<0)] = 0
    
    
    ################## Readies maps, specifies solar radius and calculates conversion value of pixel to arcsec #############
    # Not needed. They should already be SunPY maps.
    #map171 = mp.Map(data171, new_meta1)
    #map193 = mp.Map(data193, new_meta2)
    #map211 = mp.Map(data211, new_meta3)

    rs = new_meta['rsun_obs']

    if( new_meta['cdelt1'] > 1 ):
        new_meta['cdelt1'] = new_meta['cdelt1']/4
        new_meta['cdelt2'] = new_meta['cdelt2']/4
        new_meta['crpix1'] = new_meta['crpix1']*4
        new_meta['crpix2'] = new_meta['crpix2']*4

    dattoarc = new_meta['cdelt1']

    
    ###################### Get pixels with useful intensities and on disk #################
    r = new_meta['r_sun'] 
    x2y2 = (xgrid-center[1])**2 + (ygrid-center[0])**2
    wy = np.where(((data < 4000) & (x2y2 < r**2)))[0]
    wx = np.where(((data < 4000) &  (x2y2 < r**2)))[1]
    
    
    mak[np.where(sunpymap._data <= intensitylevel)] = 1

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
    deff = circ*mak
    
    # Now we take deff and contour only those points.
    ######### Contours the identified datapoints ##########
    fig = plt.figure()
    cs = plt.contour(deff, levels = [intensitylevel])
    
    
    
    
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
            
            
            if(arcar > 1000):
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
                        
                        if arcar > 50000 :
                            
                            ch_boundary_of_interest = np.zeros((s[0],s[1]))
                            
                            iy = 0
                            while(iy < len(boundary_arr)):
                                ch_boundary_of_interest[int(boundary_arr[iy][1]),int(boundary_arr[iy][0])] = 1
                                iy += 1    
                            iy = 0
                            while(iy < len(interior_boundary_arr)) :
                                ch_boundary_of_interest[int(interior_boundary_arr[iy][1]),int(interior_boundary_arr[iy][0])] = 1
                                iy += 1
                            
                            main_ch_ = True

    if main_ch_ == True:
        return iarr, onarr, bound_arr_FD, bound_arr_FD_w_holes, bound_by_index, ch_boundary_of_interest
    else:
        return iarr, onarr, bound_arr_FD, bound_arr_FD_w_holes, bound_by_index