# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:37:21 2022

@author: Landon Wells
"""
# This program creates a quick boundary of the HMI Magnetogram network boundary
# by setting a threshold over the level (strength) of the gaussian magnetic field
# strength of the HMI map using matplotlib contourf routines to define the contours
# and the collections to obtain the paths of those contours. They will be stored
# in a map the same size as the input map, and use levels that represent the 
# radially outward field and radially inward field.

import numpy as np
import matplotlib.pyplot as plt

def gauss_countours(sunpymap, gaussianlevel):
    
    min_mag = np.nanmin(sunpymap._data)
    max_mag = np.nanmax(sunpymap._data)
    mag_size = sunpymap._data.shape
    cs_positive = plt.contourf(sunpymap._data, levels = [min_mag,gaussianlevel])
    cs_negative = plt.contourf(sunpymap._data, levels = [-gaussianlevel,max_mag])
    mag_map_vis = np.zeros((mag_size[0],mag_size[1]))
    
    # Obtaining the contours for the positive (radially outward) contours
    p = cs_positive.collections[0].get_paths()[0]
    v = p.vertices
    mag_x = v[:,0]
    mag_y = v[:,1]
    for i in range(len(v)):
        mag_map_vis[int(mag_y[i]), int(mag_x[i])] = 1
        
    # Obtaining the contours for the negative (radially inwards) contours
    p = cs_negative.collections[0].get_paths()[0]
    v = p.vertices
    mag_x = v[:,0]
    mag_y = v[:,1]
    for i in range(len(v)):
        mag_map_vis[int(mag_y[i]), int(mag_x[i])] = -1
        
    return mag_map_vis