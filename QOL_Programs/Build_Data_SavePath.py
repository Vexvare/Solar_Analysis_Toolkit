# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 18:48:31 2022

@author : Landon Wells

Purpose : To build the paths which which data will be saved. Currently there are 
          three types of data that can be saved:
              - Solar_Data - RAW solar data images.
              - Classification_Data - Data obtained using classification methods.
              - Analysis_Data - Data obtained using analysis of classification method data.
"""

import os
from pathlib import Path
import numpy as np
from astropy import units as u
from os.path import expanduser
import sunpy.map as mp
import platform

# =============================================================================
# Solar Data
# =============================================================================
def buildsavepath_solardata(sunpymap_sequence):
    '''
    Function to build the save path of solar data. All solar data obtained
    with this library should use this function to build the path info that will
    be used to store this information.
    '''
    if platform.system() == 'Windows':
        path_seperator = '\\'
    if platform.system() == 'Linux':
        path_seperator = '/'
        
    # Gathering information to store the data
    # Time of data
    starttime = str(sunpymap_sequence[0].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[0].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[0].meta.get('date-obs'))[17:19]
    endtime = str(sunpymap_sequence[-1].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[-1].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[-1].meta.get('date-obs'))[17:19]
    
    # Wavelength of data
    wavelength = str(sunpymap_sequence[0].meta.get('wavelnth'))
    
    # Coordinate information of the entire map
    bot_left_x = str(sunpymap_sequence[0].bottom_left_coord.Tx)
    bot_left_y = str(sunpymap_sequence[0].bottom_left_coord.Ty)
    top_right_x = str(sunpymap_sequence[0].top_right_coord.Tx)
    top_right_y = str(sunpymap_sequence[0].top_right_coord.Ty)
    lowerleft = bot_left_x + '_'  + bot_left_y
    upperright = top_right_x + '_'  + top_right_y
    
    # Building the path that the data will be saved into.
    path = expanduser('~') + path_seperator + 'SolarToolkitData' + path_seperator +  'Solar_Data' + path_seperator + wavelength + path_seperator + starttime + '_' + endtime + path_seperator + 'Coords_AS_' + lowerleft + '_' + upperright + path_seperator 
    
    Path(path).mkdir(parents=True, exist_ok=True)
    print('The data will be saved in the following path :', path)
    
    return path

# =============================================================================
# Classification Data
# =============================================================================
def buildsavepath_classificationdata(sunpymap_sequence, classification_method):
    '''
    Function to build the save path of classification data. All classification
    data obtained with this library should use this function to build the path 
    info that will be used to store this information.
    '''
    if platform.system() == 'Windows':
        path_seperator = '\\'
    if platform.system() == 'Linux':
        path_seperator = '/'
        
    # Gathering information to store the data
    # Time of data
    starttime = str(sunpymap_sequence[0].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[0].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[0].meta.get('date-obs'))[17:19]
    endtime = str(sunpymap_sequence[-1].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[-1].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[-1].meta.get('date-obs'))[17:19]
    
    # Wavelength of data
    wavelength = str(sunpymap_sequence[0].meta.get('wavelnth'))
    
    # Coordinate information of the entire map
    bot_left_x = str(sunpymap_sequence[0].bottom_left_coord.Tx)
    bot_left_y = str(sunpymap_sequence[0].bottom_left_coord.Ty)
    top_right_x = str(sunpymap_sequence[0].top_right_coord.Tx)
    top_right_y = str(sunpymap_sequence[0].top_right_coord.Ty)
    lowerleft = bot_left_x + '_'  + bot_left_y
    upperright = top_right_x + '_'  + top_right_y
    
    # Building the path that the data will be saved into.
    path = expanduser('~') + path_seperator + 'SolarToolkitData' + path_seperator +  'Classification_Data' + path_seperator + str(classification_method) + path_seperator + wavelength + path_seperator + starttime + '_' + endtime + path_seperator + 'Coords_AS_' + lowerleft + '_' + upperright + path_seperator 
    
    Path(path).mkdir(parents=True, exist_ok=True)
    print('The data will be saved in the following path :', path)

    return path



def buildsavepath_submapclassificationdata(sunpymap_sequence, sunpysubmap_sequence, classification_method):
    '''
    Function to build the save path of analysis data. All analysis of data
    obtained with this library should use this function to build the path 
    info that will be used to store this information.
    '''
    if platform.system() == 'Windows':
        path_seperator = '\\'
    if platform.system() == 'Linux':
        path_seperator = '/'
        
    # Gathering information to store the data
    # Time of data
    starttime = str(sunpymap_sequence[0].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[0].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[0].meta.get('date-obs'))[17:19]
    endtime = str(sunpymap_sequence[-1].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[-1].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[-1].meta.get('date-obs'))[17:19]
    
    # Wavelength of data
    wavelength = str(sunpymap_sequence[0].meta.get('wavelnth'))
    
    # Coordinate information of the submap
    bot_left_x_sub = str(np.round_(sunpysubmap_sequence[0].bottom_left_coord.Tx/u.arcsec, decimals = 0))
    bot_left_y_sub = str(np.round_(sunpysubmap_sequence[0].bottom_left_coord.Ty/u.arcsec, decimals = 0))
    top_right_x_sub = str(np.round_(sunpysubmap_sequence[0].top_right_coord.Tx/u.arcsec, decimals = 0))
    top_right_y_sub = str(np.round_(sunpysubmap_sequence[0].top_right_coord.Ty/u.arcsec, decimals = 0))
    lowerleft_sub = bot_left_x_sub + '_'  + bot_left_y_sub + '_'
    upperright_sub = top_right_x_sub + '_'  + top_right_y_sub
    submap_coords = 'Submap_Coords_AS_' + lowerleft_sub + upperright_sub 
    
    # Building the path that the data will be saved into.
    path = expanduser('~') + path_seperator + 'SolarToolkitData' + path_seperator +  'Classification_Data' + path_seperator + str(classification_method) + path_seperator + wavelength + path_seperator + starttime + '_' + endtime + path_seperator + submap_coords + path_seperator
    
    Path(path).mkdir(parents=True, exist_ok=True)
    print('The data will be saved in the following path :', path)

    return path


# =============================================================================
# Analysis Data
# =============================================================================
def buildsavepath_analysisdata(sunpymap_sequence, analysis_method):
    '''
    Function to build the save path of analysis data. All analysis of data
    obtained with this library should use this function to build the path 
    info that will be used to store this information.
    '''
    if platform.system() == 'Windows':
        path_seperator = '\\'
    if platform.system() == 'Linux':
        path_seperator = '/'
    # Gathering information to store the data
    # Time of data
    starttime = str(sunpymap_sequence[0].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[0].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[0].meta.get('date-obs'))[17:19]
    endtime = str(sunpymap_sequence[-1].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[-1].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[-1].meta.get('date-obs'))[17:19]
    
    # Wavelength of data
    wavelength = str(sunpymap_sequence[0].meta.get('wavelnth'))
    
    # Coordinate information of the entire map
    bot_left_x = str(np.round_(sunpymap_sequence[0].bottom_left_coord.Tx/u.arcsec, decimals = 0))
    bot_left_y = str(np.round_(sunpymap_sequence[0].bottom_left_coord.Ty/u.arcsec, decimals = 0))
    top_right_x = str(np.round_(sunpymap_sequence[0].top_right_coord.Tx/u.arcsec, decimals = 0))
    top_right_y = str(np.round_(sunpymap_sequence[0].top_right_coord.Ty/u.arcsec, decimals = 0))
    lowerleft = bot_left_x + '_'  + bot_left_y
    upperright = top_right_x + '_'  + top_right_y
    
    # Building the path that the data will be saved into.
    path = expanduser('~') + path_seperator + 'SolarToolkitData' + path_seperator +  'Analysis_Data' + path_seperator + str(analysis_method) + path_seperator + wavelength + path_seperator + starttime + '_' + endtime + path_seperator + 'Coords_AS_' + lowerleft + '_' + upperright + path_seperator 
    
    Path(path).mkdir(parents=True, exist_ok=True)
    print('The data will be saved in the following path :', path)

    return path


def buildsavepath_submapanalysisdata(sunpymap_sequence, sunpysubmap_sequence, analysis_method):
    '''
    Function to build the save path of analysis data. All analysis of data
    obtained with this library should use this function to build the path 
    info that will be used to store this information.
    '''
    if platform.system() == 'Windows':
        path_seperator = '\\'
    if platform.system() == 'Linux':
        path_seperator = '/'
        
    # Gathering information to store the data
    # Time of data
    starttime = str(sunpymap_sequence[0].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[0].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[0].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[0].meta.get('date-obs'))[17:19]
    endtime = str(sunpymap_sequence[-1].meta.get('date-obs'))[:4] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[5:7] + '_' + str(sunpymap_sequence[-1].meta.get('date-obs'))[8:13] + str(sunpymap_sequence[-1].meta.get('date-obs'))[14:16] + str(sunpymap_sequence[-1].meta.get('date-obs'))[17:19]
    
    # Wavelength of data
    wavelength = str(sunpymap_sequence[0].meta.get('wavelnth'))
    
    # Coordinate information of the entire map
    bot_left_x = str(np.round_(sunpymap_sequence[0].bottom_left_coord.Tx/u.arcsec, decimals = 0))
    bot_left_y = str(np.round_(sunpymap_sequence[0].bottom_left_coord.Ty/u.arcsec, decimals = 0))
    top_right_x = str(np.round_(sunpymap_sequence[0].top_right_coord.Tx/u.arcsec, decimals = 0))
    top_right_y = str(np.round_(sunpymap_sequence[0].top_right_coord.Ty/u.arcsec, decimals = 0))
    lowerleft = bot_left_x + '_'  + bot_left_y
    upperright = top_right_x + '_'  + top_right_y
    
    # Coordinate information of the submap 
    bot_left_x_sub = str(np.round_(sunpysubmap_sequence[0].bottom_left_coord.Tx/u.arcsec, decimals = 0))
    bot_left_y_sub = str(np.round_(sunpysubmap_sequence[0].bottom_left_coord.Ty/u.arcsec, decimals = 0))
    top_right_x_sub = str(np.round_(sunpysubmap_sequence[0].top_right_coord.Tx/u.arcsec, decimals = 0))
    top_right_y_sub = str(np.round_(sunpysubmap_sequence[0].top_right_coord.Ty/u.arcsec, decimals = 0))
    lowerleft_sub = bot_left_x_sub + '_'  + bot_left_y_sub
    upperright_sub = top_right_x_sub + '_'  + top_right_y_sub
    submap_coords = 'Submap_Coords_AS_' + lowerleft_sub + upperright_sub 
    
    # Building the path that the data will be saved into.
    path = expanduser('~') + path_seperator + 'SolarToolkitData' + path_seperator +  'Analysis_Data' + path_seperator + str(analysis_method) + path_seperator + wavelength + path_seperator + starttime + '_' + endtime + path_seperator + 'Coords_AS_' + lowerleft + '_' + upperright + path_seperator + submap_coords + path_seperator
    
    Path(path).mkdir(parents=True, exist_ok=True)
    print('The data will be saved in the following path :', path)

    return path