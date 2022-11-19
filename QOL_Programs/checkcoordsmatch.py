# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 12:58:42 2022

@author: Landon Wells

Purpose : To determine whether or not a given set of coordinates `closely` matches
            another set of coordinates.
"""

# General Module Imports
import os
from os.path import expanduser
import re as re
from astropy import units as u

def check_coords_match(coords, inputbotleft, inputtopright, tolerance = 0):
    '''
    Purpose : To determine whether or not a given set of coordinates matches input
            coordinates.

    Parameters
    ----------
    coords : string list
        String list containing the coordinates of interest in the order 
        [lower_left_x, lower_left_y, upper_right_x, upper_right_y]. These
        are the same as those generated when saving data using Build_Data_SavePath.py
    inputbotleft : astropy Skycoordinate object
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.bottom_left_coord` on the map object.
    inputtopright : astropy Skycoordinate object
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.top_right_coord` on the map object.

    Returns
    -------
    True if the coords is a match, and false otherwise. 

    '''
    
    bot_left_x = float(coords[0])
    bot_left_y = float(coords[1])
    top_right_x = float(coords[2])
    top_right_y = float(coords[3])

    inputbotleft_x = inputbotleft.Tx/u.arcsec
    inputbotleft_y = inputbotleft.Ty/u.arcsec
    inputtopright_x = inputtopright.Tx/u.arcsec
    inputtopright_y = inputtopright.Ty/u.arcsec
    
    check = 0
    
    if inputbotleft_x - tolerance <= bot_left_x and bot_left_x <= inputbotleft_x + tolerance:
        check += 1
    if inputbotleft_y - tolerance <= bot_left_y and bot_left_y <= inputbotleft_y + tolerance:
        check += 1
    if inputtopright_x - tolerance <= top_right_x and top_right_x <= inputtopright_x + tolerance:
        check += 1
    if inputtopright_y - tolerance <= top_right_y and top_right_y <= inputtopright_y + tolerance:
        check += 1
    
    if check >= 3:
        return True
    else:
        return False
    
def grab_coords_of_data_and_compare_input(inputbotleft, inputtopright, inputstring, tolerance = 0):
    '''
    Purpose : To go through all of the files in SolarToolkitData/Classification_Data
              and see if there is a match in the coords of the data and classified
              event.
    
    Parameters
    ----------
    inputbotleft : astropy Skycoordinate object
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.bottom_left_coord` on the map object.
    inputtopright : astropy Skycoordinate object
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.top_right_coord` on the map object.
    inputstring : list
        List of files 

    Returns
    -------
    matchpathlist : List of strings
        List containing all of the coord matches to the inputtime in the directory
        \\SolarToolkitData\\Classification_Data\\ . These are events that match
        whatever data of interest and can be used to 

    '''
    
    # Walk though all of the files in home/user/SolarToolkitData/Classification_Data/ . . . 
    path = inputstring
    matchpathlist = []
    # Build the bot_left and top_right if there is one.
    coords = re.findall('(-?\d+\.\d+?)', path)
    if coords != [] :
        match = check_coords_match(coords = coords, inputbotleft = inputbotleft, inputtopright = inputtopright, tolerance = 10)
        if match is True :
           matchpathlist.append(path)
        
    return matchpathlist
    
def grab_coords_of_data_and_compare(inputbotleft, inputtopright, tolerance = 0):
    '''
    Purpose : To go through all of the files in SolarToolkitData/Classification_Data
              and see if there is a match in the coords of the data and classified
              event.
    
    Parameters
    ----------
    inputbotleft : astropy Skycoordinate object
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.bottom_left_coord` on the map object.
    inputtopright : astropy Skycoordinate object
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.top_right_coord` on the map object.

    Returns
    -------
    matchpathlist : List of strings
        List containing all of the coord matches to the inputtime in the directory
        \\SolarToolkitData\\Classification_Data\\ . These are events that match
        whatever data of interest and can be used to 

    '''
    
    # Walk though all of the files in home/user/SolarToolkitData/Classification_Data/ . . . 
    rootdir = expanduser('~') + '\\SolarToolkitData\\Classification_Data\\'
    matchpathlist = []
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            path = os.path.join(subdir, file)
            # Build the bot_left and top_right if there is one.
            coords = re.findall('(-?\d+\.\d+?)', path)
            if coords != [] :
                match = check_coords_match(coords = coords, inputbotleft = inputbotleft, inputtopright = inputtopright, tolerance = 10)
                if match is True :
                   matchpathlist.append(path)
        
    return matchpathlist