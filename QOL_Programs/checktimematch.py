# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 17:47:24 2022

@author : Landon Wells

Purpose : To determine whether or not a given time lies inbetween a set start time
          and end time for solar data.
"""

# General Module Imports
from datetime import datetime
from datetime import timedelta
import os
from os.path import expanduser
import re as re


def check_time_match(starttime, endtime, inputstarttime, inputendtime):
    '''
    Purpose : To determine whether or not a given time lies inbetween a set start time
              and end time for solar data. 

    Parameters
    ----------
    starttime : String
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.meta.get('date-obs')`, e.g. 2017-04-16T12:00:09.339063. 
        Should be formattable to a datetime object.
    endtime : String
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.meta.get('date-obs')`, e.g. 2017-04-16T12:00:09.339063.
        Should be formattable to a datetime.fromisoformat object.
    inputstarttime : String
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.meta.get('date-obs')`, e.g. 2017-04-16T12:00:09.339063.
        Should be formattable to a datetime object.
    inputendtime : String
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.meta.get('date-obs')`, e.g. 2017-04-16T12:00:09.339063.
        Should be formattable to a datetime object.

    Returns
    -------
    True if the time is a match, and false otherwise. 

    '''
    # Sometimes telescopes mismatch times, so we add a small correction to 
    # ensure that the events match
    time_change = timedelta(minutes = 10)
    
    starttime = datetime.fromisoformat(starttime)
    endtime = datetime.fromisoformat(endtime)
    inputstarttime = datetime.fromisoformat(inputstarttime) - time_change
    inputendtime = datetime.fromisoformat(inputendtime) + time_change
    
    if inputstarttime <= starttime and endtime <= inputendtime : 
        return True
    else :
        return False
    
    
def grab_times_of_data_and_compare(inputstarttime, inputendtime):
    '''
    Purpose : To go through all of the files in SolarToolkitData/Classification_Data
              and see if there is a match in the time of the data and classified
              event.
    
    Parameters
    ----------
    inputstarttime : String
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.meta.get('date-obs')`, e.g. 2017-04-16T12:00:09.339063. 
        Should be formattable to a datetime object. This time will be used to 
        compare with the times in the directory home/user/SolarToolkitData/Classification_Data 
        to see if there is a match in the solar data of interest`s time to any 
        times in those files.
    inputendtime : String
        Can be obtained from mp.MapSequences or mp.GenericMap objects by using 
        the command `.meta.get('date-obs')`, e.g. 2017-04-16T12:00:09.339063. 
        Should be formattable to a datetime object. This time will be used to 
        compare with the times in the directory home/user/SolarToolkitData/Classification_Data 
        to see if there is a match in the solar data of interest`s time to any 
        times in those files.

    Returns
    -------
    matchpathlist : List of strings
        List containing all of the date matches to the inputtime in the directory
        \\SolarToolkitData\\Classification_Data\\ . These are events that match
        whatever data of interest and can be used to 

    '''
    
    # Walk though all of the files in home/user/SolarToolkitData/Classification_Data/ . . . 
    rootdir = expanduser('~') + '\\SolarToolkitData\\Classification_Data\\'
    matchpathlist = []
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            path = os.path.join(subdir, file)
            # Build the startdate and enddate if there is one.
            date = re.findall('[0-9]{4}' + '_' + '[0-9]{2}' + '_' + '[0-9]{2}' + 'T' + '[0-9]{6}' + '_' + '[0-9]{4}' + '_' + '[0-9]{2}' + '_' + '[0-9]{2}' + 'T' + '[0-9]{6}', path)
            if date != [] :
                starttime = date[0][:13].replace('_', '-') + ':' + date[0][13:15] + ':' + date[0][15:17]
                endtime = date[0][18:31].replace('_', '-') + ':' + date[0][31:33] + ':' + date[0][33:35]
                match = check_time_match(starttime = starttime, endtime = endtime, inputstarttime = inputstarttime, inputendtime = inputendtime)
                if match is True :
                    matchpathlist.append(path)
    return matchpathlist
