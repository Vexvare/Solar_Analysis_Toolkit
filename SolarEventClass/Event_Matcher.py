# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 18:34:15 2022

@author: Landon Wells

Purpose : Function to call that matches solar events (after running and saving
          a classification method to identify events) that match the time of the 
          input data.
"""
# General Module Imports
import sunpy.map as mp

# Local Module Imports
from QOL_Programs.checktimematch import grab_times_of_data_and_compare
from QOL_Programs.checkcoordsmatch import grab_coords_of_data_and_compare_input
from .Solar_Events_Class import SolarEvent

def solar_event_matcher(sunpymap, matchwithalignment = False):
    '''
    This function assumes that sunpy maps identifying the events are saved within
    home/user/Solar_Analysis_Toolkit/Data/Classification_Data/ . . . 
    It first checks to see if the input sunpymap time matches any folder name for
    all of the classification methods using 
    '''
    
    # First we get the list of all classification_methods so we know what to search for
    solevent = SolarEvent(sunpymap)
    list_of_event_classification_methods = solevent.list_of_event_classification_methods
    
    # Now we grab the time of the data from sunpymap and compare with events that
    # are previously constructed.
    sunpymapfirst = sunpymap
    sunpymaplast = sunpymap
    if isinstance(sunpymap, mp.MapSequence):
        sunpymapfirst = sunpymap[0]
        sunpymaplast = sunpymap[-1]
    inputstarttime = sunpymapfirst.meta.get('date-obs')[:-7]
    inputendtime = sunpymaplast.meta.get('date-obs')[:-7]
    
    matchlist = grab_times_of_data_and_compare(inputstarttime = inputstarttime, 
                                               inputendtime = inputendtime)

    sunpymap_list = []
    # Check to ensure there is a match, if not nothing is returned.
    if matchlist != [] :
        # Going through the list of event classification methods
        for event in list_of_event_classification_methods:
            match_event_list = []
            # Now we need to find the indexes of matchlist where the event 
            # substring matches part of that string for that specific index.
            for idx, string in enumerate(matchlist):
                if event in string:
                    match_event_list.append(string)
            # Check match_event_list to ensure there is a item within it
            if match_event_list != []:
                # Create a sunpymap out of all filepaths that were a matched list
                sunpymaps = mp.Map([match_event_list], sequence = True)
                # Append the sunpymaps along with the event that was created from
                # them to allow the user to know what events were found at the time.
                sunpymap_list.append([sunpymaps, event])
    
    
    
    if matchwithalignment == True:
        inputbotleft = sunpymapfirst.bottom_left_coord
        inputtopright = sunpymapfirst.top_right_coord
        # We use the previously made match_event_list to go through the events and 
        # identify ones that match with the alignment. We use this list because
        # it already contains the times of the matches.
        match_list_coords = []
        for i, timematch in enumerate(match_event_list):
            matchlistalign = grab_coords_of_data_and_compare_input(inputbotleft = inputbotleft, 
                                                                   inputtopright = inputtopright, 
                                                                   inputstring = timematch)
            if matchlistalign != []:
                match_list_coords.append(matchlistalign)
                
        sunpymap_list_coords = []
        # Check to ensure there is a match, if not nothing is returned.
        if match_list_coords != [] :
             # Going through the list of event classification methods
             for event in list_of_event_classification_methods:
                 match_event_list_coords = []
                 # Now we need to find the indexes of matchlist where the event 
                 # substring matches part of that string for that specific index.
                 for idx, string in enumerate(match_list_coords):
                     if event in string[0]:
                         match_event_list_coords.append(string)
                 # Check match_event_list to ensure there is a item within it
                 if match_event_list_coords != []:
                     # Create a sunpymap out of all filepaths that were a matched list
                     sunpymaps = mp.Map([match_event_list_coords], sequence = True)
                     # Append the sunpymaps along with the event that was created from
                     # them to allow the user to know what events were found at the time.
                     sunpymap_list_coords.append([sunpymaps, event])
                     
        # From here we just return the list. The list has the labels of events with 
        # which the user can take from as they please.
        return sunpymap_list_coords
        
    # From here we just return the list. The list has the labels of events with 
    # which the user can take from as they please.
    return sunpymap_list