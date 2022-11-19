# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 13:16:30 2022

@author: Landon Wells

# =============================================================================
# # This is the class used to define the different solar events. These events 
# # can be identified using different techniques, as of right now, the supported types 
# # of solar events that are used are : Coronal Holes (CH), Active Regions
# # (AR), Loop structures, and HMI Magnetic Network, and unspecified events. 
# # The idea of the class is to be able to host the different
# # functions used to identify certain aspects of these events, for example, CH 
# # boundary thresholding techniques such as CHIMERA and CHARM. 
# # A solar event (and any event that can be physically observed),
# # needs three things, an EVENT TYPE, a classification, and a method of classification.
# # The method of classification is necessary because without it, events are not
# # able to properly be described. The SolarEvent class will be it's own standalone
# # class, but the args should be either a 'mp.GenericMap' or 'mp.MapSequence
# # object, and specific classification_methods must be called when defining an event.
# =============================================================================

"""

import sunpy.map as mp
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from sunkit_image import coalignment 
from astropy.coordinates import SkyCoord
from datetime import datetime
from datetime import timedelta

class SolarEvent:
    
    def __init__(self, *args, eventtype = None, classification = None, classification_method = None, **kwargs):
        ''' Creates a new SolarEvent object. '''
        
        # Checking to make sure that args has a MapSequence object.
        sequence_exists = False
        if isinstance(args[0], mp.MapSequence) :
            sequence_exists = True
        if sequence_exists == False:
            raise ValueError('Solar Events expect pre-constructed Sunpy MapSequence objects for the first parameter. If you have a map that you would like to pass over, define it as a sequence by calling `sequence = True` when creating the map.')

        self.mapsequence = args[0]
        
        # Currently supported solar events
        self.list_of_events = ['CoronalHoles', 
                               'ActiveRegions', 
                               'HMI_Magnetograms', 
                               'Unspecified_Event']
        
        
        self.list_of_classifications = ['CoronalHolesBoundaries', 
                                        'CoronalLoops',
                                        'Magnetogram_Boundaries',
                                        'Unspecified_Event_Classification']
        
        self.list_of_event_classification_methods = ['CHIMERA', 'CHARM', 
                                                     'TRACECoronalLoops',
                                                     'GaussianContours']
        
        self.full_list_of_classification_methods = ['CHIMERA', 'CHARM', 
                                                    'TRACECoronalLoops',
                                                    'GaussianContours',
                                                    'SolarMFDFA', 'Unspecified']
        
        self.eventtype = eventtype
        self.classification = classification
        self.classification_method = classification_method
        
    def sub_maps(self, bottom_left, top_right):
        ''' Creates a submap for each map given in args
            then returns all of the maps as a mp.MapSequence. '''
        sub_map_list = []
        # Go through the list of maps and create submaps from them, adding all maps into a list.
        for i, arg in enumerate(self.mapsequence):
            sub_map_list.append(arg.submap(bottom_left = bottom_left,
                                           top_right = top_right))
        
        map_sequence_submaps = mp.Map(sub_map_list, sequence = True)
        
        derotated_submap_sequence = coalignment.mapsequence_coalign_by_rotation(map_sequence_submaps)
        
        return derotated_submap_sequence
    
    
    
    def date_match(self, mapseq, maptocompare, index_to_consider = None, time_change = None):
        '''
        Compares temporal differences between a sunpy map sequence and a generic map.

        Parameters
        ----------
        mapseq : `sunpy.map.MapSequence`
            A `~sunpy.map.MapSequence` of shape ``(ny, nx, nt)``, where ``nt`` is the number of
            layers in the `~sunpy.map.MapSequence`. ``ny`` is the number of pixels in the
            "y" direction, ``nx`` is the number of pixels in the "x" direction. This map seq will
            be used to search for a matching date.
        
        maptocompare : `sunpy.map.GenericMap`
           A `~sunpy.map.MapSequence` of shape ``(ny, nx, nt)``, where ``nt`` is the number of
           layers in the `~sunpy.map.MapSequence`. ``ny`` is the number of pixels in the
           "y" direction, ``nx`` is the number of pixels in the "x" direction. This map seq will
           be used to find an appropriate map to overlay onto the mapseq.
           
        index_to_consider : int
            The layer, nt, of the mapseq to consider when matching times. If none is given,
            the middle index of mapseq is considered.
            
        time_change : `datetime.timedelta` object
            A small correction term to the mapseq`s date that accounts for small 
            temporal differences.

        Returns
        -------
        map_ : `sunpy.map.GenericMap`
            The map where the date matches. Or false if there is no match.

        '''
        # Sometimes telescopes mismatch times, so we add a small correction to 
        # ensure that the events match
        if time_change == None:
            time_change = timedelta(minutes = 10)
        # If index_to_consider was given, 
        if index_to_consider == None :
            index_to_consider = int(len(mapseq)/2)
        if maptocompare.meta.get('telescop') == 'SDO/AIA':
            maptocompare_date = datetime.fromisoformat(maptocompare.meta.get('date-obs')[:-8])
            
        if maptocompare.meta.get('telescop') == 'SDO/HMI' or maptocompare.meta.get('telescop') == None:
            maptocompare_date = datetime.fromisoformat(maptocompare.meta.get('date-obs')[:-4])
            
        
        # TODO make this applicable to all types of input data/sattelites
        previousdate = datetime.fromisoformat(mapseq[index_to_consider-1].meta.get('date-obs')[:-4]) - time_change
        nextdate = datetime.fromisoformat(mapseq[index_to_consider+1].meta.get('date-obs')[:-4]) + time_change

        
        if previousdate <= maptocompare_date and maptocompare_date <= nextdate:
            print(maptocompare_date)
            # First re-check to see if there is a matching date
            if previousdate == maptocompare_date:
                return mapseq[index_to_consider-1], index_to_consider - 1
            if nextdate == maptocompare_date:
                return mapseq[index_to_consider+1], index_to_consider + 1
            # Return the default otherwise
            return mapseq[index_to_consider], index_to_consider
        
            
        return False


    def plot_matched_events(self, eventlist, fig, ax):
        
        legend_elements = []
        
        # Building the correct map coordinates . . . 
        for idx, mapsequence_list in enumerate(eventlist):
            
            mapsequence = mapsequence_list[0]
            index_list = []
            coords_bot_list = []
            coords_top_list = []
            for i, map_ in enumerate(mapsequence):
                map_match = self.date_match(self.mapsequence, map_)
                
                if map_match != False:
                    
                    mapmatch = map_match[0]
                    index = map_match[1]
                    index_list.append(i)
                    bottom_left = SkyCoord(mapmatch.bottom_left_coord.Tx, 
                                           mapmatch.bottom_left_coord.Ty, 
                                           frame = self.mapsequence[index].coordinate_frame)
                    top_right = SkyCoord(mapmatch.top_right_coord.Tx, 
                                         mapmatch.top_right_coord.Ty, 
                                         frame = self.mapsequence[index].coordinate_frame)
                    coords_bot_list.append(bottom_left)
                    coords_top_list.append(top_right)
                    
                
                
            classification_method = mapsequence_list[1]
            
            if index_list != []:
                index = index_list[0]
            else:
                index = 0
                
            if coords_bot_list != []:
                bottom_left = coords_bot_list[0]
            else:
                bottom_left = self.mapsequence[0].bottom_left_coord
                
            if coords_top_list != []:
                top_right = coords_top_list[0]
            else:
                top_right = self.mapsequence[0].top_right_coord
            
            # Now we build the boundaries with the corresponding map coordinates.
            if classification_method == self.list_of_event_classification_methods[0]:
                
                self.chimera_boundary = mapsequence[index].submap(bottom_left = bottom_left, 
                                                                     top_right = top_right)
                
                boundary_chimera = plt.contourf(self.chimera_boundary._data, 
                                               levels = [0.5,1.5], alpha = 0.85,
                                               cmap = cm.winter)
                legend_elements.append(Line2D([0],[0], color = 'cornflowerblue', lw = 2, label = 'CHIMERA'))

                
            if classification_method == self.list_of_event_classification_methods[1]:
                
                self.charm_boundary = mapsequence[index].submap(bottom_left = bottom_left, top_right = top_right)
                boundary_charm = plt.contourf(self.charm_boundary._data, 
                                             levels = [0.5,1.5], alpha = 0.85,
                                             cmap = cm.cool)
                legend_elements.append(Line2D([0],[0], color = 'aqua', lw = 2, label = 'CHARM'))

        
            if classification_method == self.list_of_event_classification_methods[2]:
                
                self.loop_map = mapsequence[index].submap(bottom_left = bottom_left, top_right = top_right)
                
                trace_loop_contour = plt.contourf(self.loop_map._data, 
                                               levels = [0.5,1.5], alpha = 0.85,
                                               cmap = cm.autumn)
                legend_elements.append(Line2D([0],[0], color = 'orange', lw = 2, label = 'TRACE Loops'))
        
            if classification_method == self.list_of_event_classification_methods[3]:
                # Grab the comment that states the gaussian level . . . 
                GaussianLevel = mapsequence[index].meta.get('comment')[-2:]
                
                self.mag_boundary = mapsequence[index].submap(bottom_left = bottom_left, top_right = top_right)
                
                positive_network = plt.contourf(self.mag_boundary._data, 
                                               levels = [0.5,1.5], alpha = 0.85,
                                               cmap = cm.gray)
                negative_network = plt.contourf(self.mag_boundary._data, 
                                               levels = [-1.5,-0.5], alpha = 0.85,
                                               cmap = cm.seismic)
                legend_elements.append(Line2D([0],[0], color = 'gray', lw = 2, label = 'B LOS (+), Gaussian Level ' + str(GaussianLevel)))
                legend_elements.append(Line2D([0],[0], color = 'white', lw = 2, label = 'B LOS (-), Gaussian Level ' + str(GaussianLevel)))
        ax.legend(handles = legend_elements, fontsize = 10)