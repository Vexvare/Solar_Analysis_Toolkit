# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 21:51:33 2022

@author: Landon Wells

"""

# Class that defines the types of coronal hole boundaries that can be implemented
# ontop of the Coronal Hole object.

# General Module imports
import os
import copy
import sunpy.map as mp
import glob
import sunkit_image.trace as trace
import numpy as np
from astropy import units as u

# Local Module imports
from .ActiveRegions import ActiveRegion
from QOL_Programs.Build_Data_SavePath import buildsavepath_classificationdata

class ActiveRegionCoronalLoops(ActiveRegion):
    '''
    Class defining the loops of a coronal active region. Subclass of ActiveRegion,
    which is a child class of SolarEvent. The first parameter, *args, is assumed
    to be of type `mp.MapSequence` or `mp.GenericMap`.
    '''
    
    def __init__(self, *args, classification = 'CoronalLoops', classification_method = None, **kwargs):
        
        # Currently supported classification methods for a CH boundary.
        self.list_of_classification_methods = ['TRACECoronalLoops']
        self.mapsequence = args[0]
        
        super().__init__(*args, 
                         classification = classification, 
                         classification_method = classification_method, 
                         **kwargs)
        
        self.__dict__.update(**kwargs)
        self.kwargs = kwargs
        
        if self.classification_method == 'TRACECoronalLoops' :
            self.traceloops_maps = self.Trace_Coronal_Loops(args)


    def Trace_Coronal_Loops(self, *args, nsm1 = 3, rmin = 30, lmin = 25, nstruc = 1000, ngap = 0, qthresh1 = 0.0, qthresh2 = 3.0):
        '''
        Function to call to trace coronal loops using the TRACE algorithm, along with saving the
        data to be used in repeated callbacks. This function should be called through
        it`s parent class, `ActiveRegionCoronalLoops` using the classification_method = `CoronalLoops`.
        
        WARNING : If different parameters are used each time the same map is ran, the new map will not
                  overwrite the old map. The old map needs to be deleted if different parameters are used
                  currently. 
        
        Parameters
        ----------
        *args : `mp.MapSequence`
            A mp.MapSequence object containing all of the maps of interest.
        nsm1 : `int`, optional
            Low pass filter boxcar smoothing constant. The default is 3.
        rmin : `int`, optional
            The minimum radius of curvature of the loop to be detected in pixels. The default is 30.
        lmin : `int`, optional
            The length of the smallest loop to be detected in pixels. The default is 25.
        nstruc : `int`, optional
            Maximum limit of traced structures. The default is 1000.
        ngap : `int`, optional
            Number of pixels in the loop below the flux threshold. The default is 0.
        qthresh1 : `float`, optional
            The ratio of image base flux and median flux. All the pixels in the image below
            `qthresh1 * median` intensity value are made to zero before tracing the loops. The default is 0.0.
        qthresh2 : `float`, optional
            The factor which determines noise in the image. All the intensity values between
            `qthresh2 * median` are considered to be noise. The median for noise is chosen
            after the base level is fixed. The default is 3.0.

        Returns
        -------
        mp.Map containing identifying all of the loops on the input map.
        
        '''
        
        savepath = buildsavepath_classificationdata(sunpymap_sequence = self.mapsequence, 
                                                    classification_method = self.classification_method)
        loops_list = []
        for i, maps in enumerate(self.mapsequence):
            # Checking for the loops map . . . if it is not found in 
            # ../Solar_Analysis_Toolkit/Data/Classification_Data then the 
            # TRACE program will be ran.
            filename = 'CoronalLoops_Map_' + str(i) + '.fits'
            file_list = glob.glob(savepath + filename)
            if file_list == []:
                
                print('\nThe TRACE Coronal loops of interest was not found.\nAttempting to find the loops.\n\nCurrent time : ' + maps.meta.get('date-obs'))
    
                # Calling TRACE to identify the loops
                loops = trace.occult2(maps.data, nsm1 = nsm1, rmin = rmin, 
                                      lmin = lmin, nstruc = nstruc, ngap = ngap, 
                                      qthresh1 = qthresh1, qthresh2 = qthresh2)
                
                # Using the AIA meta data as the CHARM meta . . . lazy coding but it works.
                modified_header = copy.deepcopy(maps.meta)
                modified_header['comment'] = 'TRACE LOOPS MASK. nsm1 =' + str(nsm1) +  ', rmin = ' + str(rmin) +  ', lmin = ' + str(lmin) +  ', nstruc = ' + str(nstruc) +  ', ngap = ' + str(ngap) +  ', qthresh1 = ' + str(qthresh1) +  ', qthresh2 = ' + str(qthresh2) +  '.'
                
                # Creating a blank array to draw the loops
                coord_loop_data = np.zeros((maps.data.shape))
                for loop in loops:
                    # convert to array as easier to index `x` and `y` coordinates
                    loop = np.array(loop)
                    #coord_loops = maps.pixel_to_world(loop[:, 0] * u.pixel, loop[:, 1] * u.pixel)
                    for i, loo in enumerate(loop):
                        coord_loop_data[int(loop[i,1]),int(loop[i,0])] = 1
                traceloops_maps = mp.Map([coord_loop_data, modified_header])
                
                loops_list.append(traceloops_maps)
                
        if loops_list != []:
            trace_map = mp.Map(loops_list, sequence = True)
            trace_map.save(savepath + 'CoronalLoops_Map_{index}.fits')
        filename = 'CoronalLoops_Map_*.fits'
        # Finding the list of files that match the CHIMERA save data name.
        file_list = glob.glob(savepath + filename)
        self.trace_map = mp.Map(file_list, sequence = True)
        