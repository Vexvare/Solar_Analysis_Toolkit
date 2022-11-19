# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:10:01 2022

@author : Landon Wells
"""
# General Module imports
import os
import copy
import sunpy.map as mp
import glob

# Local Module imports
from .HMI_Magnetograms import HMIMagnetogram
from QOL_Programs.Build_Data_SavePath import buildsavepath_classificationdata
from classification_method_programs import magnetogram_gaussian_contours as gauss_countours


class HMIMagnetogramBoundary(HMIMagnetogram):
    
    def __init__(self, *args, classification = 'Magnetogram Boundary', classification_method = None, **kwargs):
        
        # Currently supported classification methods for a Magnetogram boundary.
        self.list_of_classification_methods = ['GaussianContours']
        
        self.mapsequence = args[0]
        
        super().__init__(*args, 
                         classification = classification, 
                         classification_method = classification_method,
                         **kwargs)
        
        self.__dict__.update(**kwargs)
        self.kwargs = kwargs
    
        if self.classification_method == 'GaussianContours' :
            self.boundary = self.GaussianContours(**self.kwargs)
            
            
    def GaussianContours(self, GaussianLevel):
        savepath = buildsavepath_classificationdata(sunpymap_sequence = self.mapsequence, classification_method = self.classification_method)
        GaussianLevel = self.kwargs['GaussianLevel']
        if GaussianLevel == None:
            raise ValueError('You must pass over the GaussianLevel that you would like to contour over.')
        boundary_list = []
        for i, maps in enumerate(self.mapsequence):
            filename = 'Gauss_Map_' + str(i) + '_GaussianLevel_' + str(GaussianLevel) + '.fits'
            # Checking to see if the files exist . . . 
            file_list = glob.glob(savepath + filename)
            if file_list == []:
                print('The magnetogram boundary of interest was not found.\nAttempting to find the magnetogram boundaries.\nCurrent time : ' + maps.meta.get('date-obs'))
                
                gauss_countour_map = gauss_countours.gauss_countours(sunpymap = maps, # self.mapsequence can also be used here if one were to find the mag boundaries in multiple maps
                                                    gaussianlevel = GaussianLevel)
                
                # Using the HMI meta data as the megnetogram gaussian contour meta . . . lazy coding but it works.
                modified_header = copy.deepcopy(maps.meta)
                modified_header['comment'] = 'Magnetogram boundaries mask, GaussianLevel_' + str(GaussianLevel)
                magnetogram_boundary_map = mp.Map([gauss_countour_map, modified_header])
                boundary_list.append(magnetogram_boundary_map)
        if boundary_list != []:
            magnetogram_boundary_maps = mp.Map(boundary_list, sequence = True)
            magnetogram_boundary_maps.save(savepath + 'Gauss_Map_{index}_GaussianLevel_' + str(GaussianLevel) + '.fits')
        
        filename = 'Gauss_Map_*_GaussianLevel_' + str(GaussianLevel) + '.fits'
        # Finding the list of files that match the CHIMERA save data name.
        file_list = glob.glob(savepath + filename)
        self.magnetogram_boundary_maps = mp.Map(file_list, sequence = True)