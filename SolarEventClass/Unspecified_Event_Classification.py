# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:35:40 2022

@author: Landon Wells
"""

# General Module imports
import numpy as np

# Local Module imports
from .Unspecified_Event import UnspecifiedEvent
from classification_method_programs.SolarMFDFA import SolarMFDFA
from QOL_Programs.GUI_FilePath import ObtainMapFiles
from QOL_Programs.GUI_FilePath import ObtainNpDataFiles
from QOL_Programs.Build_Data_SavePath import buildsavepath_classificationdata

class UnspecifiedClassification(UnspecifiedEvent):
    
    def __init__(self, *args, classification = 'Unspecified', classification_method = None, **kwargs):
        
        # Currently supported classification methods for unspecified events.
        self.list_of_classification_methods = ['SolarMFDFA', 'Unspecified']
        
        self.mapsequence = args[0]
        
        super().__init__(*args, 
                         classification = classification, 
                         classification_method = classification_method, 
                         **kwargs)
        
        self.__dict__.update(**kwargs)
        self.kwargs = kwargs
        
        if self.classification_method == 'SolarMFDFA' :
            self.mfdfamap = self.Solar_MFDFA(**self.kwargs)
            
        
    def Solar_MFDFA(self, qrange, orderfit, scales):
        timeseries = self.mapsequence.as_array()
        self.mfdfamap, self.qrange, self.scales, self.orderfit = SolarMFDFA(timeseries, qrange, orderfit, scales)
        return self.mfdfamap
        
    
    def save_MFDFA_map(self):
        '''
        Function to save the 4-D MFDFA map information as a np array, along
        with the first map and last map of the timeseries. These will be saved in 
        ../Data/Classification_Data/SolarMFDFA/start_---_end_---/wavelength_---/lowerleft_--AS_upperright_--AS/ . . . 
        Along side the MFDFA data, we also store the information of the qvalues, 
        and the scales used for the analysis. These are all saved as .np files 
        besides the sunpy maps, and must all be called on for future reference.
        '''
        
        # Gathering information to and cleanly store the information. 
        path = buildsavepath_classificationdata(sunpymap_sequence = self.mapsequence, 
                                                classification_method = self.classification_method)
        # Saving the data
        # MFDFA Data
        np.save(path + 'MFDFA_DATA.npy', self.mfdfamap)
        np.save(path + 'qrange.npy', self.qrange)
        np.save(path + 'scales.npy', self.scales)
        
        # We save the cut map data for future callback.
        self.mapsequence.save(path + 'map_{index}.fits')
        
    def load_MFDFA_map():
        '''
        Function that loads the previously saved data, along with the time series
        that was used to create the data.
        '''
        
        print('Please Identify the Solar MFDFA data.')
        info = ObtainNpDataFiles()
        
        data_sequence = info.listoffiles
        submap_sequence = info.map_sequence
        datafilepath = info.filepath
        
        mfdfa_map = np.load(data_sequence[0], allow_pickle=True)
        qrange = np.load(data_sequence[1], allow_pickle=True)
        scales = np.load(data_sequence[2], allow_pickle=True)
        
        #print('Identify the Sunpy Map Data that was used to create the Solar MFDFA data.')
        #sunpymaps = ObtainMapFiles().map_sequence
        
        return submap_sequence, mfdfa_map, qrange, scales, datafilepath