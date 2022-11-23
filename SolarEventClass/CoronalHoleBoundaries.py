# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:58:04 2022

@author : Landon Wells
"""

# Class that defines the types of coronal hole boundaries that can be implemented
# ontop of the Coronal Hole object.

# General Module imports
import os
import copy
import sunpy.map as mp
import glob

# Local Module imports
from .CoronalHoles import CoronalHole
from classification_method_programs import chimerapy as chimera
from classification_method_programs.charmpy import Charm
from QOL_Programs.Build_Data_SavePath import buildsavepath_classificationdata

class CoronalHoleBoundary(CoronalHole):
    
    def __init__(self, *args, classification = 'Boundary', classification_method = None, **kwargs):
        
        # Currently supported classification methods for a CH boundary.
        self.list_of_classification_methods = ['CHIMERA', 'CHARM']
        self.mapsequence = args[0]
        
        super().__init__(*args, 
                         classification = classification, 
                         classification_method = classification_method, 
                         **kwargs)
        
        self.__dict__.update(**kwargs)
        self.kwargs = kwargs
        
        if self.classification_method == 'CHIMERA' :
            self.boundary = self.Chimera(*args)
            
            
        if self.classification_method == 'CHARM' :
            self.boundary = self.Charm(**kwargs)


    def Chimera(self, *args, **kwargs):
        savepath = buildsavepath_classificationdata(sunpymap_sequence = self.mapsequence, 
                                                    classification_method = self.classification_method)
        for i, maps in enumerate(args):
            if maps[0].meta.get('telescop') == 'SDO/HMI' :
                hmimaps = maps
            if maps[0].meta.get('telescop') == 'SDO/AIA' :
                if maps[0].meta.get('wavelnth') == 211 :
                    aia211maps = maps
                if maps[0].meta.get('wavelnth') == 193 :
                    aia193maps = maps
                if maps[0].meta.get('wavelnth') == 171 :
                    aia171maps = maps
                    
        boundary_list = []
        for i, maps in enumerate(aia171maps):
            # Checking for the Chimera map . . . if it is not found in ../Solar_Analysis_Toolkit/Data/Classification_Data then the Chimera program will be ran.
            filename = 'Chimera_Map_' + str(i) + '.fits'
            
            file_list = glob.glob(savepath + filename)
            if file_list == []:
                
                print('\nThe Chimera Coronal Hole boundary of interest was not found.\nAttempting to find the Chimera boundaries.\n\nCurrent time : ' + maps.meta.get('date-obs'))
                
                chimera_boundaries = chimera.Chimera(aia171maps[i], aia193maps[i], aia211maps[i], *args)
    
                # Using the AIA meta data as the CHARM meta . . . lazy coding but it works.
                modified_header = copy.deepcopy(maps.meta)
                modified_header['comment'] = 'CHIMERA MASK.'
                
                chimera_map = mp.Map([chimera_boundaries[2], modified_header])
                
                boundary_list.append(chimera_map)
                
        if boundary_list != []:
            chimera_map = mp.Map(boundary_list, sequence = True)
            chimera_map.save(savepath + 'Chimera_Map_{index}.fits')
        filename = 'Chimera_Map_*.fits'
        # Finding the list of files that match the CHIMERA save data name.
        file_list = glob.glob(savepath + filename)
        self.chimera_map = mp.Map(file_list, sequence = True)



    def Charm(self, **kwargs):
        savepath = buildsavepath_classificationdata(sunpymap_sequence = self.mapsequence, 
                                                    classification_method = self.classification_method)
        intensitylevel = self.kwargs['intensitylevel']
        boundary_list = []
        for i, maps in enumerate(self.mapsequence):
            # Checking for the charm map . . . if it is not found in ../Solar_Analysis_Toolkit/Data/Classification_Data then the charm program will be ran.
            filename = 'Charm_Map_' + str(i) + '_Intensitylevel_' + str(intensitylevel) + '.fits'
            file_list = glob.glob(savepath + filename)
            if file_list == []:
                
                print('\nThe Charm Coronal Hole boundary of interest was not found.\nAttempting to find the Charm boundaries.\n\nCurrent time : ' + maps.meta.get('date-obs'))
    
                charm_boundaries = Charm(maps, **kwargs)
    
                # Using the AIA meta data as the CHARM meta . . . lazy coding but it works.
                modified_header = copy.deepcopy(maps.meta)
                modified_header['comment'] = 'CHARM MASK.'
                
                charm_map = mp.Map([charm_boundaries[2], modified_header])
                boundary_list.append(charm_map)
                
        if boundary_list != []:
            charm_map = mp.Map(boundary_list, sequence = True)
            charm_map.save(savepath + 'Charm_Map_{index}_Intensitylevel_' + str(intensitylevel) + '.fits')
        filename = 'Charm_Map_*_Intensitylevel_' + str(intensitylevel) + '.fits'
        # Finding the list of files that match the CHIMERA save data name.
        file_list = glob.glob(savepath + filename)
        self.charm_map = mp.Map(file_list, sequence = True)




        
        
        
        
        
        
        
        
        
        

        
    

        
        
        
        
        
        
        
        
        
        
