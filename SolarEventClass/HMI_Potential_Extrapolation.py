# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 14:54:33 2022

@author : Landon Wells
"""

# General Module imports
import os
import sunpy.map as mp
from astropy import units as u
from mayavi import mlab

# Local Module imports
from .HMI_Magnetograms import HMIMagnetogram
from QOL_Programs.Build_Data_SavePath import buildsavepath_submapclassificationdata
from classification_method_programs.solarbextrapolation.extrapolators import PotentialExtrapolator
from classification_method_programs.solarbextrapolation.visualisation_functions import visualise
from classification_method_programs.solarbextrapolation.map3dclasses import Map3D

class HMIMagnetogramExtrapolation(HMIMagnetogram):
    
    def __init__(self, *args, classification = 'Magnetogram Extrapolation', classification_method = None, **kwargs):
        
        # Currently supported classification methods for a Magnetogram potential extrapolation.
        self.list_of_classification_methods = ['Solarbextrapolation']
        
        self.mapsequence = args[0]
        
        super().__init__(*args, 
                         classification = classification, 
                         classification_method = classification_method,
                         **kwargs)
        
        self.__dict__.update(**kwargs)
        self.kwargs = kwargs
    
        if self.classification_method == 'Solarbextrapolation':
            self.solarbextrapolation(**kwargs)
    
    
    
    def solarbextrapolation(self, **kwargs):
        self.hmi_submap = kwargs.get('hmi_submap')
        self.mapbase = kwargs.get('map_base')
        self.zshape = kwargs.get('zshape')
        self.xrange = kwargs.get('xrange')
        self.yrange = kwargs.get('yrange')
        self.zrange = kwargs.get('zrange')
        self.resample_divide = kwargs.get('resample_divide')
        
        savepath = buildsavepath_submapclassificationdata(sunpymap_sequence = self.mapsequence, sunpysubmap_sequence = self.hmi_submap, classification_method = self.classification_method)
        true_savepath = savepath + 'Solarbextrapolator_resample_div_' + str(self.resample_divide) + '_zshape_' + str(self.zshape) + '_zrange_' + str(self.zrange[0]/u.arcsec) + '_' + str(self.zrange[-1]/u.arcsec) + '_arcsec' + '_Bxyz.npy'
            
        if not os.path.isfile(true_savepath):
            print('The potential extrapolation of interest was not found.\nAttempting to create the potential map . . . ')
                
            aPotExt = PotentialExtrapolator(map_magnetogram = self.hmi_submap[0], 
                                                zshape = self.zshape, 
                                                zrange = self.zrange, 
                                                xrange = self.xrange, 
                                                yrange = self.yrange)
            aMap3D = aPotExt.extrapolate()
            Map3D.save(aMap3D, filepath = true_savepath)
        # Load the results.
        aMap3D = Map3D.load(true_savepath)
        
        # We now plot the LOS magnetic field that is produced from the extrapolator above.
        # These will be plotted using Mayavi instead of matplotlib. If you do not have mayavi
        # installed, it will not work. 
        
        # Visualise the 3D vector field
        visualise(aMap3D, 
                  boundary = self.mapbase, 
                  scale=1.0*u.Mm, 
                  boundary_unit=1.0*u.arcsec, 
                  show_boundary_axes=False, 
                  show_volume_axes=True, 
                  debug=True
                  )


        mlab.show()



