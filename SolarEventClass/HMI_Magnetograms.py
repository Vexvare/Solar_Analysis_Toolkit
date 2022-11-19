# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:08:11 2022

@author : Landon Wells
"""
from .Solar_Events_Class import SolarEvent

class HMIMagnetogram(SolarEvent):
    
    def __init__(self, *args, eventtype = 'HMI Magnetogram', classification, classification_method, **kwargs):

        super().__init__(*args, eventtype = eventtype, classification = classification, 
                         classification_method = classification_method, **kwargs)
    
    
    def flux_imbalance(self, *args, **kwargs):
        raise NotImplementedError("This functionality has not yet been implemented.")
