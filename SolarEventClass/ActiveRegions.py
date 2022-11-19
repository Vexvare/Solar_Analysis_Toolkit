# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:47:37 2022

@author: Landon Wells
"""

from .Solar_Events_Class import SolarEvent

class ActiveRegion(SolarEvent):
    '''
    Subclass of SolarEvent.
    '''
    
    def __init__(self, *args, eventtype = 'Active Region', classification, classification_method, **kwargs):

        super().__init__(*args, eventtype = eventtype, classification = classification, 
                         classification_method = classification_method, **kwargs)
        