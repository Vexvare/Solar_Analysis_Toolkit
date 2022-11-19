# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:14:27 2022

@author : Landon Wells
"""

from .Solar_Events_Class import SolarEvent

class CoronalHole(SolarEvent):
    
    def __init__(self, *args, eventtype = 'Coronal Hole', classification, classification_method, **kwargs):

        super().__init__(*args, eventtype = eventtype, classification = classification, 
                         classification_method = classification_method, **kwargs)