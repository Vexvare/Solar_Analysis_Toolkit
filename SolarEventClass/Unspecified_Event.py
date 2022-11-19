# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:29:54 2022

@author: Landon Wells
"""

# Sometimes, some calculations / analysis methods will not be over any specific event 
# (or perhaps multiple events  . . .  e.g. comparison of the intensity between 
# two events, or over an entire map). These will be called 'Unspecified' events
# and the list of classification methods will be a 'GeneralClassification' class.

from .Solar_Events_Class import SolarEvent

class UnspecifiedEvent(SolarEvent):
    ''' Creates a new UnspecifiedEvent, which is a instance of SolarEvent. '''
    def __init__(self, *args, eventtype = 'UnspecifiedEvent', classification, classification_method, **kwargs):

        super().__init__(*args, eventtype = eventtype, classification = classification, 
                         classification_method = classification_method, **kwargs)
        
