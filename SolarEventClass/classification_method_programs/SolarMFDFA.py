# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 00:40:49 2022

@author : Landon Wells

Purpose : To create a Solar MFDFA map.

Inputs : timeseries - numpy 3d array with dimensions [ycoords,xcoords,image number].
         qrange - The q values with which the MFDFA analysis will use. If none is 
                  given then the q's will be [0,2] (0 because the MFDFA program was patched
                  to force the 0th moment).
         orderfit - The order of the polynomial that will be used to calculate the
                    variance across the detrended series. If none is given the order will be 1.
         scales - The scales over which the analysis will be applied. If none are given,
                  the scales will range from 4 images to the length of the time series / 3 in
                  integer values scaling in the order of 2^(1/3). 
            
Output : mfdfa_fluctuations_map - np 4d array with dimensions [ycoords,xcoords,[fluctuation functions]]
                                  where [fluctuation functions] is an array containing the values
                                  of the fluctuations across the moments for all scales.
                                  
"""

import numpy as np
from math import e
from classification_method_programs.MFDFA import MFDFA

def SolarMFDFA(timeseries: np.ndarray, 
               qrange: np.ndarray = 2, 
               orderfit: int = 1, 
               scales: np.ndarray = None):
    '''
    Function that calculates the MFDFA method over temporal data. Of particular
    interest is solar data, although any data can be used, so long as it is a 
    time series of shape [x,y,t].
    
    Parameters
    ----------
    timeseries - np.ndarray
        The time series over which the MFDFA analysis will be applied. Inputs are
         np.3Darray objects of shape [x,y,t], where x and y are the coordinates
        and t is the span of the time series, in temporal pixel number.
    qrange - np.ndarray 
        The q range that is to be considered over the data. These are the moments
        which which the MFDFA will calculate over. If none are given, q=2, the
        DFA moment is assumed (Along with q=0 due to a `bug` in the parent function).
    orderfit - int
        
    
    '''
    if scales == None:
        # Parameters to set the scales if none are given. Values are in temporal pixels.
        minval = 4.0
        maxval = len(timeseries[0,0,:])/3
        num_scales = 100
        # Using logarithmically spaced scales
        scales = np.logspace(np.log10(minval), np.log10(maxval), num = num_scales, dtype = int )


    mfdfa_fluctuations_map = np.zeros([len(timeseries[:,0,0]), len(timeseries[0,:,0]), len(scales), len(qrange)])
    area = len(timeseries[0,:,0])*len(timeseries[:,0,0])
    xcounter = 0
    ycounter = 0 
    for xc in range (len(timeseries[0,:,0])):
        xcounter += 1
        for yc in range (len(timeseries[:,0,0])):
            ycounter += 1
            area_done = (xcounter)*(ycounter)
            if area % area_done == 0:
                percent_complete = int(area_done/area * 100)
                print(str(percent_complete) + ' % completed.')
            time_series = timeseries[yc,xc,:]
            scales, fluct = MFDFA(time_series, scales, q = qrange, order = orderfit, stat = False, extensions = {'eDFA':False})
            mfdfa_fluctuations_map[yc,xc,:,:] = fluct

    return mfdfa_fluctuations_map, qrange, scales, orderfit
