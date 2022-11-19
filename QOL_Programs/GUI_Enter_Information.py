# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 11:29:21 2022

@author : Landon Wells

Purpose : To allow users to input information about data they wish to download.
"""

import os
import sys
import glob
from tkinter import *
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from astropy.io import fits
import sunpy.map as mp

# https://sdac.virtualsolar.org/cgi/show_details?keyword=PHYSOBS


class info_gather:
    
    def __init__(self,parent):
        
        self.parent = parent
        self.parent.protocol('WM_DELETE_WINDOW', self.on_close)
        self.satellite_list = ['HMI', 'AIA']
        self.wavelength_list = [1700, 4500, 1600, 304, 171, 193, 211, 335, 94, 131]
        self.phys_obs_list = ['LOS_MAGNETIC_FIELD']
        self.om_variable = tk.StringVar(self.parent)
        self.om_variable.set(self.satellite_list[0])
        self.satellite_label = tk.Label(self.parent, text = 'List of satellites.')
        self.satellite_label.grid(column = 0, row = 0)
        self.om = tk.OptionMenu(self.parent, self.om_variable, *self.satellite_list)
        self.om.grid(column = 0, row = 1)
        self.satellite_label = tk.Label(self.parent, text = 'List of observables.')
        self.satellite_label.grid(column = 1, row = 0)
        self.om_variable.trace('w', self.update_options)
        # Start Time
        self.start_label = tk.Label(self.parent, text = 'Start Time')
        self.start_label.grid(column = 0, row = 2)
        self.entry1 = tk.Entry(self.parent)
        self.entry1.grid(column = 0, row = 3)
        # End Time
        self.end_label = tk.Label(self.parent, text = 'End Time')
        self.end_label.grid(column = 1, row = 2)
        self.entry2 = tk.Entry(self.parent)
        self.entry2.grid(column = 1, row = 3)
        # jsoc Email address
        self.email_label = tk.Label(self.parent, text = 'JSOC Email')
        self.email_label.grid(column = 2, row = 2)
        self.entry3 = tk.Entry(self.parent)
        self.entry3.grid(column = 2, row = 3)
        # Label help for time format
        self.help_label = tk.Label(self.parent, text = 'Time Format Example : 2017-04-16T12:00:00')
        self.help_label.grid(column = 1, row = 4)
        # Submit info
        self.submitbutton = tk.Button(self.parent, text='Click to Submit Info', command = self.get_info)
        self.submitbutton.grid(column = 2, row = 0)
        

    def update_options(self, *args):
        if self.om_variable.get() == self.satellite_list[0]:
            self.om_var = tk.StringVar(self.parent)
            self.om_var.set('')
            self.om1 = tk.OptionMenu(self.parent, self.om_var, [])
            self.om_var.set(self.phys_obs_list[0])
            self.om1 = tk.OptionMenu(self.parent, self.om_var, *self.phys_obs_list)
            self.om1.grid(column = 1, row = 1)
            
        if self.om_variable.get() == self.satellite_list[1]:
            self.om_var = tk.IntVar(self.parent)
            self.om_var.set(self.wavelength_list[0])
            self.om1 = tk.OptionMenu(self.parent, self.om_var, [])
            self.om_var.set(self.wavelength_list[0])
            self.om1 = tk.OptionMenu(self.parent, self.om_var, *self.wavelength_list)
            self.om1.grid(column = 1, row = 1)

    def get_info(self, *args):
        
        self.info = [self.om_variable.get(), self.om_var.get(), self.entry1.get(), self.entry2.get(), self.entry3.get()]
        if self.om_var.get() == '':
            raise ValueError('Please input a physical observable or wavelength of interest.')
        if self.entry1.get() == '':
            raise ValueError('Please input a start time for the data of interest.')
        if self.entry2.get() == '':
            raise ValueError('Please input a end time for the data of interest.')
        if self.entry3.get() == '':
            raise ValueError('Please input a JSOC email address.')
        self.instrument = self.info[0]
        self.wavelength = self.info[1]
        self.observable = self.info[1]
        self.starttime = self.info[2]
        self.endtime = self.info[3]
        self.jsoc_email = self.info[4]
        
        
        if messagebox.askokcancel('Quit', 'Quit?'):
            self.parent.destroy()
        return self.info
            
    def on_close(self, *args):
        if messagebox.askokcancel('Quit', 'Quit entering information?'):
            self.parent.destroy()
            raise ValueError('The program was terminated.')
    
    
def gather_info():
    root = tk.Tk()
    info = info_gather(root)
    root.mainloop()
    return info






