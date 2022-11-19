# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 17:33:30 2022

@author: Landon Wells

Purpose : To allow the user to input info using different GUI methods. The GUIs
          are defined in the following files : GUI_Enter_Information.py, in which 
          the user inputs information based on the type of data they want, and
          GUI_FilePath.py in which the user points to the directory where their
          data is stored.
"""
# General Module Imports
import tkinter as tk
from tkinter import messagebox

# Local Module Imports
from .GUI_Enter_Information import gather_info
from .GUI_FilePath import ObtainMapFiles
from .GUI_FilePath import ObtainNpDataFiles

class choose_gui:
    
    def __init__(self, parentone):
        
        self.parentone = parentone
        self.parentone.protocol('WM_DELETE_WINDOW', self.on_close)
        self.gui_list = ['Enter Info', 'Choose Sunpy File Path', 'Chose numpy Data File Paths']
        self.om_vari = tk.StringVar(self.parentone)
        self.om_vari.set(self.gui_list[0])
        self.gui_label = tk.Label(self.parentone, text = 'Choose a way to gather data.')
        self.gui_label.grid(column = 0, row = 0)
        self.omm = tk.OptionMenu(self.parentone, self.om_vari, *self.gui_list)
        self.omm.grid(column = 0, row = 1)
        
        # Find info
        self.buttonfind = tk.Button(self.parentone, text='Click to get info', command = self.getting_info)
        self.buttonfind.grid(column = 2, row = 1)
        
        # Submit info
        self.buttonsubmit = tk.Button(self.parentone, text='Click to confirm info', command = self.confirm_info)
        self.buttonsubmit.grid(column = 3, row = 1)

    def getting_info(self, *args):
        
        self.gathermethod = self.om_vari.get()
        
        if self.gathermethod == self.gui_list[0]:
            self.getinfo = gather_info().info
            
            self.telescope = self.getinfo[0]
            self.observable = self.getinfo[1]
            self.wavelength = self.getinfo[1]
            self.starttime = self.getinfo[2]
            self.endtime = self.getinfo[3]
            self.jsocemail = self.getinfo[4]
            
        if self.gathermethod == self.gui_list[1]:
            self.getinfo = ObtainMapFiles().map_sequence
            
            self.allmapmeta = self.getinfo.all_meta()
            self.telescope = self.getinfo[0].meta.get('telescop')
            self.wavelength = self.getinfo[0].meta.get('wavelnth')
            #self.observable = self.getinfo[0].meta.get('') #TODO check hmi maps on this
            self.starttime = self.getinfo[0].meta.get('date-obs')
            self.endtime = self.getinfo[-1].meta.get('date-obs')
            


        if self.gathermethod == self.gui_list[2]:
            self.getinfo = ObtainNpDataFiles().listoffiles

        return self.getinfo
    
    def confirm_info(self, *args):
        if messagebox.askokcancel('Finished?', 'Are you sure?'):
            self.parentone.destroy()
            print('Finished gathering info.')
        
            
    def on_close(self, *args):
        if messagebox.askokcancel('Quit', 'Quit?'):
            self.parentone.destroy()
            raise ValueError('The program was terminated.')

def getting_info():
    '''
    The method to call when wanting to obtain flexable data information.
    '''
    roothome = tk.Tk()
    info = choose_gui(roothome)
    roothome.mainloop()
    #print(info.satelitte, info.observable, info.starttime, info.endtime, info.jsocemail)
    return info
