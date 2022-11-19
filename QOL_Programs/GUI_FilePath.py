# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 12:04:10 2022

@author : Landon Wells

Purpose : To call upon solar or numpy data using GUI methods.
"""
# General Module Imports
import os
import glob
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import sunpy.map as mp
from os.path import expanduser

class obtainmapfiles:
    
    def __init__(self, parent):
        
        self.parent = parent
        self.parent.protocol('WM_DELETE_WINDOW', self.on_close)
        self.button = Button(self.parent, text='Click to find solar .fits maps.', command = self.find_map_files)
        self.button.grid(column = 2, row = 0)

    def on_close(self, *args):
        if messagebox.askokcancel('Quit', 'Quit?'):
            self.parent.destroy()
            raise ValueError('The program was terminated.')
            
    def find_map_files(self, *args):
        self.filepath = filedialog.askdirectory(initialdir = expanduser('~'), 
                                       title = 'Choose the directory that contains the maps of interest.')
        self.parent.destroy()

        self.listoffiles = sorted( glob.glob(self.filepath + '\\*.fits'))
        if self.listoffiles == []:
            raise RuntimeError('No files were found.')
        print('\nFiles have been found.')
        self.map_sequence = mp.Map( [ mp.Map(str(filename)) for filename in self.listoffiles ], sequence = True)
        
        
        
def ObtainMapFiles():
    root_map = Tk()
    info = obtainmapfiles(root_map)
    root_map.mainloop()
    return info



class obtainnumpyfiles:
    
    def __init__(self, parent):
        
        self.dataparent = parent
        self.dataparent.protocol('WM_DELETE_WINDOW', self.on_close)
        self.button = Button(self.dataparent, text='Click to find .npy data.', command = self.find_data_files)
        self.button.grid(column = 2, row = 0)

    def on_close(self, *args):
        if messagebox.askokcancel('Quit', 'Quit?'):
            self.dataparent.destroy()
            raise ValueError('The program was terminated.')
            
    def find_data_files(self, *args):
        self.filepath = filedialog.askdirectory(initialdir = expanduser('~'), 
                                                    title = 'Choose the directory of the data of interest.')
        self.dataparent.destroy()
        
        self.listoffiles = sorted( glob.glob(self.filepath + '\\*.npy'))
        
        if self.listoffiles == []:
            raise RuntimeError('No files were found.')
        print('\nFiles have been found.')
        self.listoffitsfiles = sorted( glob.glob(self.filepath + '\\*.fits'))
        
        if self.listoffitsfiles != []:
            self.map_sequence = mp.Map( [ mp.Map(str(filename)) for filename in self.listoffitsfiles ], sequence = True)
            
def ObtainNpDataFiles():
    root_data = Tk()
    info = obtainnumpyfiles(root_data)
    root_data.mainloop()
    return info
























# def ObtainMapFiles():
#     '''
#     To obtain Sunpy map files from a user chosen directory. Loads in all files 
#     from that directory with the .fits extension, and outputs them as a Sunpy 
#     mapsequence if they are multiple maps, or a sunpy generic map if it is just
#     a single map. (mapsequence class has this functionality built in it).
#     '''
#     global sequence, map_sequence, destroy
#     sequence = None
#     map_sequence = None
#     destroy = False
#     def FindMapFiles():
    
#         cwd = os.getcwd()
        
#         def openMapFile():
#             global map_sequence, filepath
#             filepath = filedialog.askdirectory(initialdir = cwd, 
#                                            title = 'Choose the directory of the maps of interest.')
#             window.destroy()
    
#             listoffiles = sorted( glob.glob(filepath + '\\*.fits'))
#             map_sequence = mp.Map( [ mp.Map(str(filename)) for filename in listoffiles ], sequence = True)
        
#         def on_close():
#             if messagebox.askokcancel('Quit', 'Quit Looking for Files?'):
#                 global destroy
#                 destroy = True
#                 window.destroy()
#                 raise ValueError('The program was terminated.')
                
#         while destroy is False and map_sequence is None:
#             try:
#                 window = Tk()
#                 window.protocol('WM_DELETE_WINDOW', on_close)
#                 button = Button(text = 'Choose the directory of the maps of interest.', command = openMapFile)
#                 button.pack()
#                 window.mainloop()
#             except:
#                 return None
        
#         if destroy is False:
#             return map_sequence 
#         else: return None

#     while sequence is None and destroy is False:
#         try:
#             sequence = FindMapFiles()
#         except:
#             pass
#         if sequence is None : 
#             print('\nNo .FITS files were found.\nTry looking again.')
#         if sequence is not None : 
#             print('\nFiles have been found.')
#     return sequence, filepath
# def ObtainNpDataFiles():
#     '''
#     To obtain numpy data files from a user chosen directory. Finds in all files 
#     from that directory with the .npy extension and returns the files as a list 
#     of file paths.
#     '''
#     global destroy, listoffiles, listofmapfiles, map_sequence
#     destroy = False
#     listoffiles = None
#     listofmapfiles = None
#     map_sequence = None
#     def FindDataFiles():
    
#         cwd = os.getcwd()
        
#         def openDataFile():
#             global filepath, listoffiles, listofmapfiles, map_sequence
#             filepath = filedialog.askdirectory(initialdir = cwd, 
#                                            title = 'Choose the directory of the data of interest.')
#             window.destroy()
    
#             listoffiles = sorted( glob.glob(filepath + '\\*.npy'))
            
#             listofmapfiles = sorted( glob.glob(filepath + '\\*.fits'))
#             if listofmapfiles is not None:
#                 map_sequence = mp.Map( [ mp.Map(str(filename)) for filename in listofmapfiles ], sequence = True)
            
#         def on_close():
#             if messagebox.askokcancel('Quit', 'Quit Looking for Files?'):
#                 global destroy
#                 destroy = True
#                 window.destroy()
#                 raise ValueError('The program was terminated.')
                
#         while destroy is False and listoffiles is None:
#             try:
#                 window = Tk()
#                 window.protocol('WM_DELETE_WINDOW', on_close)
#                 button = Button(text = 'Choose the directory of the data of interest.', command = openDataFile)
#                 button.pack()
#                 window.mainloop()
#             except:
#                 return None
        
#         if destroy is False:
#             return listoffiles 
#         else: return None

#     while listoffiles is None and destroy is False:
#         try:
#             listoffiles = FindDataFiles()
#         except:
#             pass
#             if listoffiles is None : 
#                 print('\nNo .npy files were found.\nTry looking again.')
#             if listoffiles is not None : 
#                 print('\nFiles have been found.')
#     if listofmapfiles is None:
#         return listoffiles, filepath
#     if map_sequence is not None:
#         return listoffiles, map_sequence, filepath
