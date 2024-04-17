"""
Automatisation of proteomics data either from 
    - Bands cutting and proteins identification
    - Pull down experiments

-> contaminants removal
-> molecular-weight cut offs
-> comparison between datasets (bands-bands, pulldown-down)
-> plot a volcano plot highlighting the protein.s found in the band.s

Generate and save excel files, saved in the last given path :
    - Bands analysis : 
        . proteins without contaminants and after molecular weight cutoffs, one sheet per band
        . proteins in common between two bands (one sheet for the commons proteins + one sheet for each band with only the proteins specific to the band
    - Pulldown :
        . one sheet per band if protein.s in common between band / pulldown
        . volcano plot, png format

! You need to have the following columns (or adapt the code):
 - 'accession' in all 
 - 'MW' in the bands dataframes
 - 't-test_g1_vs_g2' and 'ratio_g1_vs_g2' in the pulldown dataframe


Lea Masson
19/02/2024 
"""

import tkinter as tk 
from PIL import Image, ImageTk

import base64
import io

from tkinter import messagebox, simpledialog, filedialog, ttk

import app_def as app


import pandas as pd
import itertools
import numpy as np
import matplotlib.pyplot as plt

import auto_prot_def as auto_prot



# initialize window
root = tk.Tk()
root.title('Proteomics cross-results app')
root.config(bg="#9DD4D2")
root.geometry('1300x770')


# create left frame widget : will contain subtitle, description button, requirements button, line separation and team's logo
left_frame = tk.Frame(root,width=200, height=760)
left_frame.grid_propagate(1)
left_frame.pack(side = tk.LEFT, padx=10, pady=10)

    #subtitle
app.display_subtitle(left_frame)

    #description button
description_button = tk.Button(left_frame, text = 'Description', font = ('Arial',12), bd = 3, relief="ridge",
                               command = app.description).grid(row = 3, column = 0, padx=10, pady=15, ipadx=18, ipady=7)

    #requirements button
dependencies_button = tk.Button(left_frame, text = 'Dependencies', font = ('Arial',12), bd = 3, relief="ridge",
                               command = app.dependencies).grid(row = 4, column = 0, padx=10, pady=25,ipadx=10, ipady=7)




# create right frame widget : big square in wich we'll put the code : 

# this frame is in a canvas to be scrollable
cTableContainer = tk.Canvas(root,width=930, height=750)
right_frame = tk.Frame(cTableContainer)
right_frame.pack(side = tk.RIGHT, padx=10, pady=10)

# scrollbar
sbHorizontalScrollBar = tk.Scrollbar(root)
sbVerticalScrollBar = tk.Scrollbar(root)



# Sets up the Canvas, Frame, and scrollbars for scrolling
def createScrollableContainer():
	cTableContainer.config(xscrollcommand=sbHorizontalScrollBar.set,yscrollcommand=sbVerticalScrollBar.set, highlightthickness=0)
	sbHorizontalScrollBar.config(orient=tk.HORIZONTAL, command=cTableContainer.xview)
	sbVerticalScrollBar.config(orient=tk.VERTICAL, command=cTableContainer.yview)

	sbHorizontalScrollBar.pack(fill=tk.X, side=tk.BOTTOM, expand=tk.FALSE)
	sbVerticalScrollBar.pack(fill=tk.Y, side=tk.RIGHT, expand=tk.FALSE)
	cTableContainer.pack(fill=tk.BOTH, side=tk.LEFT, expand=tk.TRUE)
	cTableContainer.create_window(0, 0, window=right_frame, anchor=tk.NW)

createScrollableContainer()



# def start : does all when button is clicked

def start():
    
    try : 
         # 1) how many bands were cut ? -> Enter the number, launch def "get_data" 
        how_many = tk.simpledialog.askinteger(' ', 'How many bands were cut ?')
                # ask for the bands, create 
        bands_dict, prints, bands_path = auto_prot.get_data(how_many)
                # display given info  
        app.display_given_data(prints, right_frame) # row 1 on grid 
        app.updateScrollRegion(cTableContainer, right_frame)
        
        current_row = 2

        # 2) automatic bands treatment -> see functions in auto_proteo_def.py
        clean_df, current_row = auto_prot.bands(bands_dict, bands_path, right_frame, current_row)
        app.updateScrollRegion(cTableContainer, right_frame)

        # 3) automanic pulldown cross-result with bands -> see functions in auto_proteo_def.py
        pulldown = tk.messagebox.askyesno('Pulldown', 'Do you have pulldown results to cross with your bands?')

        if pulldown :
            auto_prot.pulldown_treatment(clean_df, right_frame, current_row)
            app.updateScrollRegion(cTableContainer, right_frame)

        else : 
            ending = tk.Label(right_frame, 'Thanks for using this automated program!').grid(row = 8, column = 0, padx = 30, pady = 10)
            app.updateScrollRegion(cTableContainer, right_frame)
        
    except Exception  as e : 
        messagebox.showinfo('Error!', e) 


# add a start button
start_button = tk.Button(left_frame, text = 'Start', font = ('Arial',14), bd = 3, relief = 'ridge',
                         command = start).grid(row=2, column=0, padx=15, pady=10, ipadx = 10)


# start the app    
root.mainloop()

