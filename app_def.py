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
from tkinter import messagebox
from PIL import Image, ImageTk
import base64
import io



def display_subtitle(left_frame):
    subt_frame = tk.Frame(left_frame, width=220, height=50)
    subt_frame.grid(row =1, column =0, padx=10, pady=10)
    subtitle = tk.Label(subt_frame, text = "Automated processing of proteomics data from protein gel migration band cutting, which can be associated with a pull-down.",
                        justify = 'center', font=('Arial',12), wraplength=210)
    subtitle.grid(row=0, column=0, padx=15, pady=10)  



def display_given_data(prints, right_frame):
    message = 'Attention, the first band given is named B1, the second B2, etc.'
    for x in prints :
        if str(x).startswith('B'):
             message += '\n'
        message += '\n' + str(x)

    message = message.replace("(",'').replace(")",'').replace(",",'').replace("'",'')
    dialogs = tk.Label(right_frame, text = message, font = ('Arial',11), justify = 'left').grid(row = 1, column = 0, padx=10, pady=15, ipadx=18, ipady=7)
    
    

def updateScrollRegion(cTableContainer, right_frame):
	cTableContainer.update_idletasks()
	cTableContainer.config(scrollregion=right_frame.bbox())

    
    
def description():  
    description = '''This program helps you cross-results your excels from protein migration band's identification and pull down to 

From your excel files, it : 
    → Removed contaminants
    → Applies given molecular weight cut offs
    → Compares proteins between datasets (bands-bands, pulldown-bands)
    → Plots a volcano plot from the pull-down results, highlighting the protein.s found in the band.s

    → Generates and save excel files in a given path :
        → Bands analysis results :
            . Proteins without contaminants and after molecular weight cutoffs, one sheet per band
            . Proteins in common between two bands (one sheet for the commons proteins + one sheet for each band with only the proteins specific to the band)
        → Pulldown-bands cross results :
            . One sheet per band if protein.s in common between band / pulldown
            . Volcano plots (png format), and a volcano plot highlighting the proteins with the 15 best ratio

⚠ You need to have the following columns in your excels (or adapt the code) :
     → 'accession' in all with the accessions 
     → 'MW' in the bands' dataframes with the molecular weights in Dalton
     → 't-test_g1_vs_g2' and 'ratio_g1_vs_g2' in the pull down dataframe
'''
    messagebox.showinfo('Description', description) 
    
       

def dependencies():
    dependencies = ''' Welcome to this pseudo read_me.txt.

The whole thing is written in python 3.11.2.


Modules and librairies required : 
    - tkinter Tcl/Tk 8.6
    - PIL (Image, ImageTK)
    - pandas 1.5.3
    - itertools
    - numpy 1.24.2
    - matplotlib.pyplot 3.7.1
    

The code is divided in three python files :
    - auto_prot_app.py
        contains the Tkinter GUI, and hopefully will run the automation
    - app_def.py 
        contains funtions necessary to the "menu"'s buttons, can be copy / pasted to the main code without many problems
    - auto_prot_def.py
        contains all the functions to automate the given proteomics data treatment. There is no "__init__", as I have no idea what is it or how to do it. Functions are called in auto_prot_app to pass prints as tkinter labels. 

The logo is encoded in UTF-8 and display as such. 
'''
    messagebox.showinfo('Dependencies', dependencies)
    
    


