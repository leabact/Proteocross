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

import pandas as pd
import itertools
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk 
from tkinter import messagebox, simpledialog, filedialog, ttk
from PIL import Image, ImageTk



######################################## Start / Main functions ########################################


def get_data(how_many):
    
    bands_dict = {} # structure : {'B1' : [dataframe, minimal MW, maximal MW] }

    for i in range(how_many):
        file_path = filedialog.askopenfilename()
        band_df = pd.read_excel(file_path,sheet_name="Protein sets")
        band_df.name='B'+str(i+1)
        name = 'B'+str(i+1)
        file_name = file_path.split('/')[-1]
        
        def get_mini_maxi():
    
            mini_max = True
            while mini_max :
                mini_maxi = simpledialog.askstring(title="MW cut-offs", prompt="You have selected "+file_name+". \n Please write the moleculars cut-off ⚠ IN DALTON ⚠  : \n (format : minimum maximum, separated by a space)")

                try :
                    mini = int(mini_maxi.split(' ')[0])
                    maxi = int(mini_maxi.split(' ')[1])

                    if mini<maxi :          
                        mini_max= False
                        return(mini, maxi)

                    else : 
                        messagebox.showinfo('Not good!','You have to write minimum then maximum, as numbers separated by a space.')
                        pass
            
                except :  
                     messagebox.showinfo('Not good!','You have to write \n minimum maximum \n as numbers separated by a space.') 
        

        mini, maxi = get_mini_maxi()  

        bands_dict[name] = [band_df, mini, maxi]

    prints = []
    prints.append('You have given the following data :')
    for k,v in bands_dict.items() :
        prints.append(k)
        prints.append((len(v[0]),'identified proteins'))
        prints.append(('Minimal molecular weight :',v[1],', Maximal molecular weight :',v[2]))


    prints.append("Please wait while the data is processed.")
    
    return(bands_dict, prints, file_path)


######################################## Bands treatments functions ########################################


def keep_myc(df, right_frame):
    # contaminants' removal
    human=0
    mouse=0
    myc = 0
    other = 0
    
    name=df.name
    nbgenes=len(df)
    
    myc_genes = []
    
    for i in df['accession'] :
        if "HUMAN" in i:
            human += 1
            df = df.drop(df[df['accession'] == i].index[0],axis=0)
        elif "MOUSE" in i :
            mouse += 1
            df = df.drop(df[df['accession'] == i].index[0],axis=0)
        elif "MYC" in i :
            myc += 1
        else : 
            other += 1
            df = df.drop(df[df['accession'] == i].index[0],axis=0)
    
    keep_print = ''
    keep_print=str("On "+str(nbgenes)+" proteins identified in "+name+", there are "+str(human+mouse+other)+" from contaminants and "+str(myc)+" from mycobacteria.")
    
    return(df, keep_print)



def cut_MW(df,MWmin,MWmax, right_frame):
    # only keep given molecular weights 
    for i in df['MW']:
        if i > MWmax or i < MWmin :
            df = df.drop(df[df['MW'] == i].index[0],axis=0)
    return(df)    


def bands_filter(bands_dict, right_frame, current_row):
    # apply keep_myc and cut_MW for each df of bands
    # return list of accessions as sets by band and list of clean dataframes
    
    # print title of what we're doing
    label = tk.Label(right_frame, text = 'Contamintants removal and molecular weight filter :', font = ('Arial',12), justify = 'left', anchor='w')
    label.grid(row = current_row, column = 0, pady=10)
    current_row += 1
    
    clean_df = []
    access = {}
    keep_prints = []

    for df,v in bands_dict.items():
        v[0], keep_print = keep_myc(v[0], right_frame)
        v[0] = cut_MW(v[0],v[1],v[2], right_frame)
        v[0].name = df
        clean_df.append(v[0])
        keep_prints.append(keep_print)

    keep = ''
    for x in keep_prints :
        keep += str(x)+'\n'
    
    label = tk.Label(right_frame, text = keep, font = ('Arial',12), justify = 'left', anchor="w")
    label.grid(row = current_row, column = 0, pady=10)
    current_row += 1
    
    txt = 'After applying the filters eliminating contaminants and non-desired molecular wight, there is :\n'
    i=1
    for df in clean_df :
        txt += str(len(df))+' mycobacterial proteins in '+df.name+'\n'
        access[i] = set(df['accession'])
        i += 1
    
    label = tk.Label(right_frame, text = txt, font = ('Arial',12), justify = 'left', anchor='w')
    label.grid(row = current_row, column = 0, pady=10)
    current_row += 1
    
    return(access,clean_df, current_row)


def common_prot_in_bands(access,clean_df, right_frame, current_row):
    # gather common proteins between two bands
    # return a dict { 'common Bx By' : df }
    common_results = {}
    
    txt = 'Protein in common between bands :'

    for k1,k2 in itertools.combinations(access.keys(),2):

        # verify if values (accession numbers) in common between two bands, tmp is a set
        tmp = access[k1].intersection(access[k2])

        # if so, print how many & save in a dataframe
        if tmp :

            txt += "\n There are "+str(len(tmp))+" proteins in common between B"+str(k1)+" and B"+str(k2)+" ."
            txt += " Results will be saved in a specified excel sheet."

            # need to access info from whole df, stored in "clean_df"
            # k1 is the current band, k2 is the next but as index number we need to acces by -1 (start at 0)
            # clean_df[0] is 1st band, clean_df[1] is 2nd band, etc

            # register common accessions in new dataframe from first then second band analyzed
            same_prot_k1 = clean_df[k1-1].loc[clean_df[k1-1]['accession'].isin(tmp)]
            same_prot_k1 = same_prot_k1.assign(band="B1")
            same_prot_k1.reset_index(inplace=True,drop=True)

            same_prot_k2 = clean_df[k2-1].loc[clean_df[k2-1]['accession'].isin(tmp)]
            same_prot_k2 = same_prot_k2.assign(band="B2")
            same_prot_k2.reset_index(inplace=True,drop=True)

             # merge the two result dataframes 
            df_same_prot = [same_prot_k1,same_prot_k2]
            same_prot = pd.concat(df_same_prot)
            # organize it by band and accession to compare identifications between the two bands
            same_prot.sort_values(by=['accession','band'], axis=0, inplace=True)
            # give it a name
            same_full_name = 'common B'+str(k1)+' B'+str(k2)

            # get specific proteins by bands
            prot_only_k1 = clean_df[k1-1].loc[~clean_df[k1-1]['accession'].isin(tmp)]
            prot_only_k1.name = 'specific B'+str(k1)
            only_k1_name = 'specific B'+str(k1)

            prot_only_k2 = clean_df[k2-1].loc[~clean_df[k2-1]['accession'].isin(tmp)]
            prot_only_k2.name = 'specific B'+str(k2)
            only_k2_name = 'specific B'+str(k2)

            # save everthing in a dictionnary 
            common_results[same_full_name] = same_prot
            common_results[only_k1_name] = prot_only_k1
            common_results[only_k2_name] = prot_only_k2

        else : 
            txt += "\n There is no protein in common between B"+str(k1)+" and B"+str(k2)+" ."
        
    label = tk.Label(right_frame, text = txt, font = ('Arial',12), justify = 'left', anchor='w')
    label.grid(row = current_row, column = 0, pady=10)
    current_row += 1

    return(common_results, current_row)


######################################## Common functions ########################################


def path_to_save(path,new_file_name):
    # get the last path given, split on the last '/' to delete the opened file name
    # add the new file name
    # return the path to the new file
    return(path.rsplit('/',1)[0]+'/'+new_file_name)



def saving_bands(saving_path, clean_df, common_results, right_frame, current_row):
        #open an excel writer
    with pd.ExcelWriter(saving_path) as writer:

        # all proteins by band (after filter)
        for df in clean_df:
            df.to_excel(writer, sheet_name=df.name, index=False)

        # common and specific proteins between bands
        for name,df in common_results.items():
            df.to_excel(writer, sheet_name=name, index = False)
            
    label = tk.Label(right_frame, text = '\n Results are saved in the same folder as your last given band.',
                     font = ('Arial',12), justify = 'left', anchor='w')
    label.grid(row = current_row, column = 0, pady=10)
    current_row += 1
    return(current_row)

    # all done



######################################## Pulldown functions ########################################


def auto_pulldown(whole_pulld, clean_df, right_frame, current_row):
    '''
    auto pull down = remove contaminants and check if common proteins in bands
    returns a dictionary of dataframes with proteins found in bands and pulldown
    '''
    # remove contaminants
    whole_pulld.name = 'pulldown'
    txt = '\n Pull-down treatment : '
    pulldown_res, to_add_txt = keep_myc(whole_pulld, right_frame)
    txt += to_add_txt
    
    # for all myc proteins, do -log10 pval and log2 ratio for future plot
    # rename the t-test column as t_test 
    pulldown_res = pulldown_res.rename(columns={'t-test_g1_vs_g2':'t_test_g1_vs_g2'})
    # log2 ratio
    np.seterr(divide = 'ignore')
    pulldown_res = pulldown_res.assign(Ratio_Log2 = lambda x : np.log2(x.ratio_g1_vs_g2 ))
    # -log10 t-test
    pulldown_res = pulldown_res.assign(T_test_Log10 = lambda x : -(np.log10(x.t_test_g1_vs_g2)))
    

    # if ratio = 0 then replace its log2 with -7.78, if ttest = 0 then replace its -log10 with 6 (lowest and highest)
    txt += '\n'
    ttest = pulldown_res.loc[pulldown_res['t_test_g1_vs_g2']==0]
    if not ttest.empty: 
        acc_ttest = ''
        for row in ttest['accession']:
            acc_ttest += str(row)+' '
        txt += '\n'+acc_ttest+'show a t-test of exactly 0. The -Log10 value will be change to 6 (highest).'

    ratio = pulldown_res.loc[pulldown_res['ratio_g1_vs_g2']==0]
    if not ratio.empty:
        acc_ratio = ''
        for row in ratio['accession']:
            acc_ratio += row+' '
        txt += '\n'+acc_ratio+'show a ratio of exactly 0. The Log2 value will be change to -7.78 (lowest).'
    
    pulldown_res.replace([np.inf,-np.inf], [6, -7.78], inplace = True)
    
    
    # apply filters : pval <= 0.05 and ratio > 2
    pulld_access = pulldown_res.loc[(pulldown_res['t_test_g1_vs_g2'] <= 0.05) & (pulldown_res['ratio_g1_vs_g2'] > 2)]
    pulld_access = pulld_access.reset_index(drop=True)
    
    txt += "\n \n There are "+str(len(pulld_access))+" proteins identified with a p-value ≤ 0.05 and a ratio > 2. \n"

    # for each given band, check if accessions in common
    pulld_bands = {}

    for df in clean_df : 
        res_df = pulld_access.loc[pulld_access['accession'].isin(df['accession'])]
        if not res_df.empty :
            name = 'pull '+df.name
            pulld_bands[name] = res_df
            txt += str((len(res_df)))+" protein.s from "+df.name+" found enriched in pulldown.\n"
        else : 
            txt += "No proteins identified in band "+df.name+" were recovered in pulldown.\n"


    label = tk.Label(right_frame, text = txt, font = ('Arial',12), justify = 'left', anchor='w')
    label.grid(row = current_row, column = 0, pady=10)
    current_row += 1
    
    return(pulld_bands,pulld_access,pulldown_res, current_row)



def volcano_plot(pulldown_res, pulld_access, pulld_bands, path_pd):
    # draw and save the volcano plot, return its path 

    fig, ax = plt.subplots()

    plt.scatter(pulldown_res['Ratio_Log2'],pulldown_res['T_test_Log10'],color='#5DADE2',s=12)
    legend = ['All proteins identified in pulldown']

    plt.scatter(pulld_access['Ratio_Log2'],pulld_access['T_test_Log10'],color='gold',s=12)
    legend.append('p-value ≤ 0.05 and ratio > 2')

    colors = ['#8E44AD','#2ECC71','#ff3b58','#d5a0bb','#b87439','#1F618D']

    i = 0
    for name,df in pulld_bands.items():
        plt.scatter(df['Ratio_Log2'],df['T_test_Log10'],color=colors[i],s=20)
        legend.append(name.split(' ')[1])
        i +=1

    plt.legend(legend,framealpha=0,bbox_to_anchor=(0.315, 0.98))

    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.xlabel('Ratio (Log2)',loc='center')
    ax.text(0,0.5,'p-value (Log10)',rotation='vertical',transform=ax.transAxes,ha='center',va='center')

    plt.grid(color='#E5E8E8',alpha=0.3,linestyle = 'dashed')

    plt.axhline(1.3,color='#E5E8E8',alpha=0.5)
    plt.text(5.3,1.32,'p-val. ≤ 0.05',color='#566573',size=9)
    plt.axvline(1,color='#E5E8E8',alpha=0.5)
    plt.text(1.1,0,'ratio > 2',color='#566573',size=9)

    plt.title('Volcano plot of identified proteins in pull down, cross resulted with identified proteins in bands\n')

    volcano_saving_path = path_to_save(path_pd,'Volcano plot.png')
    plt.savefig(volcano_saving_path,transparent=True,format='png',bbox_inches='tight')

    plt.close(fig)
  
    return(volcano_saving_path)


def get_df_data_to_display(pulld_bands, pulld_access):
    
    def get_prot_data_to_display(df):
        # extract nice info to show
        cols_to_keep = ['accession', 'gene_name', 'description', 'protein_set_score', 'coverage',
                    'MW','ratio_g1_vs_g2', 't_test_g1_vs_g2', 'Ratio_Log2', 'T_test_Log10']
        df = df[cols_to_keep]
        df.reset_index(drop=True,inplace=True)
        #round numbers to display
        new_df = df.round({'protein_set_score':2,'Ratio_Log2':2, 'T_test_Log10':2})

        # cut description 
        new_df.insert(2,'Description', df['description'].str.split('OS=').str[0])
        new_df.drop('description', axis=1,inplace=True)

        # first letter in capitale
        new_df = new_df.rename(columns=lambda x: x.capitalize())
        new_df = new_df.rename(columns={'Mw':'MW'})

        # save the rounded df to display
        return(new_df)

    
    # from pulld_bands
    pulld_bands_info = {}
    for k,df in pulld_bands.items():
        new_df = get_prot_data_to_display(df)
        pulld_bands_info[k] = new_df

    # from pulld_access then create pulld_top15
    pulld_top15 = pd.DataFrame()
    pulld_top15 = get_prot_data_to_display(pulld_access)
    pulld_top15 = pulld_top15.sort_values(by = 'Ratio_g1_vs_g2', ascending = False).head(15)
    pulld_top15.reset_index(inplace = True, drop = True)

    #return the two df
    return(pulld_bands_info, pulld_top15)


def make_df_to_tree(df, name, current_row, color, right_frame):
    
        frame = tk.Frame(right_frame, width=600, bg=color)
        frame.grid(row = current_row, column = 0) 

        disp_name = tk.Label(frame, text=name)
        disp_name.pack()

        df_list = list(df.columns.values)
        df_rset = df.to_numpy().tolist()

        # Calculate the height of the TreeView based on the number of rows in the DataFrame
        tree_height = len(df_rset) + 1  # Add 1 to account for the header row

        df_tree = ttk.Treeview(frame, columns=df_list, height=tree_height)
        df_tree['show'] = 'headings'

        for i in df_list:
            df_tree.column(i, width=85, anchor='c')
            df_tree.heading(i, text=i.replace('_', ' '))
            if i == 'Accession':
                df_tree.column(i, width=100, anchor='c')
            elif i == 'Description':
                df_tree.column(i, width=400, anchor='c')
            elif i == 'Coverage':
                df_tree.column(i, width=60, anchor='c')
            elif i == 'MW':
                df_tree.column(i, width=55, anchor='c')

        for dt in df_rset:
            v = [r for r in dt]
            df_tree.insert('', 'end', values=v)

        df_tree.pack()


def display_info_bands(pulld_bands_info, right_frame, current_row):
    # use pulld_bands_info to display a treeview (board) of info by bands
    colors = ['#8E44AD','#2ECC71','#ff3b58','#d5a0bb','#b87439','#1F618D']

    i=0
    for name,df in pulld_bands_info.items():
        name = name.split(' ')[1]
        make_df_to_tree(df,name, current_row, colors[i], right_frame)
        i += 1
        current_row += 1

    return(current_row)


        
def top15_volcano(pulldown_res,pulld_top15,pulld_bands, path_pd):

    #plot and save the volcano
    
    fig, ax = plt.subplots()

    plt.scatter(pulldown_res['Ratio_Log2'],pulldown_res['T_test_Log10'],
                edgecolor='#5DADE2',s=10, facecolor = 'none')
    legend = ['All proteins identified in pulldown']

    plt.scatter(pulld_top15['Ratio_log2'],pulld_top15['T_test_log10'],color='gold',s=18)
    legend.append('15 best ratio with p-value ≤ 0.05')

    colors = ['#8E44AD','#2ECC71','#ff3b58','#d5a0bb','#b87439','#1F618D']

    i = 0
    for name,df in pulld_bands.items():
        plt.scatter(df['Ratio_Log2'],df['T_test_Log10'],color=colors[i],s=20)
        legend.append(name.split(' ')[1])
        i +=1

    plt.legend(legend,framealpha=0,bbox_to_anchor=(0.315, 0.98))

    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.xlabel('Ratio (Log2)',loc='center')
    ax.text(0,0.5,'p-value (Log10)',rotation='vertical',transform=ax.transAxes,ha='center',va='center')

    plt.grid(color='#E5E8E8',alpha=0.3,linestyle = 'dashed')

    plt.title('Volcano plot of identified proteins in pull down, cross resulted with identified proteins in bands\n 15 Best ratio with  p-value ≤ 0.05')

    top15_volcano_saving_path = path_to_save(path_pd,'Volcano plot - 15 best ratio.png')
    plt.savefig(top15_volcano_saving_path, transparent=True, format='png', bbox_inches='tight')
    
    plt.close(fig)
    
    return(top15_volcano_saving_path)
    


######################################## "Main" pulldown :) ########################################

def pulldown_treatment(clean_df, right_frame, current_row):
    
    # open and read pulldown file
    path_pd = filedialog.askopenfilename()
    whole_pulld = pd.read_excel(path_pd, sheet_name="Protein sets",usecols="A:E,J:S,AZ,BA")
    
    # remove contam, filter on ratio > 2 and p-value <= 0.05, check if common in bands 
    pulld_bands, pulld_access, pulldown_res, current_row = auto_pulldown(whole_pulld, clean_df, right_frame, current_row)

    # save the results
    pd_saving_path = path_to_save(path_pd,'Pulldown-bands cross results.xlsx')
    
    with pd.ExcelWriter(pd_saving_path) as writer:
        for b,df in pulld_bands.items():
            df.to_excel(writer, sheet_name=b, index=False)
    
    
    # create and display the volcano plot of common proteins between pulld and bands
    volcano_saving_path = volcano_plot(pulldown_res, pulld_access, pulld_bands, path_pd)
    vol = Image.open(volcano_saving_path)
    vol_img = ImageTk.PhotoImage(vol)
    volcan_display = tk.Label(right_frame, image = vol_img)
    volcan_display.image = vol_img
    volcan_display.grid(row = current_row, column = 0, pady = 10)
    current_row += 1
    
    # get the info corresponding to the common proteins
    # + create a df with the top 15 protein ratio-speaking 
    pulld_bands_info, pulld_top15 = get_df_data_to_display(pulld_bands, pulld_access)
    
    #save pulld_top15 as an excel file
    top15_saving_path = path_to_save(path_pd,'Top 15 Ratio.xlsx') 
    with pd.ExcelWriter(top15_saving_path) as writer:
        pulld_top15.to_excel(writer,index=False)
    
    # display the proteins in common's info
    current_row = display_info_bands(pulld_bands_info, right_frame, current_row) # row = 10 on right frame
    
    # volcano plot the pulld_top15 + display corresponding info
    top15_path = top15_volcano(pulldown_res, pulld_top15, pulld_bands, path_pd)
    top15 = Image.open(top15_path)
    top15_img = ImageTk.PhotoImage(top15)
    top15_display = tk.Label(right_frame, image = top15_img)
    top15_display.image = top15_img
    top15_display.grid(row=current_row, column=0, pady=10)
    current_row += 1
    
    make_df_to_tree(pulld_top15, 'Top 15', current_row, 'gold', right_frame)
    current_row += 1
    
    txt = '\n Results and volcano plot are saved in the same folder as your pulldown results. \n Thanks for using this automated program!'
    label = tk.Label(right_frame, text = txt, font = ('Arial',12), justify = 'left', anchor='w')
    label.grid(row = current_row, column = 0, padx=10, pady=10)


    

######################################## "Main" band :) ########################################



def bands(bands_dict, bands_path, right_frame, current_row):
    
     # contaminants' removal + MW cutoff
    access, clean_df, current_row = bands_filter(bands_dict, right_frame, current_row)
    
    # test if proteins in common between two bands by accession number
    common_results, current_row = common_prot_in_bands(access,clean_df, right_frame, current_row)

    # save band results in an excel file
    saving_path = path_to_save(bands_path,'Bands analysis results.xlsx')
    current_row = saving_bands(saving_path, clean_df, common_results, right_frame, current_row)

    return(clean_df, current_row)
        
        