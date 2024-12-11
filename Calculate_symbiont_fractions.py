###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2024##
##Calculate fractions of symbiotic OTUs in microbiomes##

#Takes input table with absolute abundances of microbial features (OTU or ASV).
#Calculates mean relative abundance and standard deviation 
#for selected OTUs in a set of selected samples.
###############################################

#######
#Packages to be imported.
#######


import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import math as math
import os

#################
### Variables to be defined.
#################

#Path to the working directory.
path = "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\\16S_metagenomics\Symbiont_plots\\"

#Path to a table with relative abundance data.
table = pd.read_excel(os.path.join(path, 'OTU_taxa.xlsx'), sheet_name='all', header=0, index_col=0)
print(table)

H_panicea_dict={'Samples' : ['H. panicea 2016 1', 'H. panicea 2016 2', 
                             'H. panicea 2018 1', 'H. panicea 2018 2', 'H. panicea 2018 3',
                             'H. panicea 2022'], 
                'OTUs' : ['OTU 4', 'OTU 23']}
M_water_dict={'Samples' : ['Marine water 2016', 
                           'Marine water 2018', 
                           'Marine water 2022'], 
                'OTUs' : ['OTU 4', 'OTU 23', 'OTU 3', 'OTU 7', 'OTU 9', 'OTU 14', 'OTU 1']}
H_sitiens_dict={'Samples' : ['H. sitiens 2016 1', 'H. sitiens 2016 2', 'H. sitiens 2016 3',
                           'H. sitiens 2018 1', 'H. sitiens 2018 2', 'H. sitiens 2018 3', 
                           'H. sitiens 2022'], 
                'OTUs' : ['OTU 3', 'OTU 7', 'OTU 9', 'OTU 14']}
I_palmata_dict={'Samples' : ['I. palmata 2016 1', 'I. palmata 2016 2', 'I. palmata 2016 3',
                           'I. palmata 2018 1', 'I. palmata 2018 2', 
                           'I. palmata 2022'], 
                'OTUs' : ['OTU 1']}

def get_mean_rel_freq_and_std(data_dict, sponge_name):

    print(f'Working with {sponge_name} data.')
    Sponge_data=table[data_dict['Samples']]
    Total_abs_freq=Sponge_data.sum()
    
    for OTU in data_dict['OTUs']:
    
        OTU_abs_freq=Sponge_data.loc[OTU]
        OTU_rel_freq=OTU_abs_freq/Total_abs_freq
        print(f'Relative frequences for {OTU} : {OTU_rel_freq}')
        Mean_OTU_rel_freq=np.mean(OTU_rel_freq)
        STD_OTU_rel_freq=np.std(OTU_rel_freq)
        print(f'Mean realtive frequency for {OTU} : {Mean_OTU_rel_freq} with STD {STD_OTU_rel_freq}')
    
    return

#H. panicea.
get_mean_rel_freq_and_std(H_panicea_dict, 'H. panicea')
#H. sitiens.
get_mean_rel_freq_and_std(H_sitiens_dict, 'H. sitiens')
#I. palmata.
get_mean_rel_freq_and_std(I_palmata_dict, 'I. palmata')
#Marine water.
get_mean_rel_freq_and_std(M_water_dict, 'Marine water')
