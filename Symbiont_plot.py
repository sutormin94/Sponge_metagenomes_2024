###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2023##
##Symbiont-plot construction##

#Takes input table with absolute abundances of microbial features (OTU or ASV) for 
#two communnities (e.g., sponge microbiome and marine water microbiome).
#Plots logs for one communnity vs logs for another, taking into account the log(relative abundance ratio).
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

#Dataset_names.
dataset_name_1 = "H. sitiens 2016 1"
dataset_name_2 = "Marine water 2016"

#Output plot name.
image_name = "H_sitiens_2016 1"

#Relative abundance threshold.
Rel_abund_thresh = 0.01

#Log ratio threshold.
Log_ratio_thresh = 1.6989

#Legend points sparsing.
sparsing=2


#Create 2 plots on one figure.
fig=plt.figure(figsize=(5, 3.5)) 
gs=gridspec.GridSpec(1, 2, width_ratios=[5, 1]) 
ax0=plt.subplot(gs[0])
ax1=plt.subplot(gs[1])

#Get data from the input table.
x = table[dataset_name_1]  #x - absolute abundance in the community 1. 
y = table[dataset_name_2]  #y - absolute abundance in the community 2.
a = list(table.index)      #OTU IDs.

#A pseudocount is added (+1) to get rid of 0 values.
x = x+1
y = y+1

#Calculate relative abundance of OTUs.
Total_counts_x = float(sum(x))
Total_counts_y = float(sum(y))
x = x/Total_counts_x
y = y/Total_counts_y

#log10 of relative abundance of OTU.
x_log = np.log10(x)
y_log = np.log10(y)

#Calculate log(ratio of relative abundances of the community 1 to the relative abundances in the community 2)
z = np.log10(x/y) 

#Calculate a radius of the OTU point. Calculated as a ((x/y)/3.14)^0.5
n = ((x/y)/3.14)**0.5

#Calculate zero values (considering the pseudocounts added).
Min_freq_log_x = np.log10(1/Total_counts_x)
Min_freq_log_y = np.log10(1/Total_counts_y)


#Prepare data for a customized legend.
z_max=math.ceil(max(z))
z_min=math.ceil(min(z))-1
print(f'Max log (ceil)={z_max}; min log (ceil)={z_min}')
if z_max-z_min<=10:
    print('Here', z_max, z_min)
    #num_points=(z_max-z_min)*2
    num_points=(z_max-z_min)
    z_legend=[]
    for i in range(z_max-z_min):
        z_legend.append(10**(z_max-i))
        #z_legend.append(10**(z_max-i-np.log10(10/3)))
    n_legend=np.sqrt(np.asarray(z_legend)/np.pi)
    y_legend=sorted(np.linspace(int(min(y))-2, int(max(y)), num_points), reverse=True)
else:
    num_points=math.ceil((z_max-z_min)/sparsing)
    z_legend=[]
    for i in range(math.ceil((z_max-z_min)/sparsing)):
        z_legend.append(10**(z_max-(sparsing*i)))  
    n_legend=np.sqrt(np.asarray(z_legend)/np.pi)
    y_legend=sorted(np.linspace(int(min(y))-2, int(max(y)), num_points), reverse=True)

x_anch=0
x_legend=np.asarray([x_anch]*len(y_legend))

z_annotate=[]
for i in range(len(z_legend)):
    if z_legend[i]<1:
        z_annotate.append(float("{:.0e}".format(z_legend[i])))
    else:
        z_annotate.append(int(float("{:.0e}".format(z_legend[i]))))
 
print(x_legend)
print(y_legend)
print(z_legend)
print(n_legend)


#Draw the legend subplot.
ax1.scatter(x_legend, y_legend, s=n_legend*100, c=np.log10(z_legend), cmap='jet', edgecolors='black', linewidth=0.3, zorder=10)
for i in range(len(y_legend)):
    ax1.annotate(z_annotate[i], (x_legend[i]+0.03, y_legend[i]), va='center', fontsize=11, zorder=20)
    ax1.annotate(z_annotate[i], (x_legend[i]+0.03, y_legend[i]), va='center', fontsize=11, zorder=20)
ax1.set_ylim([min(y_legend)-2, max(y_legend)+1])
ax1.axis('off')

#Prepare threshold y_log=x_log-THR. x should be 10 times more abundant than y if THR=1.
t=np.linspace(min(x_log)-1,max(x_log)+1,10)
tt=t-Log_ratio_thresh

#Draw data subplot.
ax0.scatter(x_log, y_log, s=n*100, c=z, cmap='jet', edgecolors='black', linewidth=0.3, zorder = 10)
ax0.axvline(np.log10(Rel_abund_thresh), linestyle='--', linewidth=1, color='black')
ax0.plot(t, tt, linestyle='--', linewidth=1, color='black')
ax0.set_xlabel(dataset_name_1, size=14)
ax0.set_ylabel(dataset_name_2, size=14)
ax0.tick_params(axis='both', labelsize=10)

#Add annotation. Select OTUs which satisfy two criteria:
#1) Relative abundance in the first community > 0.005.
#2) x should be 10 times more abundant than y (THR=1).
x_log_abundant=x_log[x_log>np.log10(Rel_abund_thresh)]
x_log_symbionts=x_log_abundant[x_log>y_log+Log_ratio_thresh]
Putative_symbionts=list(x_log_symbionts.index)
print(f'List of putative symbionts: {Putative_symbionts}')
for OTU in Putative_symbionts:
    rand=(np.random.rand(1,1)-0.5)
    print(rand)
    ax0.annotate(OTU, (x_log.loc[OTU], y_log.loc[OTU]), xytext=(x_log.loc[OTU], y_log.loc[OTU]+rand), fontsize=11, zorder=100)
ax0.axis([-5, 1, -6, -0.6])
ax0.set_xticks([Min_freq_log_x,-4, -3, -2, -1, 0])
ax0.set_yticks([Min_freq_log_y, -4, -3,-2, -1, 0])
ax0.set_xticklabels(['0', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$'], fontsize=12)
ax0.set_yticklabels(['0', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$'], fontsize=12)
ax0.set_xlim([Min_freq_log_x-0.5, 0.5])
ax0.set_ylim([Min_freq_log_y-0.5, 0.5])
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()


fig.savefig(os.path.join(path, "symbiont_plots_tables", f'{image_name}.svg'), dpi=600)
fig.savefig(os.path.join(path, "symbiont_plots_tables", f'{image_name}.png'), dpi=600)