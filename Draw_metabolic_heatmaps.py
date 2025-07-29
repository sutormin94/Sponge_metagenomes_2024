###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2024##
##Draw heatmaps##

#Takes data for pathways presence in different genomes and draws heatmaps.
###############################################

#######
#Packages to be imported.
#######


import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import seaborn as sb
import math as math
import os
from Bio import SeqIO, Align, Seq

#################
### Variables to be defined.
#################

# Path to raw data table.
Input_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\Amino_acid_and_vitamin_biosynthetic_pw_heatmaps\\Amino_acids_vitamins_sec_systems.xlsx"

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\Amino_acid_and_vitamin_biosynthetic_pw_heatmaps\\"




def draw_heatmap_aa(input_df, output_path):
    
    fig=plt.figure(figsize=(6,9), dpi=300)
    ax=fig.add_subplot(111)
    
    ax=sb.heatmap(input_df, annot=True, fmt=".1f", annot_kws={'fontsize' : 8, 'fontweight' : 'bold'}, linewidth=0.3, 
                  cbar_kws={'shrink' : 0.3, 'aspect' : 10}, cmap='Blues', square=True, xticklabels=True, yticklabels=True, ax=ax)
    cbar=ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)    
    ax.tick_params(axis='x', labelsize=15, rotation=90)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_ylabel('Biosynthetic pathway', fontsize=20) 
    ax.xaxis.tick_top()    
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'Amino_acids_simplified_heatmap.png'), dpi=300)
    plt.savefig(os.path.join(output_path, 'Amino_acids_simplified_heatmap.svg'), dpi=300)
    
    return


def draw_heatmap_aa_transp(input_df, output_path):
    
    fig=plt.figure(figsize=(11.5,5.78), dpi=300)
    ax=fig.add_subplot(111)
    
    ax=sb.heatmap(input_df, annot=True, fmt=".1f", annot_kws={'fontsize' : 8, 'fontweight' : 'bold'}, linewidth=0.3, 
                  cbar_kws={'shrink' : 0.3, 'aspect' : 10}, cmap='Blues', square=True, xticklabels=True, yticklabels=True, ax=ax)
    cbar=ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)    
    ax.tick_params(axis='x', labelsize=15, rotation=90)
    ax.tick_params(axis='y', labelsize=13)
    ax.set_ylabel('Transporter', fontsize=20) 
    ax.xaxis.tick_top()    
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'Amino_acids_transporters_simplified_heatmap.png'), dpi=300)
    plt.savefig(os.path.join(output_path, 'Amino_acids_transporters_simplified_heatmap.svg'), dpi=300)
    
    return


def draw_heatmap_vit(input_df, output_path):
    
    fig=plt.figure(figsize=(9,8.68), dpi=300)
    ax=fig.add_subplot(111)
    
    ax=sb.heatmap(input_df, annot=True, fmt=".1f", annot_kws={'fontsize' : 8, 'fontweight' : 'bold'}, linewidth=0.3, 
                  cbar_kws={'shrink' : 0.3, 'aspect' : 10}, cmap='Blues', square=True, xticklabels=True, yticklabels=True, ax=ax)
    cbar=ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)    
    ax.tick_params(axis='x', labelsize=15, rotation=90)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_ylabel('Biosynthetic pathway', fontsize=20) 
    ax.xaxis.tick_top()    
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'Vitamins_simplified_heatmap.png'), dpi=300)
    plt.savefig(os.path.join(output_path, 'Vitamins_simplified_heatmap.svg'), dpi=300)
    
    return


def draw_heatmap_vit_transp(input_df, output_path):
    
    fig=plt.figure(figsize=(5.13,5.03), dpi=300)
    ax=fig.add_subplot(111)
    
    ax=sb.heatmap(input_df, annot=True, fmt=".1f", annot_kws={'fontsize' : 8, 'fontweight' : 'bold'}, linewidth=0.3, 
                  cbar_kws={'shrink' : 0.3, 'aspect' : 10}, cmap='Blues', square=True, xticklabels=True, yticklabels=True, ax=ax)
    cbar=ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)    
    ax.tick_params(axis='x', labelsize=15, rotation=90)
    ax.tick_params(axis='y', labelsize=13)
    ax.set_ylabel('Transporter', fontsize=20) 
    ax.xaxis.tick_top()    
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'Vitamins_transporters_simplified_heatmap.png'), dpi=300)
    plt.savefig(os.path.join(output_path, 'Vitamins_transporters_simplified_heatmap.svg'), dpi=300)
    
    return


def draw_heatmap_transp_sec(input_df, output_path):
    
    fig=plt.figure(figsize=(5.13,2.7), dpi=300)
    ax=fig.add_subplot(111)
    
    ax=sb.heatmap(input_df, annot=True, fmt=".0f", annot_kws={'fontsize' : 11, 'fontweight' : 'bold'}, linewidth=0.3, 
                  cbar_kws={'shrink' : 0.5, 'aspect' : 10}, cmap='Reds', square=False, xticklabels=True, yticklabels=True, ax=ax)
    #mask=(input_df==0),
    cbar=ax.collections[0].colorbar
    cbar.set_ticks([0,2,4])
    cbar.set_ticklabels([0,2,4])    
    cbar.ax.tick_params(labelsize=15)    
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=13, rotation=0)
    #ax.set_ylabel('', fontsize=20) 
    ax.xaxis.tick_top()    
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'Transport_and_secretion_systems_simplified_heatmap.png'), dpi=300)
    plt.savefig(os.path.join(output_path, 'Transport_and_secretion_systems_simplified_heatmap.svg'), dpi=300)    
    
    return


def wrapper_func(input_data_path, output_path):
    
    # Read input data for amino acids biosynthesis.
    Input_aa_df=pd.read_excel(input_data_path, sheet_name="KEGG_aa_symb_strains_simpl", index_col=0)
    
    # Draw heatmap for amino acids biosynthesis.
    draw_heatmap_aa(Input_aa_df, output_path)
    
    # Read input data for amino acids transporters.
    Input_aa_transp_df=pd.read_excel(input_data_path, sheet_name="KEGG_aa_symb_strains_transp", index_col=0)
    
    # Draw heatmap for amino acids transporters.
    draw_heatmap_aa_transp(Input_aa_transp_df, output_path)    
    
    # Read input data for vitamins biosynthesis.
    Input_vit_df=pd.read_excel(input_data_path, sheet_name="KEGG_vit_symb_strains_simpl", index_col=0)
    
    # Draw heatmap for vitamins biosynthesis.
    draw_heatmap_vit(Input_vit_df, output_path)    
    
    # Read input data for vitamins transporters.
    Input_vit_transp_df=pd.read_excel(input_data_path, sheet_name="KEGG_vit_symb_strains_transp", index_col=0)
    
    # Draw heatmap for vitamins transporters.
    draw_heatmap_vit_transp(Input_vit_transp_df, output_path)   
    
    # Read input data for transport and secretion systems.
    Input_transp_sec_df=pd.read_excel(input_data_path, sheet_name="Transport_and_secretion_systems", index_col=0)
    
    # Draw heatmap for transport and secretion systems.
    draw_heatmap_transp_sec(Input_transp_sec_df, output_path)      
    
    return


wrapper_func(Input_data_path, Output_path)