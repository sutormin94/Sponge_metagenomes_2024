###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2024##
##Trim gff and fasta files for visualization in clinker##

#Takes input from InterProsacn annotation of MAGs and extracts predicted ELP-containing proteins.
###############################################

#######
#Packages to be imported.
#######


import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
import numpy as np
import pandas as pd
import math as math
import scipy as sp
import os
from Bio import SeqIO, Align, Seq

#################
### Variables to be defined.
#################

# Path to initial gff file.
Path_to_gff="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Metabolism_analysis\Sulfoacetate_clusters\\to_trimm\I_palmata_scaffold_9.gff3"

# Path to initial fasta file.
Path_to_fasta="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Metabolism_analysis\Sulfoacetate_clusters\\to_trimm\I_palmata_scaffold_9.fasta"

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Metabolism_analysis\Sulfoacetate_clusters\\to_trimm\\"



#######
# Read fasta.
#######

def read_fasta(path_to_fasta):
    
    with open(path_to_fasta) as fasta_f:
        fasta_records=list(SeqIO.parse(fasta_f, "fasta"))
    
    return fasta_records


#######
# Read gff3.
#######

def read_gff3(path_to_gff):
    
    gff_filein=open(path_to_gff, 'r')
    
    region_info=[]
    Genes_CDS_info=[]
    
    for line in gff_filein:
        
        if (line[0]=='#' and ('sequence-region' not in line)) or (len(line)<5):
            continue
        elif line[0]=='#' and ('sequence-region' in line):
            line=line.rstrip().split(' ')
            accession=line[1]
            region_start=int(line[2])
            region_end=int(line[3])
            
            region_info=[accession, region_start, region_end]
            
        else:
            line=line.rstrip().split('\t')
            line[3]=int(line[3])
            line[4]=int(line[4])
            Genes_CDS_info.append(line)
            
    gff_filein.close()
    
    return region_info, Genes_CDS_info


#######
# Trim gff3.
#######

def trim_gff(region_info, Genes_CDS_info):
    
    Flanks=150
    
    Region_hard_start=Genes_CDS_info[0][3]
    Region_hard_end=Genes_CDS_info[-1][4]
    
    Region_soft_start=Region_hard_start-Flanks
    if Region_soft_start>=region_info[1]:
        Region_new_start=Region_soft_start
    else:
        Region_new_start=region_info[1]
        
    Region_soft_end=Region_hard_end+Flanks
    if Region_soft_end<=region_info[2]:
        Region_new_end=Region_soft_end
    else:
        Region_new_end=region_info[2]   
        
    region_info[1]=1
    region_info[2]=Region_new_end-Region_new_start
    
    for i in range(len(Genes_CDS_info)):
        
        Genes_CDS_info[i][3]=Genes_CDS_info[i][3]-Region_new_start
        Genes_CDS_info[i][4]=Genes_CDS_info[i][4]-Region_new_start
    
    return Region_new_start, Region_new_end, region_info, Genes_CDS_info


#######
# Write trimmed gff3.
#######

def write_gff3(region_info, Genes_CDS_info, output_path):
    
    gff3_outfile=open(output_path, 'w')
    
    gff3_outfile.write('##gff-version 3\n')
    gff3_outfile.write(f'##sequence-region {region_info[0]} {region_info[1]} {region_info[2]}\n')
    
    for Gene_CDS_line in Genes_CDS_info:
        
        gff3_outfile.write(f'{Gene_CDS_line[0]}\t{Gene_CDS_line[1]}\t{Gene_CDS_line[2]}\t{Gene_CDS_line[3]}\t{Gene_CDS_line[4]}\t{Gene_CDS_line[5]}\t{Gene_CDS_line[6]}\t{Gene_CDS_line[7]}\t{Gene_CDS_line[8]}\n')
    
    gff3_outfile.close()
    
    return


#######
# Trim fasta.
#######

def trim_fasta(fasta_records, Region_new_start, Region_new_end):
    
    for record in fasta_records:
        
        record.seq=record.seq[Region_new_start-1:Region_new_end]
        
    return fasta_records


#######
# Write trimmed fasta.
#######

def write_fasta(fasta_records, output_path):
    
    fasta_outfile=open(output_path, 'w')
    
    for record in fasta_records:
        
        fasta_outfile.write(f'{record.id}\n{record.seq}\n')
        
    fasta_outfile.close()
    
    return


#######
# Wrapper function.
#######

def wrapper_func(path_to_fasta, path_to_gff, output_path):
    
    # Read gff.
    region_info, Genes_CDS_info=read_gff3(path_to_gff)
    
    # Read fasta
    fasta_records=read_fasta(path_to_fasta)
    
    # Trim gff.
    Region_new_start, Region_new_end, region_info_upd, Genes_CDS_info_upd=trim_gff(region_info, Genes_CDS_info)
    
    # Trim fasta.
    fasta_records_upd=trim_fasta(fasta_records, Region_new_start, Region_new_end)
    
    # Write gff.
    path_to_gff_suf=path_to_gff.split('.')[0]
    output_path_gff=f'{path_to_gff_suf}_trimmed.gff3'
    write_gff3(region_info_upd, Genes_CDS_info_upd, output_path_gff)
    
    # Write fasta.
    path_to_fasta_suf=path_to_fasta.split('.')[0]
    output_path_fasta=f'{path_to_fasta_suf}_trimmed.fasta'    
    write_fasta(fasta_records_upd, output_path_fasta)
    
    return

wrapper_func(Path_to_fasta, Path_to_gff, Output_path)


