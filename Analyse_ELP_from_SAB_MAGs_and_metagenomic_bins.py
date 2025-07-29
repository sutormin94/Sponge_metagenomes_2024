###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2025##
##Analyse abundance of ELP-domains and ELP-encoding genes from bacterial MAGs:
##compares SAB MAGs and other obtained metagenomic bins.

#Takes input from InterProScan annotation of MAGs and extracts predicted ELP-containing proteins.
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
from pca import pca
import colourmap

#################
### Variables to be defined.
#################

# Path to ELP definition file.
ELP_def_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\ELPs_in_SAB_MAGs_and_other_bins\ELPs_list.tsv"

# Path to folder with InterProScan output.
Interproscan_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\ELPs_in_SAB_MAGs_and_other_bins\Pfam_annotation_results\\"

# Path to folder with MGM2 output.
MGM2_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\ELPs_in_SAB_MAGs_and_other_bins\Mgm2_protein_annotation\\"

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\ELPs_in_SAB_MAGs_and_other_bins\ELP_statistics_and_figures\\"




#######
## Read ELP definitions file.
#######

def read_elp_def(elp_def_path):
    
    ELP_types_dict={}
    
    filein=open(elp_def_path, 'r')
    
    for line in filein:
        line=line.rstrip().split('\t')
        ELP_type=line[0]
        PFAM_ID=line[1]
        EKP_class=line[2]
        
        ELP_types_dict[PFAM_ID]=[ELP_type, EKP_class]
          
    filein.close()
    
    return ELP_types_dict


#######
## Read Interproscan annotation.
#######

def read_interproscan_annot(interproscan_path, ELP_types_dict):
    
    All_genes_dict={}
    
    Bact_MAGs_path_list=os.listdir(interproscan_path)
    
    for bac_file in Bact_MAGs_path_list:
        Bac_MAG_full_path=os.path.join(interproscan_path, bac_file)
        
        Bac_MAG_ar=bac_file.split('_')
        Bac_MAG_name=f'{Bac_MAG_ar[0]}_{Bac_MAG_ar[1]}_{Bac_MAG_ar[2]}_{Bac_MAG_ar[3]}'
        
        print(f'Now processing {Bac_MAG_name}')
        
        Bac_MAG_dict={}
        
        filein=open(Bac_MAG_full_path, 'r')
        for line in filein:
            line=line.rstrip().split('\t')
            gene_id=line[0]
            Pfam_ID=line[4]
            Domain_description=line[5]
            E_value=float(line[8])
            
            if gene_id in Bac_MAG_dict:
                Bac_MAG_dict[gene_id].append([Pfam_ID, Domain_description, E_value])
            else:
                Bac_MAG_dict[gene_id]=[[Pfam_ID, Domain_description, E_value]]
                
        All_genes_dict[Bac_MAG_name]=Bac_MAG_dict
        
        filein.close()
     
     
    Bac_all_MAGs_ELP_dict={}
    
    for bact_genome_name, bact_genome_dict in All_genes_dict.items():
        
        Bac_MAG_ELP_dict={}
        
        for gene_id, gene_domains_ar in bact_genome_dict.items():
            
            ELP_check=0
            
            Gene_full_info=[]
            
            for domain_info in gene_domains_ar:
                
                PFAM_ID=domain_info[0]
                Domain_desc=domain_info[1]
                E_value=domain_info[2]
                
                if PFAM_ID in ELP_types_dict:
                    
                    if E_value<=0.01:
                        
                        ELP_check=1
                        
                        ELP_type=ELP_types_dict[PFAM_ID][0]
                        ELP_class=ELP_types_dict[PFAM_ID][1]
                        
                        Domain_full_info=[PFAM_ID, Domain_desc, E_value, ELP_type, ELP_class]
                        
                    else:
                        
                        Domain_full_info=[PFAM_ID, Domain_desc, E_value, '-', '-']
                        
                else:
                    
                    Domain_full_info=[PFAM_ID, Domain_desc, E_value, '-', '-']
                    
                Gene_full_info.append(Domain_full_info)
                                          
            if ELP_check==1:
                
                Bac_MAG_ELP_dict[gene_id]=Gene_full_info
                
        Bac_all_MAGs_ELP_dict[bact_genome_name]=Bac_MAG_ELP_dict
                        
    return Bac_all_MAGs_ELP_dict


#######
## Read proteins annotation.
#######

def read_magm_annot(mgm2_path, Bac_all_MAGs_ELP_dict, output_path):
    
    Bact_mgm2_path_list=os.listdir(mgm2_path)
    
    Bac_all_MAGs_ELP_seq_dict={}
    
    for bac_file in Bact_mgm2_path_list:
        Bac_mgm2_full_path=os.path.join(mgm2_path, bac_file)  
        
        Bac_MAG_ar=bac_file.split('_')
        Bac_MAG_name=f'{Bac_MAG_ar[0]}_{Bac_MAG_ar[1]}_{Bac_MAG_ar[2]}_{Bac_MAG_ar[3]}' 
        
        print(f'Now processing {Bac_MAG_name}')
        
        Bac_all_MAG_ELP_seq_dict={}
           
        for handle in SeqIO.parse(Bac_mgm2_full_path, 'fasta'):
            
            gene_name=str(handle.name)
            gene_seq=str(handle.seq)  
            
            if gene_name in Bac_all_MAGs_ELP_dict[Bac_MAG_name]:
                
                Bac_all_MAG_ELP_seq_dict[gene_name]=gene_seq
                Bac_all_MAGs_ELP_dict[Bac_MAG_name][gene_name].append(gene_seq)
                
        Bac_all_MAGs_ELP_seq_dict[Bac_MAG_name]=Bac_all_MAG_ELP_seq_dict
            
    return Bac_all_MAGs_ELP_seq_dict, Bac_all_MAGs_ELP_dict


#######
## Write ELP annotation.
#######

def write_ELP_annotation(Bac_all_MAGs_ELP_dict, output_path):
    
    All_data_outpath=os.path.join(output_path, 'All_ELPs_found_annotation.tsv')
    All_data_fileout=open(All_data_outpath, 'w')
    
    for bac_MAG_name, bac_MAG_genes_dict in Bac_all_MAGs_ELP_dict.items():
        
        Bac_MAG_outpath=os.path.join(output_path, f'{bac_MAG_name}_ELPs_found_annotation.tsv')
        Bac_MAG_fileout=open(Bac_MAG_outpath, 'w')        
        
        for gene_id, gene_info_ar in bac_MAG_genes_dict.items():
            
            gene_seq=gene_info_ar[-1]
            
            for i in range(len(gene_info_ar)-1):
                
                domain_info=gene_info_ar[i]
                
                All_data_fileout.write(f'{bac_MAG_name}\t{gene_id}\t{domain_info[0]}\t{domain_info[1]}\t{domain_info[2]}\t{domain_info[3]}\t{domain_info[4]}\t{gene_seq}\n')
                Bac_MAG_fileout.write(f'{bac_MAG_name}\t{gene_id}\t{domain_info[0]}\t{domain_info[1]}\t{domain_info[2]}\t{domain_info[3]}\t{domain_info[4]}\t{gene_seq}\n')
                
        Bac_MAG_fileout.close()
    
    All_data_fileout.close()
    
    return


#######
## Write ELP sequences.
#######

def write_ELP_sequences(Bac_all_MAGs_ELP_seq_dict, output_path):
    
    All_data_outpath=os.path.join(output_path, 'All_ELPs_found_seq.fasta')
    All_data_fileout=open(All_data_outpath, 'w')    
    
    for bac_MAG_name, bac_MAG_ELP_seq_dict in Bac_all_MAGs_ELP_seq_dict.items():
        
        Bac_ELP_seq_outpath=os.path.join(output_path, f'{bac_MAG_name}_ELPs_found_seq.fasta')
        Bac_ELP_fileout=open(Bac_ELP_seq_outpath, 'w')        
        
        for gene_id, gene_seq in bac_MAG_ELP_seq_dict.items():
            
            Bac_ELP_fileout.write(f'>{gene_id}\n{gene_seq}\n')
            
            All_data_fileout.write(f'>{bac_MAG_name}__{gene_id}\n{gene_seq}\n')
            
    All_data_fileout.close()
    
    return
    

def wrapper_function_aggregation(elp_def_path, interproscan_path, mgm2_path, output_path):
    
    # Read ELP definitions file.
    ELP_types_dict=read_elp_def(elp_def_path)
    
    # Read Interproscan annotation.
    Bac_all_MAGs_ELP_dict=read_interproscan_annot(interproscan_path, ELP_types_dict)
    
    # Read proteins annotation.
    Bac_all_MAGs_ELP_seq_dict, Bac_all_MAGs_ELP_dict=read_magm_annot(mgm2_path, Bac_all_MAGs_ELP_dict, output_path)
    
    # Write ELP annotation.
    write_ELP_annotation(Bac_all_MAGs_ELP_dict, output_path)
    
    # Write ELP sequences.
    write_ELP_sequences(Bac_all_MAGs_ELP_seq_dict, output_path)
    
    return

wrapper_function_aggregation(ELP_def_path, Interproscan_path, MGM2_path, Output_path)



##########
########## Run subsequent analysis of found ELPs and their signal sequences.
##########

# Path to table with ELP annotation.
ELPs_in_MAGs_annotation_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\\Source_data\ELPs_in_SAB_MAGs_and_other_bins\All_ELPs_found_annotation.tsv"

# Path to SignalP output.
SignalP_results_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\\Source_data\ELPs_in_SAB_MAGs_and_other_bins\All_ELPs_found_seq_SignalP6.txt"

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\ELPs_in_SAB_MAGs_and_other_bins\ELP_statistics_and_figures\\"



#######
## Read ELPs annotation.
#######

def read_elp_annot(elps_in_mags_annotation_path):
    
    MAGs_ELPs_domains_dict={}
    
    filein=open(elps_in_mags_annotation_path, 'r')
    
    for line in filein:
        line=line.rstrip().split('\t')
        Bac_MAG_name=line[0]
        Gene_name=line[1]
        E_value=float(line[4])
        ELP_class=line[5]
        ELP_type=line[6]
        
        if Bac_MAG_name not in MAGs_ELPs_domains_dict:
            
            MAGs_ELPs_domains_dict[Bac_MAG_name]={Gene_name : [[E_value, ELP_class, ELP_type]]}
            
        else:
            
            if Gene_name not in MAGs_ELPs_domains_dict[Bac_MAG_name]:
                
                MAGs_ELPs_domains_dict[Bac_MAG_name][Gene_name]=[[E_value, ELP_class, ELP_type]]
                
            else:
                
                MAGs_ELPs_domains_dict[Bac_MAG_name][Gene_name].append([E_value, ELP_class, ELP_type])
                
    filein.close()
    
    return MAGs_ELPs_domains_dict


#######
## Read SignalP results path.
#######

def read_signalP_res(signalP_results_path):
    
    MAGs_ELPs_genes_signalP_dict={}
    
    filein=open(signalP_results_path, 'r')
    
    for line in filein:
        
        if line[0]!='#':
            
            line=line.rstrip().split('\t')
            
            MAG_gene_name=line[0].replace('___', '_+_')
            MAG_gene_name=MAG_gene_name.replace('gene_type_', 'gene_type=')
            MAG_gene_name=MAG_gene_name.replace('partial_', 'partial=')
            MAG_gene_name=MAG_gene_name.replace('_1_scaffold', '_1+scaffold')
            MAG_gene_name=MAG_gene_name.replace('rc_scaffold', 'rc+scaffold')
            MAG_gene_name=MAG_gene_name.replace(',1_scaffold', ',1+scaffold')
            
            MAG_gene_name_ar=MAG_gene_name.split('__')
            
            Bac_MAG_name=MAG_gene_name_ar[0]
            Gene_name=MAG_gene_name_ar[1]
            
            Signal_type=line[1]
            
            if Bac_MAG_name not in MAGs_ELPs_genes_signalP_dict:
                
                MAGs_ELPs_genes_signalP_dict[Bac_MAG_name]={Gene_name : Signal_type}
                
            else:
                
                if Gene_name not in MAGs_ELPs_genes_signalP_dict[Bac_MAG_name]:
                    
                    MAGs_ELPs_genes_signalP_dict[Bac_MAG_name][Gene_name]=Signal_type
                    
                else:
                    
                    print(f'Mixed signal for {Bac_MAG_name} {Gene_name} !')
                    MAGs_ELPs_genes_signalP_dict[Bac_MAG_name][Gene_name]=Signal_type
             
    filein.close()
    
    return MAGs_ELPs_genes_signalP_dict

#######
## Add signalP data to ELPs dict.
#######

def add_signalP_data(MAGs_ELPs_domains_dict, MAGs_ELPs_genes_signalP_dict):
    
    MAGs_ELPs_domains_signal_only_dict={}
    
    for Bac_MAG_name, Genes_info_dict in MAGs_ELPs_domains_dict.items():
        
        MAGs_ELPs_domains_signal_only_dict[Bac_MAG_name]={}
        
        for Gene_name, gene_domains_ar in Genes_info_dict.items():          
            
            MAGs_ELPs_domains_signal_only_dict[Bac_MAG_name][Gene_name]=[]
            
            Signal_type=MAGs_ELPs_genes_signalP_dict[Bac_MAG_name][Gene_name]
            
            for i in range(len(gene_domains_ar)):
                
                MAGs_ELPs_domains_dict[Bac_MAG_name][Gene_name][i].append(Signal_type)
                
            if Signal_type!='OTHER':
                
                MAGs_ELPs_domains_signal_only_dict[Bac_MAG_name][Gene_name]=gene_domains_ar
    
    return MAGs_ELPs_domains_dict, MAGs_ELPs_domains_signal_only_dict


#######
## Write stat for ELP types and classes for genes and domains. For all ELPs.
#######

def write_ELPs(MAGs_ELPs_domains_signalP_dict, dataset_name, output_path):
    
    # List of bacterial genomes.
    Genomes_names_ar=['I_palmata_OTU1_sp', 'Persebacter_sydneyensis_OTU1_sp', 'Beroebacter_blanensis_OTU1_sp', 'GCA_020027345_OTU1_fl', 'H_sitiens_OTU3_sp', 'GCA_007570945_OTU3_sp', 'GCA_009842505_OTU3_sp', 'GCA_002897635_OTU3_fl', 'GCA_024638545_OTU3_fl', 'H_sitiens_OTU7_sp', 'GCA_018607685_OTU7_fl', 'GCA_902567095_OTU7_fl', 'H_sitiens_OTU9_sp', 'GCA_028286085_OTU9_sp', 'GCA_014238615_OTU9_sp', 'GCA_002454015_OTU9_sp', 'H_sitiens_OTU14_sp', 'GCA_007570905_OTU14_sp', 'GCA_014323925_OTU14_sp', 'GCA_002631715_OTU14_sp', 'GCA_028286055_OTU14_sp', 'GCA_007571105_OTU14_sp', 'H_panicea_OTU4_sp', 'Halichondribacter_symbioticus_OTU4_sp', 'GCA_002742285_OTU4_fl', 'GCA_913061275_OTU4_fl', 'H_panicea_OTU23_sp', 'GCA_021296375_OTU23_sp', 'GCA_913056245_OTU23_fl']
    
    # For ELP types.
    ELP_types_ar=['TPR', 'Ank', 'CAD', 'CUB', 'NHL', 'Big', 'Fn3', 'WD40', 'PQQ', 'Sel1', 'LRR']
    
    ELP_types_stat_genes={}
    ELP_types_stat_domains={}
    
    for Bac_MAG_name, Genes_info_dict in MAGs_ELPs_domains_signalP_dict.items():
        
        ELP_types_stat_genes[Bac_MAG_name]=[0]*len(ELP_types_ar)
        ELP_types_stat_domains[Bac_MAG_name]=[0]*len(ELP_types_ar)
        
        for Gene_name, Gene_info_ar in  Genes_info_dict.items():
            
            domains_types_dict={}
            
            for domain_info in Gene_info_ar:
                
                ELP_evalue=domain_info[0]
                ELP_type=domain_info[2]
                
                if (ELP_evalue<=0.01) & (ELP_type in ELP_types_ar):
                    
                    ELP_domain_index=ELP_types_ar.index(ELP_type)    
                    ELP_types_stat_domains[Bac_MAG_name][ELP_domain_index]+=1                    
                    
                    if ELP_type not in domains_types_dict:
                
                        domains_types_dict[ELP_type]=1
                        
                    else:
                        
                        domains_types_dict[ELP_type]+=1
            
            domain_freq_init=0
            domain_type_selected=''
            
            for domain_type, domain_freq in domains_types_dict.items():
                
                if domain_freq>domain_freq_init:
                    
                    domain_freq_init=domain_freq
                    domain_type_selected=domain_type
            
            if domain_type_selected!='':   
                
                ELP_domain_index=ELP_types_ar.index(domain_type_selected)    
                ELP_types_stat_genes[Bac_MAG_name][ELP_domain_index]+=1
            
    ELP_types_stat_genes_df=pd.DataFrame.from_dict(ELP_types_stat_genes, orient='index', columns=ELP_types_ar)
    ELP_types_stat_genes_df=ELP_types_stat_genes_df.reindex(Genomes_names_ar)
    ELP_types_stat_genes_df.to_excel(os.path.join(output_path, f'{dataset_name}_ELPs_types_genes_counts.xlsx'))
    
    ELP_types_stat_domains_df=pd.DataFrame.from_dict(ELP_types_stat_domains, orient='index', columns=ELP_types_ar)
    ELP_types_stat_domains_df=ELP_types_stat_domains_df.reindex(Genomes_names_ar)
    ELP_types_stat_domains_df.to_excel(os.path.join(output_path, f'{dataset_name}_ELPs_types_domains_counts.xlsx'))    
    
    
    # For ELP classes.
    ELP_classes_ar=['TPR_1', 'TPR_10', 'TPR_11', 'TPR_12', 'TPR_14', 'TPR_15', 'TPR_16', 'TPR_17', 'TPR_18', 'TPR_19', 'TPR_2', 'TPR_20', 'TPR_21', 'TPR_3', 'TPR_4', 'TPR_6', 'TPR_7', 'TPR_8', 'TPR_9', 'Ank', 'Ank_2', 'Ank_3', 'Ank_4', 'Ank_5', 'Cadherin', 'Cadherin_3', 'Cadherin_4', 'Cadherin_5', 'Cadherin-like', 'CUB', 'NHL', 'Big_3_4', 'Big_3_2', 'Big_3_3', 'Big_1', 'Big_2', 'Big_3', 'Big_5', 'Fn3', 'Fn3_like', 'WD40', 'PQQ', 'PQQ_2', 'PQQ_3', 'Sel1', 'LRR_1', 'LRRNT_2']
    
    ELP_classes_stat_genes={}
    ELP_classes_stat_domains={}
    
    for Bac_MAG_name, Genes_info_dict in MAGs_ELPs_domains_signalP_dict.items():
        
        ELP_classes_stat_genes[Bac_MAG_name]=[0]*len(ELP_classes_ar)
        ELP_classes_stat_domains[Bac_MAG_name]=[0]*len(ELP_classes_ar)
        
        for Gene_name, Gene_info_ar in  Genes_info_dict.items():
            
            domains_classes_dict={}
            
            for domain_info in Gene_info_ar:
                
                ELP_evalue=domain_info[0]
                ELP_class=domain_info[1]
                
                if (ELP_evalue<=0.01) & (ELP_class in ELP_classes_ar):
                    
                    ELP_domain_index=ELP_classes_ar.index(ELP_class)    
                    ELP_classes_stat_domains[Bac_MAG_name][ELP_domain_index]+=1  
                    
                    if ELP_class not in domains_classes_dict:
                
                        domains_classes_dict[ELP_class]=1
                        
                    else:
                        
                        domains_classes_dict[ELP_class]+=1
            
            domain_freq_init=0
            domain_class_selected=''
            
            for domain_class, domain_freq in domains_classes_dict.items():
                
                if domain_freq>domain_freq_init:
                    
                    domain_freq_init=domain_freq
                    domain_class_selected=domain_class
            
            if domain_class_selected!='':   
                
                ELP_domain_index=ELP_classes_ar.index(domain_class_selected)    
                ELP_classes_stat_genes[Bac_MAG_name][ELP_domain_index]+=1
                    
                    
    ELP_classes_stat_genes_df=pd.DataFrame.from_dict(ELP_classes_stat_genes, orient='index', columns=ELP_classes_ar)
    ELP_classes_stat_genes_df=ELP_classes_stat_genes_df.reindex(Genomes_names_ar)
    ELP_classes_stat_genes_df.to_excel(os.path.join(output_path, f'{dataset_name}_ELPs_classes_genes_counts.xlsx'))
    
    ELP_classes_stat_domains_df=pd.DataFrame.from_dict(ELP_classes_stat_domains, orient='index', columns=ELP_classes_ar)
    ELP_classes_stat_domains_df=ELP_classes_stat_domains_df.reindex(Genomes_names_ar)
    ELP_classes_stat_domains_df.to_excel(os.path.join(output_path, f'{dataset_name}_ELPs_classes_domains_counts.xlsx'))    
    
    
    return



def wrapper_function_signalP(elps_in_mags_annotation_path, signalP_results_path, output_path):
    
    # Read ELPs annotation.
    MAGs_ELPs_domains_dict=read_elp_annot(elps_in_mags_annotation_path)
    
    # Read SignalP results path.
    MAGs_ELPs_genes_signalP_dict=read_signalP_res(signalP_results_path)
    
    # Add signalP data to ELPs dict.
    MAGs_ELPs_domains_signalP_dict, MAGs_ELPs_domains_signal_only_dict=add_signalP_data(MAGs_ELPs_domains_dict, MAGs_ELPs_genes_signalP_dict)
    
    # Write stat for ELP types and classes for genes. For all ELPs.
    write_ELPs(MAGs_ELPs_domains_signalP_dict, 'All', output_path)
    
    # Write stat for ELP types and classes for genes. For signal-containing ELPs.
    write_ELPs(MAGs_ELPs_domains_signal_only_dict, 'Exported', output_path)
    
    
    return

wrapper_function_signalP(ELPs_in_MAGs_annotation_path, SignalP_results_path, Output_path)




##########
########## Compare ELP frequences between free-living and sponge-associated MAGs.
##########


# Path to tables with ELP frequences for free-living MAGs.
ELPs_in_FL_MAGs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\\ELPs_in_SAB_MAGs_and_other_bins\\ELP_data_all_FL_genomes\\"

# Path to tables with ELP frequences for sponge-associated MAGs.
ELPs_in_SA_MAGs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\\ELPs_in_SAB_MAGs_and_other_bins\\ELP_data_all_SA_genomes\\"

# Path to CheckM info for bins.
CheckM_bins_info="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\\ELPs_in_SAB_MAGs_and_other_bins\\MAGs_and_bins_CheckM.xlsx"

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\\ELPs_in_SAB_MAGs_and_other_bins\\ELP_statistics_and_figures\\"



def norm_ELP_counts(FL_ELPs_df, SP_MAGs_ELPs_df, CheckM_info_df):
    
    Genome_size_column=CheckM_info_df['Genome_Size']/1000000
        
    FL_ELPs_df=FL_ELPs_df.div(Genome_size_column, axis='index')
    FL_ELPs_df=FL_ELPs_df.dropna(axis='index')
    
    SP_MAGs_ELPs_df=SP_MAGs_ELPs_df.div(Genome_size_column, axis='index')
    SP_MAGs_ELPs_df=SP_MAGs_ELPs_df.dropna(axis='index')    
    
    return FL_ELPs_df, SP_MAGs_ELPs_df


def draw_heatmap_all(ELP_cat, input_df, norm, output_path):
    
    # Draw heatmap with rounded to integer values.
    fig=plt.figure(figsize=(7,16), dpi=300) # 13*18 for All classes; 7*16 for types;
    ax=fig.add_subplot(111)
    
    ax=sb.heatmap(input_df, annot=True, fmt=".0f", annot_kws={'fontsize' : 7, 'fontweight' : 'normal'}, linewidth=0.3, 
                  cbar_kws={'shrink' : 0.2, 'aspect' : 20}, cmap='Blues', square=False, xticklabels=True, yticklabels=True, ax=ax)
    cbar=ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)    
    ax.tick_params(axis='x', labelsize=11, rotation=90)
    ax.tick_params(axis='y', labelsize=9)
    ax.set_ylabel('Metagenomic MAGs', fontsize=20) 
    ax.set_xlabel('ELP classes', fontsize=20)
    ax.xaxis.tick_top()    
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_All_MAGs_heatmap_{norm}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_All_MAGs_heatmap_{norm}.svg'), dpi=300)
    
    fig=plt.figure(figsize=(7,16), dpi=300) # 13*18 for All classes; 7*16 for types;
    ax=fig.add_subplot(111)
    
    # Draw heatmap with less rounded values.
    ax=sb.heatmap(input_df, annot=True, fmt=".1f", annot_kws={'fontsize' : 7, 'fontweight' : 'normal'}, linewidth=0.3, 
                  cbar_kws={'shrink' : 0.2, 'aspect' : 20}, cmap='Blues', square=False, xticklabels=True, yticklabels=True, ax=ax)
    cbar=ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=13)    
    ax.tick_params(axis='x', labelsize=11, rotation=90)
    ax.tick_params(axis='y', labelsize=9)
    ax.set_ylabel('Metagenomic MAGs', fontsize=20) 
    ax.set_xlabel('ELP classes', fontsize=20)
    ax.xaxis.tick_top()    
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_All_MAGs_heatmap_{norm}_not_rounded.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_All_MAGs_heatmap_{norm}_not_rounded.svg'), dpi=300)    
    
    return


#######
# Run PCA analysis on gene/domain frequencies.
#######

def draw_PCA_all(ELP_cat, All_source_ELPs_df, norm, Sponge_assoc_MAGs, output_path):
    
    model=pca() 
    results=model.fit_transform(All_source_ELPs_df) 
    
    All_genome_names=list(All_source_ELPs_df.index)
    
    Genome_names_to_label_ar=[]
    Genome_colors_to_label_ar=[]
    Genome_alpha_to_label_ar=[]
    
    for genome_name in All_genome_names:
        
        if genome_name in Sponge_assoc_MAGs:
            
            Genome_names_to_label_ar.append(genome_name.split('_')[2])
            Genome_colors_to_label_ar.append('black')
            Genome_alpha_to_label_ar.append(1)
            
        else:
            
            Genome_names_to_label_ar.append('')
            Genome_colors_to_label_ar.append('blue')
            Genome_alpha_to_label_ar.append(0.5) 
            
    Genome_names_to_label_ser=np.array(Genome_names_to_label_ar)
    
    fig, ax=model.biplot(SPE=True, HT2=True, n_feat=5, labels=Genome_names_to_label_ser, 
                         alpha=Genome_alpha_to_label_ar, c=Genome_colors_to_label_ar, linewidths=0.2,
                         legend=False, s=20, dpi=300, grid=False, 
                         arrowdict={'fontsize': 5}, title=None, fontsize=5)
    ax.set_title(None)
    
    ax.tick_params(axis='both', labelsize=10)
    
    xlabel=ax.get_xlabel()
    ax.set_xlabel(xlabel, fontsize=10)
    
    ylabel=ax.get_ylabel()
    ax.set_ylabel(ylabel, fontsize=10)
    
    fig.set_size_inches(3, 3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_All_MAGs_biplot_{norm}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_All_MAGs_biplot_{norm}.svg'), dpi=300)    
    
    return


def plot_ELP_distrib(FL_freqs, SP_freqs, ELP_cat, ELP_name, norm, output_path):
    
    Full_ELPs_data=FL_freqs+SP_freqs
    
    fig=plt.figure(figsize=(4,2), dpi=300) 
    ax=fig.add_subplot(111)
    
    ax.hist(Full_ELPs_data, histtype=u'step', color='#84c2ff', alpha=0.9, linewidth=3)
      
    for SP_freq in SP_freqs:
        ax.axvline(SP_freq, color='black', linestyle='--', linewidth=1.5, alpha=0.6)

    ax.set_xlabel(f'Number of ELPs ({ELP_name})', size=12)
    ax.set_ylabel('Number of MAGs', size=12)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
   
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_dist_for_{ELP_name}_all_MAGs_{norm}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_dist_for_{ELP_name}_all_MAGs_{norm}.svg'), dpi=300)    
    
    
    return


def calc_stat_draw_distribs(ELP_cat, FL_ELPs_df, SP_MAGs_ELPs_df, norm, output_path):
    
    print(f'Statistical tests for {ELP_cat}')
    
    Average_ELP_freq_FL=list(FL_ELPs_df.mean(axis=1))
    Average_ELP_freq_SP=list(SP_MAGs_ELPs_df.mean(axis=1))
    FL_SP_stat=sp.stats.ttest_ind(Average_ELP_freq_FL, Average_ELP_freq_SP)
    print(f'T-test p-value for {ELP_cat} FL ({np.mean(Average_ELP_freq_FL)}) vs SP ({np.mean(Average_ELP_freq_SP)}) : {FL_SP_stat[1]}')
    
    ELP_names_ar=list(FL_ELPs_df.columns)
    
    for ELP_name in ELP_names_ar:
        FL_freqs=list(FL_ELPs_df[ELP_name])
        SP_freqs=list(SP_MAGs_ELPs_df[ELP_name])
        FL_SP_stat=sp.stats.ttest_ind(FL_freqs, SP_freqs)
        
        print(f'T-test p-value for {ELP_name} FL ({np.mean(FL_freqs)}) vs SP ({np.mean(SP_freqs)}) : {FL_SP_stat[1]}')
        if FL_SP_stat[1]<(0.05):
            plot_ELP_distrib(FL_freqs, SP_freqs, ELP_cat, ELP_name, norm, output_path)
        
    
    return


def plot_bins_stats(CheckM_info_df, Sponge_assoc_MAGs, output_path):
    
    # Plot completeness vs contamination.
    fig=plt.figure(figsize=(4.5, 4.5))
    gs=fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                        left=0.1, right=0.9, bottom=0.1, top=0.9,
                        wspace=0.05, hspace=0.05) 
    
    ax=fig.add_subplot(gs[1, 0])
    ax_histx=fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy=fig.add_subplot(gs[1, 1], sharey=ax)
    
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)   
    
    ax.scatter(CheckM_info_df['Completeness'], CheckM_info_df['Contamination'], color='k', label='Free-living')
    CheckM_info_SA_MAGs_df=CheckM_info_df.loc[Sponge_assoc_MAGs]
    ax.scatter(CheckM_info_SA_MAGs_df['Completeness'], CheckM_info_SA_MAGs_df['Contamination'], s=45, color='r', label='Sponge-associated')
    ax.set_xlabel('Completeness, %')
    ax.set_ylabel('Contamination, %')
    
    ax_histx.hist(CheckM_info_df['Completeness'], bins=15, color='grey')
    for SP_MAG_name in Sponge_assoc_MAGs:
        ax_histx.axvline(CheckM_info_df.loc[SP_MAG_name]['Completeness'], color='red', linestyle='--', linewidth=1.5, alpha=0.6)
        
    ax_histy.hist(CheckM_info_df['Contamination'], bins=15, orientation='horizontal', color='grey')  
    for SP_MAG_name in Sponge_assoc_MAGs:
        ax_histy.axhline(CheckM_info_df.loc[SP_MAG_name]['Contamination'], color='red', linestyle='--', linewidth=1.5, alpha=0.6)    
    
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax_histx.spines["top"].set_visible(False)
    ax_histx.spines["right"].set_visible(False)
    ax_histy.spines["top"].set_visible(False)
    ax_histy.spines["right"].set_visible(False)    
 
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'Sponge_associated_and_free_living_bins_and_MAGs_comp_vs_cont.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'Sponge_associated_and_free_living_bins_and_MAGs_comp_vs_cont.svg'), dpi=300)        
    
    
    # Plot distribution of bin/MAG length.
    fig=plt.figure(figsize=(4,2), dpi=300) 
    ax=fig.add_subplot(111)
    
    ax.hist(CheckM_info_df['Genome_Size']/1000000, histtype=u'step', color='k', alpha=0.8, linewidth=2.5)
      
    for SP_MAG_name in Sponge_assoc_MAGs:
        ax.axvline(CheckM_info_df.loc[SP_MAG_name]['Genome_Size']/1000000, color='red', linestyle='--', linewidth=1.5, alpha=0.6)

    ax.set_xlabel(f'Genome length, Mb', size=12)
    ax.set_ylabel('Number of MAGs', size=12)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
   
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'Sponge_associated_and_free_living_bins_and_MAGs_size.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'Sponge_associated_and_free_living_bins_and_MAGs_size.svg'), dpi=300)    
    
    
    return




def wrapper_function_ELP_vis(elps_in_FL_MAGs_path, elps_in_SA_MAGs_path, checkM_bins_info, output_path):
    
    if os.path.isdir(output_path)!=True:
        os.mkdir(output_path)
    
    # Select the category of ELPs to compare. Options:
    ELP_cat="All_ELPs_types_genes"           # Counts for ELP-encoding genes.
    #ELP_cat="All_ELPs_types_domains"          # Counts fro ELP domains.
    #ELP_cat="Exported_ELPs_types_genes"      # Counts for ELP-encoding genes with predicted high-confident exporting signals.
    #ELP_cat="Exported_ELPs_types_domains"    # Counts for ELP domains with predicted high-confident exporting signals.
    
    
    # Read FL data.
    FL_ELPs_df=pd.read_excel(os.path.join(elps_in_FL_MAGs_path, f'{ELP_cat}_counts.xlsx'), sheet_name='Sheet1', header=0, index_col=0)
    #print(FL_ELPs_df)
    
    # Read SP data.
    SP_ELPs_df=pd.read_excel(os.path.join(elps_in_SA_MAGs_path, f'{ELP_cat}_counts.xlsx'), sheet_name='Sheet1', header=0, index_col=0)
    Sponge_assoc_MAGs=['I_palmata_OTU1_sp', 'H_sitiens_OTU3_sp', 'H_sitiens_OTU7_sp', 'H_sitiens_OTU9_sp', 'H_sitiens_OTU14_sp', 'H_panicea_OTU4_sp', 'H_panicea_OTU23_sp']
    SP_MAGs_ELPs_df=SP_ELPs_df.loc[Sponge_assoc_MAGs]
    #print(SP_MAGs_ELPs_df)  
    
    # Concatenate data and draw heatmap.
    norm='frequences'
    All_source_ELPs_df=pd.concat([FL_ELPs_df, SP_MAGs_ELPs_df])
    draw_heatmap_all(ELP_cat, All_source_ELPs_df, norm, output_path)
    
    # Run PCA analysis on gene/domain frequencies.
    draw_PCA_all(ELP_cat, All_source_ELPs_df, norm, Sponge_assoc_MAGs, output_path)
    
    # Calculate statistics for the frequency of ELPs, draw distributions.
    calc_stat_draw_distribs(ELP_cat, FL_ELPs_df, SP_MAGs_ELPs_df, norm, output_path)
    
    # Read CheckM info for metagenomic bins.
    CheckM_info_df=pd.read_excel(checkM_bins_info, sheet_name='Sheet1', header=0, index_col=0)
    
    # Draw bins and MAGs stats.
    plot_bins_stats(CheckM_info_df, Sponge_assoc_MAGs, output_path)
    
    # Normalize ELP counts by bins\MAGs length.
    FL_ELPs_df_norm, SP_MAGs_ELPs_df_norm=norm_ELP_counts(FL_ELPs_df, SP_MAGs_ELPs_df, CheckM_info_df)
    
    # Concatenate data and draw heatmap with normalized values.
    All_source_ELPs_df_norm=pd.concat([FL_ELPs_df_norm, SP_MAGs_ELPs_df_norm])
    List_of_bins=list(FL_ELPs_df.index)+Sponge_assoc_MAGs
    All_source_ELPs_df_norm=All_source_ELPs_df_norm.reindex(List_of_bins)
    norm='normalized'
    draw_heatmap_all(ELP_cat, All_source_ELPs_df_norm, norm, output_path)   
    
    # Run PCA analysis on gene/domain frequencies.
    draw_PCA_all(ELP_cat, All_source_ELPs_df_norm, norm, Sponge_assoc_MAGs, output_path)    
    
    # Calculate statistics for the normalized frequency of ELPs, draw distributions.
    calc_stat_draw_distribs(ELP_cat, FL_ELPs_df_norm, SP_MAGs_ELPs_df_norm, norm, output_path)    
    
    return



wrapper_function_ELP_vis(ELPs_in_FL_MAGs_path, ELPs_in_SA_MAGs_path, CheckM_bins_info, Output_path)