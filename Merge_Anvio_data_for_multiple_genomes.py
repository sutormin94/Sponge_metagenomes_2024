###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2025##
##Prepares data obtained from Anvio for downstream analysis.

# Takes Anvio output from multiple datasets and merges them into one large dataframe, 
# also adds information on taurine and sulfoacetate metabolic pathways, annotated separately.
###############################################

#######
#Packages to be imported.
#######

import os
import pandas as pd
import numpy as np



def read_datasets_merge(anvio_path, tau_sulf_path, tau_sulf_sheet, output_path):
    
    # Initiate dataframes for merging.
    
    all_genomes_metab_stepwise_df=pd.DataFrame()
    all_genomes_metab_pathwise_df=pd.DataFrame()
    
    description_columns=['module_name', 'module_class', 'module_category', 'module_subcategory', 'module_definition']
    all_genomes_descr_col_df=pd.DataFrame()
    
    # Read and merge Anvio results.
    
    for anvio_set_name in os.listdir(anvio_path):
        
        anvio_set_full_path=os.path.join(anvio_path, anvio_set_name)
        genome_name=anvio_set_name.split('.modules')[0]
        
        print(f'Merging {genome_name}')
        
        anvio_full_df=pd.read_csv(anvio_set_full_path, sep='\t', header=0, index_col="module")
        
        # Merge description columns.
        anvio_descr_col_df=anvio_full_df[description_columns] 
        new_rows=anvio_descr_col_df[~anvio_descr_col_df.index.isin(all_genomes_descr_col_df.index)]
        all_genomes_descr_col_df=pd.concat([all_genomes_descr_col_df, new_rows], ignore_index=False)       
        
        # Merge completeness column.
        anvio_stepwise_df=anvio_full_df[['stepwise_module_completeness']]
        anvio_stepwise_df=anvio_stepwise_df.rename(columns={'stepwise_module_completeness': genome_name})
        all_genomes_metab_stepwise_df=pd.merge(all_genomes_metab_stepwise_df, anvio_stepwise_df, left_index=True, right_index=True, how='outer') 
        
        anvio_pathwise_df=anvio_full_df[['pathwise_module_completeness']]
        anvio_pathwise_df=anvio_pathwise_df.rename(columns={'pathwise_module_completeness': genome_name})
        all_genomes_metab_pathwise_df=pd.merge(all_genomes_metab_pathwise_df, anvio_pathwise_df, left_index=True, right_index=True, how='outer')   
        
    
    # Merge completeness columns with description columns.
    all_genomes_metab_stepwise_df=pd.merge(all_genomes_metab_stepwise_df, all_genomes_descr_col_df, left_index=True, right_index=True, how='outer')
    all_genomes_metab_pathwise_df=pd.merge(all_genomes_metab_pathwise_df, all_genomes_descr_col_df, left_index=True, right_index=True, how='outer')
    
    # Add sulfoacetate and taurine metabolic pathways completeness.
    
    tau_sulf_df=pd.read_excel(tau_sulf_path, sheet_name=tau_sulf_sheet, header=0, index_col=0)
    tau_sulf_df=tau_sulf_df.T
    all_genomes_metab_stepwise_df=pd.concat([all_genomes_metab_stepwise_df, tau_sulf_df], ignore_index=False) 
    all_genomes_metab_pathwise_df=pd.concat([all_genomes_metab_pathwise_df, tau_sulf_df], ignore_index=False)
    
    # Write resultant merged dataframes.
    all_genomes_metab_stepwise_df=all_genomes_metab_stepwise_df.fillna(0)
    all_genomes_metab_stepwise_df=all_genomes_metab_stepwise_df.T
    all_genomes_metab_stepwise_df.index.name='Genome'
    Stepwise_outpath=os.path.join(output_path, 'Anvio_and_tau_sulf_metab_stepwise_completeness_data.xlsx')
    all_genomes_metab_stepwise_df.to_excel(Stepwise_outpath, index=True)
    
    all_genomes_metab_pathwise_df=all_genomes_metab_pathwise_df.fillna(0)
    all_genomes_metab_pathwise_df=all_genomes_metab_pathwise_df.T
    all_genomes_metab_pathwise_df.index.name='Genome'
    Pathwise_outpath=os.path.join(output_path, 'Anvio_and_tau_sulf_metab_pathwise_completeness_data.xlsx')
    all_genomes_metab_pathwise_df.to_excel(Pathwise_outpath, index=True)
    
    return


# Path to anvio datasets.
Anvio_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\\Source_data\\Metabolic_analysis_SAB_MAGs_and_related_species\\Anvio_data\\"

# Path to sulfoacetate and taurine metabolic pathways annotation.
Tau_sulf_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\\Source_data\\Metabolic_analysis_SAB_MAGs_and_related_species\\Taurine_sulfoacetate_metab_pathways_completeness.xlsx"
# Spreadsheet with sulfoacetate and taurine metabolic pathways annotation.
Tau_sulf_sheet="Sheet1"

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\\Source_data\\Metabolic_analysis_SAB_MAGs_and_related_species\\"


read_datasets_merge(Anvio_path, Tau_sulf_path, Tau_sulf_sheet, Output_path)