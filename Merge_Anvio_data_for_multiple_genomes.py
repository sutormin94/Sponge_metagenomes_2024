###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2025##
##Prepares data obtained from Anvio for downstream analysis.

# Takes Anvio output from multiple datasets and merges them into one large dataframe, 
# also adds taurine and sulfoacetate metabolic pathways, annotated separately.
###############################################

#######
#Packages to be imported.
#######

import os
import pandas as pd
import numpy as np


# Path to anvio datasets.
Anvio_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\\Source_data\\Metabolic_analysis_SAB_MAGs_and_related_species\\Anvio_data\\"

# Path to sulfoacetate and taurine metabolic pathways annotation.
Tau_sulf_path="C:\\Users\\sutor\\OneDrive\\ThinkPad_working\\Sutor\\Science\\Spongy\\Reports_and_burocracy\\Manuscripts\\Manuscript_2023\\Review\\Rusanova,2025,Supplementary_Tables.xlsx"
# Spreadsheet with sulfoacetate and taurine metabolic pathways annotation.
Tau_sulf_sheet="S13"

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\\Source_data\\Metabolic_analysis_SAB_MAGs_and_related_species\\"

