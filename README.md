# Sponge_metagenomes_2024
 Scripts used for analysis and visualization of sponge metagenomic data.
 
## Preparation of symbiont plots

### Description 

 **`Symbiont_plot.py`**
  - Was used to make symbiont plots based on 16S counts.

### Usage  

 Provide OTU table with 16S counts (Additional_data\OTU_counts_and_taxa.xlsx) and a pair of 16S datasets for comparison (e.g., "H. sitiens 2016 1" and "Marine water 2016").


## Script for preparation of full-length 16S sequences for multiple alignment and their annotation for iTol

 **`Prepare_16S_dataset_for_multiple_alignment.py`**
   - Was used to merge blast results obtained from different sources (NCBI nt, NCBI rRNA, SILVA), cluster with mmseqs2 (module 1), and add biome information from NCBI records (module 2).
 
### Usage  

 Use module 1 for merging and clustering of initial sequences identified by blastn. Use module 2 to prepare iTol metadata and decorate phylogenetic trees.
16S alignments (fasta), metadata (csv), and resultant phylogenetic trees (newick) can be found in `Additional_data\Full_length_16S_trees\`.


## Scripts for visualization of metabolic features

### Description  

 **`Draw_metabolic_heatmaps.py`**
   - Was used to make heatmaps for amino acids and vitamins biosynthsis and transport, and for secretion and transport systems abundance.
   
### Usage  

 Provide script with the excel table containing pathways completeness `Additional_data\Amino_acids_vitamins_sec_systems.xlsx`
 
 
## Scripts for ELP annotation and analysis 

Scripts used for identifying **Eukaryotic-Like Proteins (ELPs)** in genome annotation files.  

### Description  

1. **`Run_IPS_and_find_ELPs.sh`**  
   - Executes **MetaGeneMark-2** to predict coding sequences (CDS) from genome files.  
   - Uses **InterProScan** to annotate files with amino acid sequences.  

2. **`Analyse_ELP_form_MAGs.py`** 
   - Contains three modules for a) aggregation of data predicted by InterProScan, b) adding the signalP annotation to predicted ALP-containing proteins,
c) statistical analysis and visualization of ELP counts.

### Usage  

1. **Predict and annotate genes**  
   Begin with `Run_IPS_and_find_ELPs.sh` to process genome files (genomic fasta) and generate annotated outputs (`Additional_data\All_ELPs_found_annotation.tsv`).  

2. **Search for ELPs and analyze**  
   Use `Analyse_ELP_form_MAGs.py` to identify potential ELPs (domains listed in `Additional_data\ELPs_list.tsv`) based on the annotation results, add signalP data (`Additional_data\All_ELPs_found_seq_SignalP6.txt`) and visualize the resulting matrices and distribution.  
   Provide the script with ELP count files (xlsx) in `Additional_data\ELP_data_all_FL_genomes` for free-living bins and `Additional_data\ELP_data_all_SA_genomes` for sponge-associated bins.

These scripts streamline the process of identifying and analyzing ELPs across multiple genomes.


## Trimm sequence annotation for clinker

### Description  

**`Trim_gff_for_clinker.py`** was used to trim sequnce (fasta) and its annotation (gff) to a particular region specified in the gff and to
convert coordinates to relative coordinates for the extracted fragment. Useful to remove unecessary flanks when using clinker.

### Usage  

 Provide script with a sequence (fasta) and the annotation of a region of interest (gff).
