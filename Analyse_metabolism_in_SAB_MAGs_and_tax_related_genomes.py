###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2025##
##Analyse completeness of metabolic pathways in bacterial MAGs:
##compares SAB MAGs and taxonomically related genomes.

###############################################

#######
#Packages to be imported.
#######

import os
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import scipy as sp
from pca import pca
from sklearn.decomposition import PCA as sklearnPCA


# -------------------- Parameters --------------------
Metabolism_folder_path = "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Scripts\Sponge_metagenomes_2024\Source_data\Metabolic_analysis_SAB_MAGs_and_related_species\\"  # Set this to the folder containing your 4 Excel files

# The 2 Excel files to process:
Metabolism_filenames = [
    "Anvio_and_tau_sulf_metab_pathwise_completeness_data.xlsx",
]

CheckM_bins_info = os.path.join(Metabolism_folder_path, "MAGs_and_bins_CheckM.xlsx")

Output_path = os.path.join(Metabolism_folder_path, "Metabolism_stat_and_pictures")
os.makedirs(Output_path, exist_ok=True)

if not os.path.exists(Output_path):
    os.makedirs(Output_path)

# OTU main genomes and order (as before)
OTU_main_genomes = {
    "OTU4": "Ca_Halichondribacter_symbioticus_OTU4",
    "OTU23": "Ca_Yagmuria_paniceus_OTU23",
    "OTU3": "Ca_Ahtobacter_symbioticus_OTU3",
    "OTU7": "Ca_Vellamobacter_salmiensis_OTU7",
    "OTU9": "Ca_Vienanmeria_sitiensis_OTU9",
    "OTU14": "Ca_Sampovibrio_pertsovi_OTU14",
    "OTU1": "Ca_Eurynomebacter_symbioticus_OTU1"
}
OTU_order = ["OTU4", "OTU23", "OTU3", "OTU7", "OTU9", "OTU14", "OTU1"]
Sponge_assoc_MAGs = list(OTU_main_genomes.values())


# -------------------- Functions --------------------

def load_and_prepare_metab_table(filepath):
    """
    Load a dataset table and add 'Group' column.
    """
    df = pd.read_excel(filepath, index_col=0, header=0)
    df["Genome"] = df.index
    df["Group"] = np.where(df["Habitat"].str.strip() == "Free living", "FL", "SP")
    module_description_df = df.loc[['module_name', 'module_class', 'module_category', 'module_subcategory', 'module_definition']]
    df = df.drop(['module_name', 'module_class', 'module_category', 'module_subcategory', 'module_definition'])
    return df, module_description_df


def get_sorted_genomes(df, otu, main_genome):
    genomes = []
    mask_otu = df["Relative"] == otu
    genomes += df.loc[mask_otu & (df["Genome"] == main_genome), "Genome"].tolist()
    genomes += df.loc[mask_otu & (df["Group"] == "SP") & (df["Genome"] != main_genome), "Genome"].tolist()
    genomes += df.loc[mask_otu & (df["Group"] == "FL"), "Genome"].tolist()
    return genomes


def make_final_genome_order(df):
    genomes = []
    for otu in OTU_order:
        main_genome = OTU_main_genomes[otu]
        genomes += get_sorted_genomes(df, otu, main_genome)
    return genomes


def draw_heatmap(input_df, row_colors, title, options, output_path):
    """
    Plot heatmap with seaborn clustermap, with colored rows by group.
    Clustering disabled to preserve ordering.
    """   
    
    cg = sb.clustermap(
        input_df,
        row_colors=row_colors,
        cmap='YlGnBu',
        linewidths=0.2,
        row_cluster=False,
        col_cluster=options[0],
        yticklabels=True,
        xticklabels=True,
        cbar_kws={'aspect' : 20, 'shrink' : 0.5}
    )  
    cg.fig.set_size_inches(options[1][0], options[1][1])
    
    plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=options[2])
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), fontsize=options[2])
    cg.ax_heatmap.set_xlabel("Metabolic pathways", fontsize=options[3])
    cg.ax_heatmap.set_ylabel("Genomes", fontsize=options[3])

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{title}_heatmap.png'), dpi=800)
    plt.savefig(os.path.join(output_path, f'{title}_heatmap.svg'), dpi=800)
    plt.close()
    
    return


def draw_PCA_all(ELP_cat, all_data_df, all_data_df_norm, sponge_assoc_MAGs, habitat_color_dict, output_path):
    
    fig, plot_pca=plt.subplots(1,2,figsize=(6,3), dpi=300)
    
    All_genome_names=list(all_data_df.index)
    
    model=pca() 
    results=model.fit_transform(all_data_df_norm)        
    
    Genome_names_to_label_ar=[]
    Genome_colors_to_label_ar=[]
    Genome_alpha_to_label_ar=[]
    Genome_tax_colors_to_label=[]
    
    for genome_name in All_genome_names:
        
        if genome_name in sponge_assoc_MAGs:
            
            Genome_names_to_label_ar.append(genome_name.split('_')[-1])
            point_source_color=habitat_color_dict[all_data_df.loc[genome_name, "Habitat"]]
            Genome_colors_to_label_ar.append(point_source_color)
            point_tax_color=all_data_df.loc[genome_name, "Tax group color"]
            Genome_tax_colors_to_label.append(point_tax_color)            
            Genome_alpha_to_label_ar.append(1)
            
        else:
            
            Genome_names_to_label_ar.append("")
            point_source_color=habitat_color_dict[all_data_df.loc[genome_name, "Habitat"]]
            Genome_colors_to_label_ar.append(point_source_color)
            point_tax_color=all_data_df.loc[genome_name, "Tax group color"]
            Genome_tax_colors_to_label.append(point_tax_color)            
            Genome_alpha_to_label_ar.append(1) 
            
    Genome_names_to_label_ser=np.array(Genome_names_to_label_ar)
    
    print(f'Making PCA for all genomes, number of genomes/features {all_data_df_norm.shape}')
    
    # Color points by isolation source.
    model.biplot(SPE=False, HT2=False, n_feat=5, labels=Genome_names_to_label_ser, 
                 alpha=Genome_alpha_to_label_ar, c=Genome_colors_to_label_ar, linewidths=0.2,
                 legend=False, s=20, dpi=300, grid=False, 
                 arrowdict={'fontsize': 5}, title=None, fontsize=5, ax=plot_pca[0])
    plot_pca[0].set_title(None)
    plot_pca[0].tick_params(axis='both', labelsize=10)
    xlabel=plot_pca[0].get_xlabel()
    plot_pca[0].set_xlabel(xlabel, fontsize=10)
    ylabel=plot_pca[0].get_ylabel()
    plot_pca[0].set_ylabel(ylabel, fontsize=10)          
    
    # Color points by taxonomy group.
    model.biplot(SPE=False, HT2=False, n_feat=5, labels=Genome_names_to_label_ser, 
                 alpha=Genome_alpha_to_label_ar, c=Genome_tax_colors_to_label, linewidths=0.2,
                 legend=False, s=20, dpi=300, grid=False, 
                 arrowdict={'fontsize': 5}, title=None, fontsize=5, ax=plot_pca[1])
    plot_pca[1].set_title(None)
    plot_pca[1].tick_params(axis='both', labelsize=10)
    xlabel=plot_pca[1].get_xlabel()
    plot_pca[1].set_xlabel(xlabel, fontsize=10)
    ylabel=plot_pca[1].get_ylabel()
    plot_pca[1].set_ylabel(ylabel, fontsize=10)        
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_all_PCA_bipl.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_all_PCA_bipl.svg'), dpi=300)     
    
    return


def draw_PCA_groupwise(ELP_cat, all_data_df, all_data_df_norm, otu_main_genomes, otu_order, habitat_color_dict, output_path):
    
    num_plots=len(otu_order)
    num_columns=3
    if num_plots%num_columns==0:
        num_rows=num_plots//num_columns
    else:
        num_rows=(num_plots//num_columns)+1
    
    fig, plot_pca=plt.subplots(num_rows,num_columns,figsize=(6,6), dpi=300)
    
    for i in range(len(otu_order)):
        
        otu=otu_order[i]
        row_id=i//num_columns
        column_id=i%num_columns
        
        data_df_otu=all_data_df[all_data_df["Relative"]==otu]
        All_genome_names=list(data_df_otu.index)
        data_df_norm_otu=all_data_df_norm.loc[All_genome_names]
        
        model=pca() 
        results=model.fit_transform(data_df_norm_otu)        
        
        main_genome=otu_main_genomes[otu]
        
        Genome_names_to_label_ar=[]
        Genome_colors_to_label_ar=[]
        Genome_alpha_to_label_ar=[]
        
        for genome_name in All_genome_names:
            
            if genome_name==main_genome:
                
                Genome_names_to_label_ar.append(otu)
                point_color=habitat_color_dict[all_data_df.loc[genome_name, "Habitat"]]
                Genome_colors_to_label_ar.append(point_color)
                Genome_alpha_to_label_ar.append(1)
                
            else:
                
                Genome_names_to_label_ar.append("")
                point_color=habitat_color_dict[all_data_df.loc[genome_name, "Habitat"]]
                Genome_colors_to_label_ar.append(point_color)
                Genome_alpha_to_label_ar.append(1) 
                
        Genome_names_to_label_ser=np.array(Genome_names_to_label_ar)
        
        print(f'Making PCA for {otu}, number of genomes/features {data_df_norm_otu.shape}')
        
        model.biplot(SPE=False, HT2=False, n_feat=5, labels=Genome_names_to_label_ser, 
                     alpha=Genome_alpha_to_label_ar, c=Genome_colors_to_label_ar, linewidths=0.2,
                     legend=False, s=20, dpi=300, grid=False, 
                     arrowdict={'fontsize': 5}, title=None, fontsize=5, ax=plot_pca[row_id,column_id])
        plot_pca[row_id,column_id].set_title(otu)
        
        plot_pca[row_id,column_id].tick_params(axis='both', labelsize=10)
        
        xlabel=plot_pca[row_id,column_id].get_xlabel()
        plot_pca[row_id,column_id].set_xlabel(xlabel, fontsize=10)
        
        ylabel=plot_pca[row_id,column_id].get_ylabel()
        plot_pca[row_id,column_id].set_ylabel(ylabel, fontsize=10)        
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_MAGs_relwise_bipl.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_MAGs_relwise_bipl.svg'), dpi=300)  
    
    return


def plot_module_compl_distrib(compl_df, dataset_name, otu_ID, compl_metrics, output_path):
    
    # Plot distribution of module completeness.
    module_compl_ar=list(compl_df[compl_df['Relative']==otu_ID][f"{compl_metrics} module completeness"])
    otu_module_compl=compl_df[compl_df['Relative']==otu_ID].loc[OTU_main_genomes[otu_ID], f"{compl_metrics} module completeness"]    

    plt.figure(figsize=(4, 2), dpi=300)
    ax = plt.gca()
    ax.hist(module_compl_ar, bins=30, histtype='step', color='#84c2ff', alpha=0.9, linewidth=3)
    ax.axvline(otu_module_compl, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.set_xlabel(f'{compl_metrics} module completeness', fontsize=12)
    ax.set_ylabel('Number of MAGs', fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    safe_name = f'{dataset_name}_{compl_metrics}_dist_{otu_ID}'.replace(" ", "_")
    plt.savefig(os.path.join(output_path, f'{safe_name}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{safe_name}.svg'), dpi=300)
    plt.close()
    
    # Plot genome completeness vs module completeness.
    genome_compl_ar=list(compl_df[compl_df['Relative']==otu_ID][f"Genome completeness"])
    otu_genome_compl=compl_df[compl_df['Relative']==otu_ID].loc[OTU_main_genomes[otu_ID], f"Genome completeness"]    
    
    fig = plt.figure(figsize=(4.5, 4.5))
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.15, right=0.95, bottom=0.15, top=0.95,
                          wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    ax.scatter(genome_compl_ar, module_compl_ar, color='k')
    ax.scatter(otu_genome_compl, otu_module_compl, s=100, color='r', label=otu_ID, zorder=10)

    ax.set_xlabel('Genome completeness, %', fontsize=10)
    ax.set_ylabel(f'{compl_metrics} module completeness, %', fontsize=10)
    ax.tick_params(axis='both', labelsize=10)

    ax_histx.hist(genome_compl_ar, bins=15, color='grey')
    ax_histx.axvline(otu_genome_compl, color='red', linestyle='--', linewidth=1.5, alpha=0.6)
    ax_histx.tick_params(axis='both', labelsize=10)

    ax_histy.hist(module_compl_ar, bins=15, orientation='horizontal', color='grey')
    ax_histy.axhline(otu_module_compl, color='red', linestyle='--', linewidth=1.5, alpha=0.6)
    ax_histy.tick_params(axis='both', labelsize=10)

    for axis in [ax, ax_histx, ax_histy]:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{dataset_name}_{compl_metrics}_cg_cm_{otu_ID}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{dataset_name}_{compl_metrics}_cg_cm_{otu_ID}.svg'), dpi=300)
    plt.close()
    
    return


def plot_ELP_distrib(fl_freqs, sp_freqs, ELP_cat, ELP_name, norm, output_path):
    all_freqs = fl_freqs + sp_freqs
    plt.figure(figsize=(4, 2), dpi=300)
    ax = plt.gca()
    ax.hist(all_freqs, bins=30, histtype='step', color='#84c2ff', alpha=0.9, linewidth=3)
    for sp_val in sp_freqs:
        ax.axvline(sp_val, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.set_xlabel(f'Number of ELPs ({ELP_name})', fontsize=12)
    ax.set_ylabel('Number of MAGs', fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    safe_name = f'{ELP_cat}_dist_for_{ELP_name}_all_MAGs_{norm}'.replace(" ", "_")
    plt.savefig(os.path.join(output_path, f'{safe_name}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{safe_name}.svg'), dpi=300)
    plt.close()
    return


def calc_stat_draw_distribs(ELP_cat, fl_df, sp_df, norm, output_path):
    features = fl_df.columns.tolist()
    for feature in features:
        fl_freqs = fl_df[feature].dropna().tolist()
        sp_freqs = sp_df[feature].dropna().tolist()
        stat_result = sp.stats.ttest_ind(fl_freqs, sp_freqs, equal_var=False)  # Welch's t-test
        p_val = stat_result.pvalue
        print(f"T-test p-value for {feature} {norm} FL (mean={np.mean(fl_freqs):.2f}) vs SP (mean={np.mean(sp_freqs):.2f}): {p_val:.4f}")
        if p_val < 0.05:
            plot_ELP_distrib(fl_freqs, sp_freqs, ELP_cat, feature, norm, output_path)
    return


def plot_bins_stats(checkm_df, sponge_genomes, output_path):
    fig = plt.figure(figsize=(4.5, 4.5))
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    fl_genomes = checkm_df.index.difference(sponge_genomes)
    ax.scatter(checkm_df.loc[fl_genomes, 'Completeness'], checkm_df.loc[fl_genomes, 'Contamination'],
               color='k', label='Free-living')
    ax.scatter(checkm_df.loc[sponge_genomes, 'Completeness'], checkm_df.loc[sponge_genomes, 'Contamination'],
               s=45, color='r', label='Sponge-associated')

    ax.set_xlabel('Completeness, %')
    ax.set_ylabel('Contamination, %')

    ax_histx.hist(checkm_df['Completeness'], bins=15, color='grey')
    for sp_genome in sponge_genomes:
        ax_histx.axvline(checkm_df.loc[sp_genome]['Completeness'], color='red', linestyle='--', linewidth=1.5, alpha=0.6)

    ax_histy.hist(checkm_df['Contamination'], bins=15, orientation='horizontal', color='grey')
    for sp_genome in sponge_genomes:
        ax_histy.axhline(checkm_df.loc[sp_genome]['Contamination'], color='red', linestyle='--', linewidth=1.5, alpha=0.6)

    for axis in [ax, ax_histx, ax_histy]:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'Completeness_vs_Contamination.png'), dpi=300)
    plt.savefig(os.path.join(output_path, 'Completeness_vs_Contamination.svg'), dpi=300)
    plt.close()

    plt.figure(figsize=(4, 2), dpi=300)
    ax = plt.gca()
    genome_sizes_mb = checkm_df['Genome_Size'] / 1e6
    ax.hist(genome_sizes_mb, bins=30, histtype='step', color='k', alpha=0.8, linewidth=2.5)
    for sp_genome in sponge_genomes:
        ax.axvline(genome_sizes_mb.loc[sp_genome], color='red', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.set_xlabel('Genome Size (Mb)', fontsize=12)
    ax.set_ylabel('Number of MAGs', fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'Genome_Size_Distribution.png'), dpi=300)
    plt.savefig(os.path.join(output_path, 'Genome_Size_Distribution.svg'), dpi=300)
    plt.close()
    return


# -------------------- Main wrapper handling multiple files --------------------

def wrapper_metabolism_vis_multiple(metabolism_folder_path, metabolism_filenames, checkm_file, output_dir):
    """
    Main wrapper function that processes several Excel files in a folder.
    It loads each file, assigns groups based on Habitat, orders genomes by OTU and Relative, 
    performs PCA, statistical tests, and draws heatmaps.
    Output plots and results are saved to a directory which is created if missing.

    Parameters:
    -----------
    metabolism_folder_path : str
        Path to the folder containing the 2 Excel files with completeness of metabolic pathways.
    metabolism_filenames : list of str
        List of the 2 Excel filenames to process.
    checkm_file : str
        Path to the CheckM file containing genome statistics.
    output_dir : str
        Directory path where output results and figures will be saved.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkm_data = pd.read_excel(checkm_file)
    # Ensure column names are standardized
    checkm_data.rename(columns={"Name": "Genome", "GenomeSize": "Genome_Size"}, inplace=True)
    checkm_data.set_index("Genome", inplace=True)

    for filename in metabolism_filenames:
        print(f"Processing file: {filename}")
        fullpath = os.path.join(metabolism_folder_path, filename)
        metab_df, metab_desc_df = load_and_prepare_metab_table(fullpath)

        meta_cols = ["Genome", "Relative", "Habitat", "Group"]
        heatmap_cols = [c for c in metab_df.columns if c not in meta_cols]
        
        for col in metab_df.columns:
            
            if col not in meta_cols:
                
                metab_df[col] = metab_df[col].astype(float)

        sorted_gens = make_final_genome_order(metab_df)

        # Filter genomes missing in the dataframe and print warning
        filtered_sorted_gens = [g for g in sorted_gens if g in metab_df.index]
        missing = set(sorted_gens) - set(filtered_sorted_gens)
        
        if missing:
            print(f"Warning: these genomes are missing in data and will be skipped: {missing}")

        metab_df = metab_df.loc[filtered_sorted_gens]

        # Assign row colors based on Group
        row_color_mask = metab_df["Group"].replace({"SP": "#ff9999", "FL": "#66b3ff"})
        
        Habitat_color_dict={"Sponge": "#FBA276", "Free living": "#66b3ff", 
                            "Seagrass" : "#027EB2", "Foraminifera" : "#DDE130", 
                            "Coral" : "#E39CB2", "Algae" : "#50C647",
                            "Tunicata" : "#DC868D"}
        row_color_mask = metab_df["Habitat"].replace(Habitat_color_dict)
        
        # Assign colors based on taxonomy group.
        Tax_group_color_dict={"OTU1" : "#EB665F", "OTU4" : "#ECA42D", "OTU23" : "#EBE384",
                              "OTU3" : "#B1AA67", "OTU7" : "#D8AC7C", "OTU9" : "#CDE6DF",
                              "OTU14" : "#344541"}
        metab_df["Tax group color"] = metab_df["Relative"].replace(Tax_group_color_dict)

        # Debug information to check shapes and alignments
        print(f"File: {filename}")

        base_title = filename.replace(".xlsx", "")
        # Draw heatmaps for completeness data
        options=[False, [30,10], 5, 8]
        draw_heatmap(metab_df[heatmap_cols], row_color_mask, base_title, options, output_dir)
        
        # Perform PCA on normalized data
        draw_PCA_all(base_title, metab_df, metab_df[heatmap_cols], Sponge_assoc_MAGs, Habitat_color_dict, output_dir)  
        draw_PCA_groupwise(base_title, metab_df, metab_df[heatmap_cols], OTU_main_genomes, OTU_order, Habitat_color_dict, output_dir)
        
        # Compare mean and median completeness of modules in genomes.
        metab_df['Median module completeness']=metab_df[heatmap_cols].median(axis=1)
        metab_df['Mean module completeness']=metab_df[heatmap_cols].mean(axis=1)
        metab_df['Genome completeness']=checkm_data['Completeness']
        metab_df['Genome contamination']=checkm_data['Contamination']
        
        for OTU_name in OTU_order:
            
            # Median completeness.
            plot_module_compl_distrib(metab_df, base_title, OTU_name, 'Median', output_dir)
            
            # Mean completeness.
            plot_module_compl_distrib(metab_df, base_title, OTU_name, 'Mean', output_dir)      
            
        
        # Identify modules unique for SAB MAGs.
        for OTU_name in OTU_order:
            
            OTU_group_df=metab_df[metab_df['Relative']==OTU_name]
            print(f'{OTU_name} modules with completeness higher than in other genomes.')
            print(OTU_group_df.loc[:, (OTU_group_df.loc[OTU_main_genomes[OTU_name]] > OTU_group_df.loc[OTU_group_df.index != OTU_main_genomes[OTU_name]].max()) & (OTU_group_df.columns.isin(heatmap_cols))].to_string())
            
            print(f'{OTU_name} modules with completeness higher than in other free-living genomes.')                  
            print(OTU_group_df.loc[:, OTU_group_df.loc[OTU_main_genomes[OTU_name]] > OTU_group_df.loc[OTU_group_df['Group']=='FL'].max()].to_string())
               
        
        
        
        
        
        ## General metabolism pathways.
        #General_metabolis_mod=["M00001", "M00002", "M00003", "M00307", "M00009", "M00010", "M00011", "M00004", "M00006", "M00007", 
        #                       "M00005", "M00008", "M00308", "M00012", "M00165", "M00173", "M00376", "M00375", "M00374", "M00377",
        #                       "M00144", "M00149", "M00151", "M00155", "M00157"]
        ## Sulfur metabolism pathways.
        #Sulfur_metabolism_mod=["M00596", "M00595", "M00176"]  
        ## Amino acid biosynthetic pathway.
        #Amino_acid_biosyn_mod=["M00020", "M00018", "M00021", "M00609", "M00015", "M00017", "M00535", "M00570", "M00019", "M00432", 
        #                       "M00016", "M00525", "M00526", "M00527", "M00844", "M00845", "M00026", "M00023", "M00024", "M00025", 
        #                       "M00040"]
        ## Vitamin biosynthetic pathway.
        #Vitamin_biosynthe_mod=["M00125", "M00911", "M00119", "M00913", "M00120", "M00572", "M00123", "M00950", "M00573", "M00577", 
        #                       "M00126", "M00840", "M00841", "M00924", "M00122", "M00117", "M00846", "M00868", "M00121", 
        #                       "M00926", "M00124", "M00115", "M00912", "M00881", "M00882", "M00883", "M00884", "M00842", "M00880", 
        #                       "M00127", "M00895", "M00896", "M00897", "M00116"]
        ## Taurine and sulfoacetate catabolic pathways.
        #Taurine__sulfoacetate=["Taurine", "Sulfoacetate"]
        ## All analyzed features.
        #Long_list_of_features=General_metabolis_mod+Sulfur_metabolism_mod+Amino_acid_biosyn_mod+Vitamin_biosynthe_mod+Taurine__sulfoacetate
        ## All symbiotic features.
        #Symbiotic_met_featurs=["M00018", "M00017", "M00019", "M00535", "M00570", "M00432", "M00016", "M00525", "M00526", "M00527", 
        #                       "M00844", "M00845", "M00026", "M00023", "M00024", "M00025", "M00040", "M00127", "M00895", "M00896", 
        #                       "M00897", "M00125", "M00911", "M00119", "M00913", "M00123", "M00950", "M00573", "M00577", "M00924", 
        #                       "M00122", "M00572", "M00117", "M00846", "M00116", "Taurine", "Sulfoacetate"]
        #
        #meta_cols_SF=["Genome", "Relative", "Habitat", "Group", "Tax group color"]
        #
        ## Run analysis for symbiotic metabolic features only.
        #Symbiotic_met_featurs+=meta_cols_SF
        #metab_SF_df = metab_df[Symbiotic_met_featurs]
        #heatmap_SF_cols = [c for c in metab_SF_df.columns if c not in meta_cols_SF]
        #
        #base_title_SF = f'{base_title}_SF'
        ## Draw heatmaps for completeness data
        #options=[False, [15,20], 5, 8]
        #draw_heatmap(metab_SF_df[heatmap_SF_cols], row_color_mask, base_title_SF, options, output_dir)
        #
        ## Perform PCA on normalized data
        #draw_PCA_all(base_title_SF, metab_SF_df, metab_SF_df[heatmap_SF_cols], Sponge_assoc_MAGs, Habitat_color_dict, output_dir)  
        #draw_PCA_groupwise(base_title_SF, metab_SF_df, metab_SF_df[heatmap_SF_cols], OTU_main_genomes, OTU_order, Habitat_color_dict, output_dir)
        #
        ## Run for all analyzed metabolic features.
        #Long_list_of_features+=meta_cols_SF
        #metab_AF_df = metab_df[Long_list_of_features]
        #heatmap_AF_cols = [c for c in metab_AF_df.columns if c not in meta_cols_SF]
        #
        #base_title_AF = f'{base_title}_AF'
        ## Draw heatmaps for completeness data
        #options=[False, [10,10], 5, 8]
        #draw_heatmap(metab_AF_df[heatmap_AF_cols], row_color_mask, base_title_AF, options, output_dir)
        #
        ## Perform PCA on normalized data
        #draw_PCA_all(base_title_AF, metab_AF_df, metab_AF_df[heatmap_AF_cols], Sponge_assoc_MAGs, Habitat_color_dict, output_dir)
        #draw_PCA_groupwise(base_title_AF, metab_AF_df, metab_AF_df[heatmap_AF_cols], OTU_main_genomes, OTU_order, Habitat_color_dict, output_dir)     
        
        
        

        ## Split data by group for statistical analyses
        #fl_df = elp_df[elp_df["Group"] == "FL"][heatmap_cols]
        #sp_df = elp_df[elp_df["Group"] == "SP"][heatmap_cols]
        #fl_norm_df = norm_counts.loc[fl_df.index]
        #sp_norm_df = norm_counts.loc[sp_df.index]
        #
        ## Calculate statistics and plot distributions (raw and normalized)
        #print(f'Comparing raw ELP counts for {base_title}, all genomes considered.')
        #calc_stat_draw_distribs(base_title, fl_df, sp_df, "raw", output_dir)
        #print(f'Comparing normalized per Mb ELP counts for {base_title}, all genomes considered.')
        #calc_stat_draw_distribs(base_title, fl_norm_df, sp_norm_df, "normalized", output_dir)
        #
        ## Analyse only for SAB MAGs with free-living relatives: OTU4, OTU23, OTU3, OTU7.
        #elp_with_fl_df = elp_df[elp_df["Relative"].isin(["OTU4", "OTU23", "OTU3", "OTU7"])]
        #SABMAG_rel_with_fl_df = elp_with_fl_df[~elp_with_fl_df.index.isin(["Ca. Halichondribacter symbioticus", "Ca. Yagmuria paniceus", "Ca. Ahtobacter symbioticus", "Ca. Vellamobacter salmiensis"])][heatmap_cols]
        #SABMAG_with_fl_df = elp_with_fl_df[elp_with_fl_df.index.isin(["Ca. Halichondribacter symbioticus", "Ca. Yagmuria paniceus", "Ca. Ahtobacter symbioticus", "Ca. Vellamobacter salmiensis"])][heatmap_cols]
        #SABMAG_rel_with_fl_norm_df = norm_counts.loc[SABMAG_rel_with_fl_df.index]
        #SABMAG_with_fl_norm_df = norm_counts.loc[SABMAG_with_fl_df.index] 
        #
        ## Calculate statistics and plot distributions (raw and normalized)
        #base_title+='_OTUs_with_FLR'
        #print(f'Comparing raw ELP counts for {base_title}, genomes, related to OTU4, OTU23, OTU3, OTU7, considered.')
        #calc_stat_draw_distribs(base_title, SABMAG_rel_with_fl_df, SABMAG_with_fl_df, "raw", output_dir)
        #print(f'Comparing normalized per Mb ELP counts for {base_title}, genomes, related to OTU4, OTU23, OTU3, OTU7, considered.')
        #calc_stat_draw_distribs(base_title, SABMAG_rel_with_fl_norm_df, SABMAG_with_fl_norm_df, "normalized", output_dir)        

    ## Plot CheckM statistics once for all data
    #plot_bins_stats(checkm_data, Sponge_assoc_MAGs, output_dir)
    
    return


# -------------------- Run --------------------

wrapper_metabolism_vis_multiple(Metabolism_folder_path, Metabolism_filenames, CheckM_bins_info, Output_path)
