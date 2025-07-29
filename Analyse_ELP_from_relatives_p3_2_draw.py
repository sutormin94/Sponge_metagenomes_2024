###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2025##
##Analyse abundance of ELP-domains and ELP-encoding genes from bacterial MAGs:
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
ELP_folder_path = "C:\\Users\\sutor\\OneDrive\\ThinkPad_working\\Sutor\\Science\\Spongy\\Scripts\\Sponge_metagenomes_2024\\Source_data\\ELPs_in_SAB_MAGS_and_related_species\\"  # Set this to the folder containing your 4 Excel files

# The 4 Excel files to process:
ELP_filenames = [
    "All_ELPs_types_domains_counts.xlsx",
    "All_ELPs_types_genes_counts.xlsx",
    "Exported_ELPs_types_domains_counts.xlsx",
    "Exported_ELPs_types_genes_counts.xlsx"
]

CheckM_bins_info = "C:\\Users\\sutor\\OneDrive\\ThinkPad_working\\Sutor\\Science\\Spongy\\Scripts\\Sponge_metagenomes_2024\\Source_data\\ELPs_in_SAB_MAGS_and_related_species\\MAGs_and_bins_CheckM.xlsx"
checkm_dir = os.path.dirname(os.path.abspath(CheckM_bins_info))
Output_path = os.path.join(checkm_dir, "ELP_stat_and_pictures")
os.makedirs(Output_path, exist_ok=True)

if not os.path.exists(Output_path):
    os.makedirs(Output_path)

# OTU main genomes and order (as before)
OTU_main_genomes = {
    "OTU4": "Ca. Halichondribacter symbioticus",
    "OTU23": "Ca. Yagmuria paniceus",
    "OTU3": "Ca. Ahtobacter symbioticus",
    "OTU7": "Ca. Vellamobacter salmiensis",
    "OTU9": "Ca. Vienanmeria sitiensis",
    "OTU14": "Ca. Sampovibrio pertsovi",
    "OTU1": "Ca. Eurynomebacter symbioticus"
}
OTU_order = ["OTU4", "OTU23", "OTU3", "OTU7", "OTU9", "OTU14", "OTU1"]
Sponge_assoc_MAGs = list(OTU_main_genomes.values())


# -------------------- Functions --------------------

def load_and_prepare_ELP_table(filepath):
    """
    Load one ELP data table and add 'Group' column.
    """
    df = pd.read_excel(filepath)
    df["Group"] = np.where(df["Habitat"].str.strip() == "Free living", "FL", "SP")
    return df


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


def norm_ELP_counts(elp_df_sub, checkm_df, heatmap_cols):
    joined = elp_df_sub.join(checkm_df[["Genome_Size"]], how='left')
    factor = joined["Genome_Size"] / 1e6
    normed = joined[heatmap_cols].div(factor, axis=0)
    return normed


import seaborn as sb
import matplotlib.pyplot as plt
import os

def draw_heatmap(input_df, row_colors, title, norm, output_path):
    """
    Plot heatmap with seaborn clustermap, with colored rows by group.
    Clustering disabled to preserve ordering.
    """
    cg = sb.clustermap(
        input_df,
        row_colors=row_colors,
        cmap='YlGnBu',
        linewidths=1.0,
        row_cluster=False,
        col_cluster=False,
        yticklabels=True,
        xticklabels=True,
    )  
    cg.fig.set_size_inches(10, 20)
    
    plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=11)
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), fontsize=9)
    cg.ax_heatmap.set_xlabel("ELP classes", fontsize=16)
    cg.ax_heatmap.set_ylabel("Genomes", fontsize=16)
    cg.fig.suptitle(title, fontsize=18)

    safe_title = title.replace(" ", "_")
    #plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{safe_title}_heatmap_{norm}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{safe_title}_heatmap_{norm}.svg'), dpi=300)
    plt.close()
    return


def draw_PCA_all(ELP_cat, all_data_df, norm, sponge_assoc_genomes, output_path):
    pca = sklearnPCA(n_components=None)
    pcs = pca.fit_transform(all_data_df)

    genome_names = all_data_df.index.tolist()
    labels, colors, alphas = [], [], []
    for genome in genome_names:
        if genome in sponge_assoc_genomes:
            labels.append(genome)
            colors.append('black')
            alphas.append(1)
        else:
            labels.append('')
            colors.append('blue')
            alphas.append(0.5)

    plt.figure(figsize=(6, 6))
    ax = plt.gca()
    ax.scatter(pcs[:, 0], pcs[:, 1], c=colors, alpha=alphas, s=30, edgecolors='k', linewidth=0.2)
    for i, label in enumerate(labels):
        if label:
            ax.text(pcs[i, 0], pcs[i, 1], label, fontsize=6)
    ax.set_title(f"PCA of {ELP_cat} ({norm})", fontsize=14)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=12)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_PCA_{norm}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_PCA_{norm}.svg'), dpi=300)
    plt.close()
    return

def draw_PCA_groupwise(ELP_cat, all_data_df, all_data_df_norm, norm, otu_main_genomes, otu_order, habitat_color_dict, output_path):
    
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
        
        model=pca(n_components=3) 
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
    
    
    #fig.set_size_inches(3, 3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_MAGs_relative_wise_biplot_{norm}.png'), dpi=300)
    plt.savefig(os.path.join(output_path, f'{ELP_cat}_MAGs_relative_wise_biplot_{norm}.svg'), dpi=300)  
    
    
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

def wrapper_ELP_vis_multiple(ELP_folder_path, ELP_filenames, checkm_file, output_dir):
    """
    Main wrapper function that processes multiple ELP Excel files in a folder.
    It loads each file, assigns groups based on Habitat, orders genomes by OTU and Relative,
    normalizes counts by genome size (CheckM data), performs PCA, statistical tests, and draws heatmaps.
    Output plots and results are saved to a directory which is created if missing.

    Parameters:
    -----------
    ELP_folder_path : str
        Path to the folder containing the 4 ELP Excel files.
    ELP_filenames : list of str
        List of the four Excel filenames to process.
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

    for filename in ELP_filenames:
        print(f"Processing file: {filename}")
        fullpath = os.path.join(ELP_folder_path, filename)
        elp_df = load_and_prepare_ELP_table(fullpath)

        meta_cols = ["Genome", "Relative", "Habitat", "Group"]
        heatmap_cols = [c for c in elp_df.columns if c not in meta_cols]

        sorted_gens = make_final_genome_order(elp_df)

        elp_df.set_index("Genome", inplace=True)

        # Filter genomes missing in the dataframe and print warning
        filtered_sorted_gens = [g for g in sorted_gens if g in elp_df.index]
        missing = set(sorted_gens) - set(filtered_sorted_gens)
        
        if missing:
            print(f"Warning: these genomes are missing in data and will be skipped: {missing}")

        elp_df = elp_df.loc[filtered_sorted_gens]

        # Normalize counts by genome size
        norm_counts = norm_ELP_counts(elp_df[heatmap_cols], checkm_data, heatmap_cols)
        norm_counts = norm_counts.loc[filtered_sorted_gens]

        # Assign row colors based on Group
        row_color_mask = elp_df["Group"].replace({"SP": "#ff9999", "FL": "#66b3ff"})
        
        Habitat_color_dict={"Sponge": "#FBA276", "Free living": "#66b3ff", 
                            "Seagrass" : "#027EB2", "Foraminifera" : "#DDE130", 
                            "Coral" : "#E39CB2", "Algae" : "#50C647",
                            "Tunicata" : "#DC868D"}
        row_color_mask = elp_df["Habitat"].replace(Habitat_color_dict)

        # Debug information to check shapes and alignments
        print(f"File: {filename}")

        base_title = filename.replace(".xlsx", "").replace("_", " ")
        # Draw heatmaps for normalized and raw data
        draw_heatmap(norm_counts, row_color_mask, base_title + " (normalized)", "normalized", output_dir)
        draw_heatmap(elp_df[heatmap_cols], row_color_mask, base_title + " (raw)", "raw", output_dir)
        
        # Perform PCA on normalized data
        #draw_PCA_all(base_title, norm_counts, "normalized", Sponge_assoc_MAGs, output_dir)  
        draw_PCA_groupwise(base_title, elp_df, norm_counts, "normalized", OTU_main_genomes, OTU_order, Habitat_color_dict, output_dir)

        # Split data by group for statistical analyses
        fl_df = elp_df[elp_df["Group"] == "FL"][heatmap_cols]
        sp_df = elp_df[elp_df["Group"] == "SP"][heatmap_cols]
        fl_norm_df = norm_counts.loc[fl_df.index]
        sp_norm_df = norm_counts.loc[sp_df.index]

        # Calculate statistics and plot distributions (raw and normalized)
        print(f'Comparing raw ELP counts for {base_title}, all genomes considered.')
        calc_stat_draw_distribs(base_title, fl_df, sp_df, "raw", output_dir)
        print(f'Comparing normalized per Mb ELP counts for {base_title}, all genomes considered.')
        calc_stat_draw_distribs(base_title, fl_norm_df, sp_norm_df, "normalized", output_dir)
        
        # Analyse only for SAB MAGs with free-living relatives: OTU4, OTU23, OTU3, OTU7.
        elp_with_fl_df = elp_df[elp_df["Relative"].isin(["OTU4", "OTU23", "OTU3", "OTU7"])]
        SABMAG_rel_with_fl_df = elp_with_fl_df[~elp_with_fl_df.index.isin(["Ca. Halichondribacter symbioticus", "Ca. Yagmuria paniceus", "Ca. Ahtobacter symbioticus", "Ca. Vellamobacter salmiensis"])][heatmap_cols]
        SABMAG_with_fl_df = elp_with_fl_df[elp_with_fl_df.index.isin(["Ca. Halichondribacter symbioticus", "Ca. Yagmuria paniceus", "Ca. Ahtobacter symbioticus", "Ca. Vellamobacter salmiensis"])][heatmap_cols]
        SABMAG_rel_with_fl_norm_df = norm_counts.loc[SABMAG_rel_with_fl_df.index]
        SABMAG_with_fl_norm_df = norm_counts.loc[SABMAG_with_fl_df.index] 
        
        # Calculate statistics and plot distributions (raw and normalized)
        base_title+='_OTUs_with_FLR'
        print(f'Comparing raw ELP counts for {base_title}, genomes, related to OTU4, OTU23, OTU3, OTU7, considered.')
        calc_stat_draw_distribs(base_title, SABMAG_rel_with_fl_df, SABMAG_with_fl_df, "raw", output_dir)
        print(f'Comparing normalized per Mb ELP counts for {base_title}, genomes, related to OTU4, OTU23, OTU3, OTU7, considered.')
        calc_stat_draw_distribs(base_title, SABMAG_rel_with_fl_norm_df, SABMAG_with_fl_norm_df, "normalized", output_dir)        

    # Plot CheckM statistics once for all data
    plot_bins_stats(checkm_data, Sponge_assoc_MAGs, output_dir)
    return


# -------------------- Run --------------------

wrapper_ELP_vis_multiple(ELP_folder_path, ELP_filenames, CheckM_bins_info, Output_path)
