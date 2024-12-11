###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2024##
##Prepare dataset for multiple alignment##

#Takes an input multi-fasta file with sequences collected from different sources: a querry, NCBI nt, NCBI rRNA, SILVA.
#Tests the orientation of sequences and makes reverse-complement if needed.
#Runs mmseqs2 with properly oriented sequences with 99.5% identity level to make a non-redundant dataset.
#Reads mmseqs output:
#For a cluster chooses a representative record from an NCBI database (if avaliable).
#Concatenates NCBI record name with shortened SILVA name (if avaliable).
#Adds information about host and isolation source from full genebank records for NCBI nt and NCBI rRNA databases.
#Returns a new multi-fasta file ready for multiple alignment.
###############################################

#######
#Packages to be imported.
#######

from Bio import SeqIO, Align, Seq
import os
from os import listdir, mkdir
from Bio.Seq import Seq
import time
import subprocess
import copy

#################
### Variables to be defined.
#################

#Path to working directory.
PWD_base="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\\16S_metagenomics\Symbiotic_16S\OTU_3_full_length\\"

#Dataset name.
Dataset_name="OTU_3_500_nt_100_NCBI_16S_100_SILVA_plus_GtDB_auto"

#Query sequence.
Query_seq_path=PWD_base + "OTU_3.fasta"

#nt BLAST results.
NT_blast_res_path=PWD_base + "OTU_3_500_nt.fasta"

#nt BLAST results, full Genbank files.
NT_blast_res_gb_path=PWD_base + "OTU_3_500_nt_full_genbank.gb"

#NCBI 16S rRNA BLAST results.
NCBI_16S_blast_res_path=PWD_base + "OTU_3_100_rRNA_NCBI.fasta"

#NCBI 16S rRNA BLAST results, full Genbank files.
NCBI_16S_blast_res_gb_path=PWD_base + "OTU_3_100_rRNA_NCBI_full_genbank.gb"

#SILVA BLAST results.
SILVA_blast_res_path=PWD_base + "OTU_3_200_SILVA_nr99.fasta"



#################
### Read input data.
#################

########
# Read input fasta file.
########

def read_multi_fasta(inpath, seq_dict):
    
    input_fasta=SeqIO.parse(inpath, "fasta")
        
    for record in input_fasta:
        rec_full_name=record.description
        rec_seq=record.seq
        if rec_full_name not in seq_dict:
            seq_dict[rec_full_name]=rec_seq
        else:
            print(f'The sequence "{rec_full_name}" was already present, not included now.')
    
    return seq_dict


########
# Test sequence orientation.
# Returns reverse compelement if a sequence is a reverse complement.
########

def test_sequence_orientation(seq_dict):
    
    aligner=Align.PairwiseAligner()
    aligner.mode='local'
    aligner.match_score=2
    aligner.mismatch_score=-3
    aligner.open_gap_score=-7
    aligner.extend_gap_score=-2    
    
    Test_seq="GGAGGCAGCAGTGGGGAATCT"
    
    rc_count=0
    for record_name, record_seq in seq_dict.items():
        
        record_seq=Seq(str(record_seq).replace('-', '')) #Remove gap symbols.
        alignment=aligner.align(record_seq, Test_seq)
        
        if len(alignment)>0 and alignment[0].score>22:
            seq_dict[record_name]=record_seq
        else:
            record_seq_rc=record_seq.reverse_complement()
            alignment_rc=aligner.align(record_seq_rc, Test_seq)  
            
            if len(alignment_rc)>0 and alignment_rc[0].score>22:
                seq_dict[record_name]=record_seq_rc
                
                rc_count+=1
            else:
                print(f'Warning! Record {record_name} was not aligned with an rRNA-specific sequence in any orientation.')
    
    print(f'A number of reverse-complement sequences detected: {rc_count}')
    
    return seq_dict


########
# Read GenBank file.
########

def read_GenBank(gb_inpath, gb_dict):
    
    input_genbank=SeqIO.parse(gb_inpath, "genbank")
        
    for record in input_genbank: 
        record_info=''
        
        record_source_info=record.features[0].qualifiers
        
        if 'organism' in record_source_info:
            record_info+=record_source_info['organism'][0].replace(' ', '_')
        if 'isolation_source' in record_source_info:
            record_info+='##'
            record_info+=record_source_info['isolation_source'][0].replace(' ', '_')
        if 'host' in record_source_info:
            record_info+='##'
            record_info+=record_source_info['host'][0].replace(' ', '_')
        if 'type_material' in record_source_info:
            record_info+='##'
            record_info+=record_source_info['type_material'][0].replace(' ', '_')
        
        record_end=record.features[0].location.end
        record_start=record.features[0].location.start
        record_length=int(record_end)-int(record_start)+1
        
        if record_length>1000000:
            record_info+='//'
            record_info+=f'{record.id}:{record_length}'
        
        print(record_info) 
        
        if record.id not in gb_dict:
            gb_dict[record.id]=record_info
    
    return gb_dict


########
# Write fasta file.
########

def write_fasta(seq_dict, outfile_path):
    
    outfile=open(outfile_path, 'w')
    
    for record_name, record_seq in seq_dict.items():
        record_name=record_name.replace(' ', '_')
        outfile.write(f'>{record_name}\n{str(record_seq)}\n')
        
    outfile.close()
    
    return


########
# Add GenBank info to fasta records.
########

def add_GenBank_info(seq_dict, gb_dict):
    
    seq_dict_upd={}
    
    for fasta_record_name, fasta_seq in seq_dict.items():
        check_add=0
        
        for gb_id, gb_info in gb_dict.items():
            if gb_id in fasta_record_name:
                fasta_record_name_upd=f'{fasta_record_name}##{gb_info}'
                seq_dict_upd[fasta_record_name_upd]=fasta_seq
                check_add+=1
        
        if check_add==0:
            seq_dict_upd[fasta_record_name]=fasta_seq
                
    return seq_dict_upd


########
# Run mmseqs2.
########

def run_mmseqs2(fasta_outpath, dataset_name, pwd_base): #For linux only.
    
    #Create temporary dir.
    command_line_mkdir_tmp=['mkdir', f'{pwd_base}/tmp_dir']
    subprocess.run(command_line_mkdir_tmp) 
    
    #Create output dir.
    command_line_mkdir_res=['mkdir', f'{pwd_base}/mmseq']
    subprocess.run(command_line_mkdir_res)    
    
    #Converting FASTA file to MMSEQS database
    command_line_createdb=['/usr/bin/mmseqs', 'createdb', fasta_outpath, f'{dataset_name}_DB']
    subprocess.run(command_line_createdb) 
    
    #Then execute the clustering.
    command_line_clustering=['/usr/bin/mmseqs', 'cluster', f'{dataset_name}_DB', f'{dataset_name}_clu', 
                             f'{pwd_base}/tmp_dir', '--min-seq-id', '0.999', '-c', '0.999',
                             '--single-step-clustering', '0', '--cluster-mode', '1', '--remove-tmp-files', '1',
                             '-v', '1', '--split-memory-limit', '1000000']
    subprocess.run(command_line_clustering) 
    
    #Generate a TSV-style formatted output file from the ffindex output file.
    command_line_createtsv=['/usr/bin/mmseqs', 'createtsv', f'{dataset_name}_DB', f'{dataset_name}_DB', 
                            f'{dataset_name}_clu', f'{dataset_name}_clu.tsv'] 
    subprocess.run(command_line_createtsv) 

    return


########
# Read mmseqs2 output.
########

def read_mmseqs_output(fasta_seq_dict, mmseqs_data_path):
    
    #Read mmseqs data.
    mmseqs_infile=open(mmseqs_data_path, 'r')
    
    Cluster_dict={}
    
    for line in mmseqs_infile:
        line=line.rstrip().split('\t')
        main_name=line[0]
        dep_name=line[1]
        if main_name not in Cluster_dict:
            Cluster_dict[main_name]=[main_name]
        else:
            Cluster_dict[main_name].append(dep_name)
    
    mmseqs_infile.close()
    
    print(f'Number of clusters detected: {len(Cluster_dict)}')
    
    #Polish cluster content.
    Full_names_ar=[]
    
    for cluster_name, cluster_content in Cluster_dict.items():
        base_name_ar=[]
        for record_name in cluster_content:
            record_name_base=record_name.split('.')[0]
            base_name_ar.append(record_name_base)
        base_name_ar=list(set(base_name_ar))
        
        full_name_ar=[]
        
        for record_name_base in base_name_ar:
            combined_full_name=''
            for record_name in cluster_content:
                if record_name_base in record_name:
                    if combined_full_name=='':
                        combined_full_name+=record_name
                    else:   
                        if '##' in record_name:
                            record_name_pars_ar=record_name.split('##')[1:]
                            for rec_add_info in record_name_pars_ar:
                                if rec_add_info not in combined_full_name:
                                    combined_full_name+=f'##{rec_add_info}'
            
            full_name_ar.append(combined_full_name)
            
        cluster_combined_name=''
        for comb_name in full_name_ar:
            if cluster_combined_name=='':
                cluster_combined_name+=comb_name  
            else:
                if '##' in comb_name:
                    record_name_pars_ar=comb_name.split('##')[1:]
                    for rec_add_info in record_name_pars_ar:
                        if rec_add_info not in cluster_combined_name:
                            cluster_combined_name+=f'##{rec_add_info}'
                            
        #Organize full names.
        if '//' in cluster_combined_name:
            rec_genome_info_ar=cluster_combined_name.split('//')
            cluster_combined_name_base=rec_genome_info_ar[0]
            cluster_combined_name_ar=rec_genome_info_ar[1:]
            for i in range(len(cluster_combined_name_ar)):
                if '##' in cluster_combined_name_ar[i]:
                    info_ele_ar=cluster_combined_name_ar[i].split('##')
                    cluster_combined_name_ar[i]=info_ele_ar
            
            Genome_info_str=""
            Additn_info_str=""
            if len(cluster_combined_name_ar)==1:
                if isinstance(cluster_combined_name_ar[0], str)==True:
                    Genome_info_str+=f'//{cluster_combined_name_ar[0]}'
                else:
                    Genome_info_str+=f'//{cluster_combined_name_ar[0][0]}'
                    for j in range(len(cluster_combined_name_ar[0])-1):
                        Additn_info_str+=f'##{cluster_combined_name_ar[0][j+1]}'
            else:
                for ele in cluster_combined_name_ar:
                    if isinstance(ele, str)==True:
                        Genome_info_str+=f'//{ele}'
                    else:
                        Genome_info_str+=f'//{ele[0]}'
                        for j in range(len(ele)-1):
                            Additn_info_str+=f'##{ele[j+1]}'
            
            cluster_combined_name=cluster_combined_name_base+Additn_info_str+Genome_info_str  
        
        Full_names_ar.append(cluster_combined_name)
        
    Final_seq_dict={}
    
    for cluster_combined_name in Full_names_ar:
        record_name_base=cluster_combined_name.split('.')[0]
        take_first=0
        for record_id, record_seq in fasta_seq_dict.items():
            if record_name_base in record_id:
                if take_first==0:
                    Final_seq_dict[cluster_combined_name]=record_seq
                take_first=1
                
    #Take sequences with full genomes avaliable only.
    Finel_seq_dict_fg={}
    for full_name, record_seq in Final_seq_dict.items():
        if '//' in full_name:
            Finel_seq_dict_fg[full_name]=record_seq
           
    return Final_seq_dict, Finel_seq_dict_fg



########
# Wrapper function.
########

def wrapper_func(query_seq_path, nt_blast_res_path, ncbi_16s_blast_res_path, 
                 silva_blast_res_path, nt_blast_res_gb_path, ncbi_16s_blast_res_gb_path, 
                 pwd_base, dataset_name):
    
    #Initialize a timer.
    starttime=time.time()
    
    Seq_dict={} #A dictionary to keep fasta sequences.
    
    #Read input query.
    Seq_dict=read_multi_fasta(query_seq_path, Seq_dict)
    print(f'Reading a query sequence. Working time since start, min: {(time.time()-starttime)/60}')
    print(f'Total number of sequences: {len(Seq_dict)}')
    
    Seq_dict=read_multi_fasta(nt_blast_res_path, Seq_dict)
    print(f'Reading nt blast results. Working time since start, min: {(time.time()-starttime)/60}')
    print(f'Total number of sequences: {len(Seq_dict)}')
    
    Seq_dict=read_multi_fasta(ncbi_16s_blast_res_path, Seq_dict)
    print(f'Reading NCBI 16S blast results. Working time since start, min: {(time.time()-starttime)/60}') 
    print(f'Total number of sequences: {len(Seq_dict)}')
    
    Seq_dict=read_multi_fasta(silva_blast_res_path, Seq_dict)
    print(f'Reading SILVA blast results. Working time since start, min: {(time.time()-starttime)/60}')  
    print(f'Total number of sequences: {len(Seq_dict)}')
    
    #Test sequence orientation and make reverse-complement if needed.
    Seq_dict=test_sequence_orientation(Seq_dict)
    print(f'Polishign sequence orientation. Working time since start, min: {(time.time()-starttime)/60}')  
    print(f'Total number of sequences: {len(Seq_dict)}')
    
    #Read NT blast results Genbank file.
    Gb_dict={}
    Gb_dict=read_GenBank(nt_blast_res_gb_path, Gb_dict)
    
    #Read NCBI 16S blast results Genbank file.
    Gb_dict=read_GenBank(ncbi_16s_blast_res_gb_path, Gb_dict)
    print(f'Parsing GenBank files. Working time since start, min: {(time.time()-starttime)/60}')
    print(f'Total number of GenBank records: {len(Gb_dict)}')
    
    #Add GenBank info to fasta records.
    Seq_dict=add_GenBank_info(Seq_dict, Gb_dict)
    
    #Write polished fasta set.
    fasta_outpath=pwd_base+dataset_name+'.fasta'
    write_fasta(Seq_dict, fasta_outpath)
    print(f'Writing sequences. Working time since start, min: {(time.time()-starttime)/60}') 
    
    #Run mmseqs2 to cluster identical sequences.
    run_mmseqs2(fasta_outpath, dataset_name, pwd_base)
    print(f'Clustering sequences with MMseqs2. Working time since start, min: {(time.time()-starttime)/60}') 
    
    #Read mmseqs2 output.
    Combined_clust_full_dict, Combined_clust_full_fg_dict=read_mmseqs_output(Seq_dict, pwd_base+dataset_name+'_clu.tsv')
    
    #Write clustered fasta set.
    fasta_outpath_clust=pwd_base+dataset_name+'_clust.fasta'
    write_fasta(Combined_clust_full_dict, fasta_outpath_clust)
    
    fasta_outpath_clust_fg=pwd_base+dataset_name+'_fg_clust.fasta'
    write_fasta(Combined_clust_full_fg_dict, fasta_outpath_clust_fg)
    
    return

#wrapper_func(Query_seq_path, NT_blast_res_path, NCBI_16S_blast_res_path, 
#             SILVA_blast_res_path, NT_blast_res_gb_path, NCBI_16S_blast_res_gb_path, 
#             PWD_base, Dataset_name)



#################
### Prepare data and metadata for iTOL.
#################

def read_final_fasta_and_reformat(fasta_final_inpath, fasta_final_outpath):
    
    input_fasta=SeqIO.parse(fasta_final_inpath, "fasta")
    
    fasta_iTOL_out=open(fasta_final_outpath+'.fasta', 'w')
    meta_iTOL_out=open(fasta_final_outpath+'.tsv', 'w')
    
    Records_dict={}
    
    meta_iTOL_out.write('ID\tBiome\tBiome comment\tFull genome\tFull genome comment\tRecord additional info\tSequence\n')
        
    for record in input_fasta:
        Record_info_dict={'ID' : '', 'Biome' : '', 'Biome comment' : '',
                          'Full genome' : '', 'Full genome comment' : '',
                          'Record additional info' : '', 'Sequence' : '',}
        
        record_full_name=str(record.id)
        record_seq=str(record.seq)
        Record_info_dict['Sequence']=record_seq
        
        if '//' in record_full_name:
            record_full_name_ar_gen=record_full_name.split('//')[1:]
            Record_info_dict['Full genome']=1
            gen_id_str=''
            for ele in record_full_name_ar_gen:
                gen_id_str+=f'//{ele}'
            print(gen_id_str)
            gen_id_str=gen_id_str.lstrip('//')
            Record_info_dict['Full genome comment']=gen_id_str
            
            record_full_name_nogen=record_full_name.split('//')[0]
            if '##' in record_full_name_nogen:
                record_full_name_nogen_ar_src=record_full_name_nogen.split('##')[1:]
                src_str=''
                for ele in record_full_name_nogen_ar_src:
                    src_str+=f'//{ele}'
                src_str=src_str.lstrip('//')
                Record_info_dict['Biome comment']=src_str  
                
                record_full_name_base=record_full_name_nogen.split('##')[0]
                record_full_name_base_ar=record_full_name_base.split(':')
                record_full_name_id=record_full_name_base_ar[0]
                record_full_name_comment=record_full_name_base_ar[1]
                Record_info_dict['ID']=record_full_name_id
                Record_info_dict['Record additional info']=record_full_name_comment
                
            else:
                Record_info_dict['Biome']=-1
                
                record_full_name_base_ar=record_full_name_nogen.split(':')
                record_full_name_id=record_full_name_base_ar[0]
                record_full_name_comment=record_full_name_base_ar[1]
                Record_info_dict['ID']=record_full_name_id
                Record_info_dict['Record additional info']=record_full_name_comment                
                
                
        else:
            Record_info_dict['Full genome']=-1
            
            if '##' in record_full_name:
                record_full_name_nogen_ar_src=record_full_name.split('##')[1:]
                src_str=''
                for ele in record_full_name_nogen_ar_src:
                    src_str+=f'//{ele}'
                src_str=src_str.lstrip('//')
                Record_info_dict['Biome comment']=src_str  
                
                record_full_name_base=record_full_name.split('##')[0]
                record_full_name_base_ar=record_full_name_base.split(':')
                record_full_name_id=record_full_name_base_ar[0]
                #print(record_full_name_base_ar)
                record_full_name_comment=record_full_name_base_ar[1]
                Record_info_dict['ID']=record_full_name_id
                Record_info_dict['Record additional info']=record_full_name_comment
                
            else:
                Record_info_dict['Biome']=-1 
                
                record_full_name_base_ar=record_full_name.split(':')
                record_full_name_id=record_full_name_base_ar[0]
                record_full_name_comment=record_full_name_base_ar[1]
                Record_info_dict['ID']=record_full_name_id
                Record_info_dict['Record additional info']=record_full_name_comment 
        
        fasta_iTOL_out.write(f'>{Record_info_dict["ID"]}\n{Record_info_dict["Sequence"]}\n')
        meta_iTOL_out.write(f'{Record_info_dict["ID"]}\t{Record_info_dict["Biome"]}\t{Record_info_dict["Biome comment"]}\t{Record_info_dict["Full genome"]}\t{Record_info_dict["Full genome comment"]}\t{Record_info_dict["Record additional info"]}\t{Record_info_dict["Sequence"]}\n')
    
    
    fasta_iTOL_out.close()
    meta_iTOL_out.close()
    
    return

#Fasta_final_inpath=PWD_base+Dataset_name+'_fg_clust_polyshed_for_iTOL.fasta'
#Fasta_final_outpath=PWD_base+Dataset_name+'_fg_clust_iTOL_fin'
#read_final_fasta_and_reformat(Fasta_final_inpath, Fasta_final_outpath)


#################
### Prepare metadata for iTOL from a string of leaves names.
#################


def read_leaves_names_make_metadata(leaves_names_path, metadata_csv_outpath):
    
    leaves_names_in=open(leaves_names_path, 'r')
    metadata_out=open(metadata_csv_outpath, 'w')
    
    for line in leaves_names_in:
        line=line.rstrip().split('\t')
        for leaf_name in line:
            try:
                float(leaf_name)
            except:
                leaf_name=leaf_name.replace('**', '\t')
                metadata_out.write(f'{leaf_name}\n')
    
    
    leaves_names_in.close()
    metadata_out.close()
    
    return

#Leaves_names_path=PWD_base+Dataset_name+'_names_OTU_1_leaves_with_GtDB_final_sub_tree.txt'
#Metadata_csv_outpath=PWD_base+Dataset_name+'_names_OTU_1_leaves_with_GtDB_final_sub_tree_metadata.csv'
#read_leaves_names_make_metadata(Leaves_names_path, Metadata_csv_outpath)



#################
### Guess biome type based on keywords from sample metadata.
### Choose coolor code based on biome name.
#################

def decipher_color_code_and_biome(name_str_descr):
    
    Biome_color_name_dict={"Seawater" : ["#86b8dc", ["seawater", "Seawater", "marine_water"]], 
                           "Marine sediments" : ["#b3b1b0", ["marine_sediment", "sea_sediment", "sediments_sea", "sediment_sea", "nodule"]], 
                           "Hydrotermal vent": ["#eb1f59", ["hydrotermal", "vent", "chimney"]], 
                           "Deep-sea cold seep carbonates" : ["#a6ccd3", ["seep", "seep_carbonates"]],
                           "Marine coral" : ["#e39cb2", ["coral"]], 
                           "Marine sponge" : ["#fba276", ["sponge"]], 
                           "Marine polychaete" : ["#ff8d00", ["polychaete"]], 
                           "Marine mollusc" : ["#dc868d", ["mollusc", "bivalve"]], 
                           "Vestimentiferan tubeworm" : ["#ff8d00", ["tubeworm", "vestimentifer", "trophosome", "pogonophor"]],
                           "Plant" : ["#3fe312", ["plant", "rice", "wheat", "parsley"]], 
                           "Mangrove sediments" : ["#a67a65", ["mangrove"]], 
                           "Seaweed" : ["#4da131", ["Seaweed"]],
                           "Marine fish" : ["#141384", ["marine_fish", "fish_marine", "sea_fish", "fish_sea"]], 
                           "Biofilm" : ["#f2de49", ["biofilm"]], 
                           "Oolitic sand" : ["#f7f1e3", ["Oolitic", "oolitic"]], 
                           "Oil spill": ["#f2faad", ["oil", "hydrocarbon", "petroleum"]], 
                           "Salt marsh" : ["#3b7524", ["marsh"]], 
                           "Saltern" : ["#ccc8a1", ["Saltern", "saltern"]], 
                           "Salt lake" : ["#5486e3", ["salt_lake", "lake_salt"]], 
                           "Soda lake" : ["#d97deb", ["soda_lake", "alkaline_lake", "lake_soda", "lake_alkaline"]],	
                           "Freshwater" : ["#50c0cf", ["pond", "river", "lake", "spring", "creek"]], 
                           "Freshwater sediments" : ["#b3b1b0", ["river_sediment", "pond_sediment", "lake_sediment"]], 
                           "Soil" : ["#c46935", ["soil"]], 
                           "Hot spring" : ["#e02e5a", ["hot_spring"]], 
                           "Activated sludge" : ["#b87272", ["activated", "sludge"]],}
    
    color_code=""
    biome_name=""
    
    for Biome_name, Biome_data in Biome_color_name_dict.items():
        
        for keyword in Biome_data[1]:
            
            if keyword in name_str_descr:
                
                if len(color_code)==0:
                    color_code+=Biome_data[0]
                    biome_name+=Biome_name
                else:
                    color_code+=f'.{Biome_data[0]}'
                    biome_name+=f'.{Biome_name}'            
    
    
    return color_code, biome_name


#################
### Reads newik format, extracts leaves names, makes leaves names shorter.
### For bootstrap consensus trees.
#################

def read_fix_nwk(nwk_tree_path, metadata_csv_outpath, nwk_tree_outpath):
    
    filein=open(nwk_tree_path, 'r')
    meta_fileout=open(metadata_csv_outpath, 'w')
    tree_fileout=open(nwk_tree_outpath, 'w')
    
    tree_ar=[]
    
    for line in filein:
        
        for i in range(len(line)):
            if (line[i]==",") and (line[i-1].isdigit()) and (line[i+1].isdigit()):
                line=line[:i] + '.' + line[i+1:]
                
        line_copy=copy.deepcopy(line)
        
        state="init"
        start=0
        end=0
        quotes_indicator=0
        
        for i in range(len(line)):
            if (line[i]=="'") and (state=="init"):
                start=i
                state="start_name"
                quotes_indicator=1
                continue
                
            elif (line[i]=="'") and (state=="start_name") and (quotes_indicator==1):
                end=i
                name_str=line[start:end+1]
                
                name_str_orig=copy.deepcopy(name_str)
                
                if ':' in name_str:
                    name_str=name_str.replace(':', '.')
                    
                name_str_short=name_str[1:-1].split('.')[0]
                
                if ' ' in name_str_short:
                    name_str_short=name_str_short.replace(' ', '_')
                    
                if ".1" not in name_str_short:
                    name_str_short=f'{name_str_short}.1'
                    
                name_str_descr=name_str[1:-1].split(':')[1].replace('**', '\t')
                
                print(name_str_short, name_str_descr)
                color_code, biome_name=decipher_color_code_and_biome(name_str_descr)
                meta_fileout.write(f'{name_str_short}\t{color_code}\t{biome_name}\t{name_str_descr}\n')
                
                line_copy=line_copy.replace(name_str_orig, name_str_short)
                
                state="init"
                quotes_indicator=0
                continue
                
            elif (line[i].isalpha()) and (state=="init"):
                start=i
                state="start_name"
                continue  
            
            elif (line[i] in [',', '(', ')']) and (state=="start_name") and (quotes_indicator==0):
                end=i
                name_str=line[start:end]
                
                name_str_orig=copy.deepcopy(name_str)
                
                if ':' not in name_str:
                    name_str=name_str.replace('.', ':', 1) 
                
                name_str_short=name_str.split(':')[0]
                
                if ' ' in name_str_short:
                    name_str_short=name_str_short.replace(' ', '_') 
                    
                if ".1" not in name_str_short:
                    name_str_short=f'{name_str_short}.1'
                    
                name_str_descr=name_str.split(':')[1].replace('**', '\t')
                
                print(name_str_short, name_str_descr)
                color_code, biome_name=decipher_color_code_and_biome(name_str_descr)
                meta_fileout.write(f'{name_str_short}\t{color_code}\t{biome_name}\t{name_str_descr}\n')
                
                line_copy=line_copy.replace(name_str_orig, name_str_short)
                
                state="init"
                continue    
            
        #print(line_copy)
        tree_fileout.write(f'{line_copy}')
                

    
    filein.close()
    meta_fileout.close()
    tree_fileout.close()  
    
    return


#Nwk_tree_path=PWD_base+Dataset_name+'_clust_MUSCLE_alignment_full_random_tree.nwk'
#Metadata_csv_outpath=PWD_base+Dataset_name+'_leaves_names_final_full_random_tree_metadata.csv'
#Nwk_tree_outpath=PWD_base+Dataset_name+'_clust_MUSCLE_alignment_full_random_tree_iTOL.nwk'
#read_fix_nwk(Nwk_tree_path, Metadata_csv_outpath, Nwk_tree_outpath)




#################
### Reads newik format, extracts leaves names, makes leaves names shorter.
### For random trees with branch length and consensus values.
#################

def read_fix_rand_nwk(nwk_tree_path, metadata_csv_outpath, nwk_tree_outpath):
    
    filein=open(nwk_tree_path, 'r')
    meta_fileout=open(metadata_csv_outpath, 'w')
    tree_fileout=open(nwk_tree_outpath, 'w')
    
    tree_ar=[]
    
    for line in filein:
        
        for i in range(len(line)):
            if (line[i]==",") and (line[i-1].isdigit()) and (line[i+1].isdigit()):
                line=line[:i] + '.' + line[i+1:]
                
        line_copy=copy.deepcopy(line)
        
        state="init"
        start=0
        end=0
        quotes_indicator=0
        
        for i in range(len(line)):
            if (line[i]=="'") and (state=="init"):
                start=i
                state="start_name"
                quotes_indicator=1
                continue
                
            elif (line[i]=="'") and (state=="start_name") and (quotes_indicator==1):
                end=i
                name_str=line[start:end+1]
                
                name_str_orig=copy.deepcopy(name_str)
                
                if ':' in name_str:
                    name_str=name_str.replace(':', '.')
                    
                name_str_short=name_str[1:-1].split('.')[0]
                
                if ' ' in name_str_short:
                    name_str_short=name_str_short.replace(' ', '_')
                    
                if ".1" not in name_str_short:
                    name_str_short=f'{name_str_short}.1'
                    
                name_str_descr=name_str[1:-1].replace('**', '\t')
                
                print(name_str_short, name_str_descr)
                color_code, biome_name=decipher_color_code_and_biome(name_str_descr)
                meta_fileout.write(f'{name_str_short}\t{color_code}\t{biome_name}\t{name_str_descr}\n')
                
                line_copy=line_copy.replace(name_str_orig, name_str_short)
                
                state="init"
                quotes_indicator=0
                continue
                
            elif (line[i].isalpha()) and (state=="init"):
                start=i
                state="start_name"
                continue  
            
            elif (line[i] in [',', '(', ')']) and (state=="start_name") and (quotes_indicator==0):
                end=i
                name_str=line[start:end]
                
                name_str_orig=copy.deepcopy(name_str)
                
                if ':' in name_str:
                    name_str=name_str.replace(':', '.') 
                
                name_str_short=name_str.split('.')[0]
                
                if ' ' in name_str_short:
                    name_str_short=name_str_short.replace(' ', '_') 
                    
                if ".1" not in name_str_short:
                    name_str_short=f'{name_str_short}.1'
                    
                name_str_descr=name_str.replace('**', '\t')
                
                print(name_str_short, name_str_descr)
                color_code, biome_name=decipher_color_code_and_biome(name_str_descr)
                meta_fileout.write(f'{name_str_short}\t{color_code}\t{biome_name}\t{name_str_descr}\n')
                
                line_copy=line_copy.replace(name_str_orig, name_str_short)
                
                state="init"
                continue    
            
        #print(line_copy)
        tree_fileout.write(f'{line_copy}')
                

    
    filein.close()
    meta_fileout.close()
    tree_fileout.close()  
    
    return


Nwk_tree_path=PWD_base+Dataset_name+'_clust_MUSCLE_alignment_full_random_tree.nwk'
Metadata_csv_outpath=PWD_base+Dataset_name+'_leaves_names_final_full_random_tree_metadata.csv'
Nwk_tree_outpath=PWD_base+Dataset_name+'_clust_MUSCLE_alignment_full_random_tree_iTOL.nwk'
read_fix_rand_nwk(Nwk_tree_path, Metadata_csv_outpath, Nwk_tree_outpath)