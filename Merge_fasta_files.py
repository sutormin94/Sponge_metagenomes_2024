###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2024##
##Merge fasta files##

#Merge fasta files in a single file. Remove duplicated sequences.
###############################################

#######
#Packages to be imported.
#######

from Bio import SeqIO, Align, Seq
import os


def merge_fasta_files(dataset_name, fasta_path_dict, output_path):
    
    Total_seq_dict={}
    
    count_overlap=0
    
    for file_ID, filepath in fasta_path_dict.items():
        
        for handle in SeqIO.parse(filepath, 'fasta'):
            
            seq_name=str(handle.name)
            seq_seq=str(handle.seq)
            
            if seq_name in Total_seq_dict:
                count_overlap+=1
                continue
            
            else:
                Total_seq_dict[seq_name]=seq_seq
    
    print(f'{count_overlap} sequences were duplicated.')
    print(f'{len(Total_seq_dict)} unique sequences were collected.')
    
    
    outpath=os.path.join(output_path, f'{dataset_name}.fasta')
    fileout=open(outpath, 'w')
    
    for Seq_name, Seq_seq in Total_seq_dict.items():
        fileout.write(f'>{Seq_name}\n{Seq_seq}\n')
        
    fileout.close()
    
    return


# Dataset name.
Dataset_name="OTU_9_OTU_14_500_nt_100_NCBI_16S_100_SILVA_auto_clust"

# Dictionary with pathes to input fasta files to be merged.
Fasta_path_dict={'OTU_9' : "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\\16S_metagenomics\Symbiotic_16S\OTU_9_full_length\OTU_9_500_nt_100_NCBI_16S_100_SILVA_auto_clust.fasta",
                 'OTU_14' : "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\\16S_metagenomics\Symbiotic_16S\OTU_14_full_length\OTU_14_500_nt_100_NCBI_16S_100_SILVA_auto_clust.fasta",}

# Output path.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\\16S_metagenomics\Symbiotic_16S\OTU_9_and_14_full_length\\"

merge_fasta_files(Dataset_name, Fasta_path_dict, Output_path)