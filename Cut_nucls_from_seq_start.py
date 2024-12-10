#!/usr/bin/env python3

import argparse
import os
import csv
import subprocess

parser = argparse.ArgumentParser(description="Different HEADCROP with Trimmomatic")
parser.add_argument("-f", "--file", help="File with sample names and seqs")
parser.add_argument("-i", "--initial", help="Initial directory")
parser.add_argument("-o", "--output", help="Output dictionary")
args, unknown = parser.parse_known_args()

#Creating dictionary with sample names
names_dict = {}

#Reading the input directory content.
Input_dir_path=args.initial
Input_dir_file_list=os.listdir(Input_dir_path)

#Creating the output directory.
Output_dir_path=args.output
if os.path.exists(Output_dir_path)==False:
    os.mkdir(Output_dir_path)

#Moving files.
with open(args.file,"r") as file:
    for line in csv.reader(file, delimiter='\t'):
        names_dict[line[0]] = line[1]
    
    for names, nucls in names_dict.items():    
        file_name_prefix=names
        file_cut_nucls=nucls
        file_name=f'{file_name_prefix}_1.fastq.gz'
        file_new_name=f'{file_name_prefix}_1.fastq.gz'

        
        if file_name in Input_dir_file_list:
            print(f'File {file_name} is found!')
            if file_cut_nucls != 'No':
                cut_nucl=len(file_cut_nucls) + 17
                print(f'At the beginning of the sequence {file_name}, {cut_nucl} nucleotides need to be cut off')
                command_line = ['java', '-jar', '/home/niagara/Progs/Trimmomatic-0.39/trimmomatic-0.39.jar', 'SE', '-threads', '20', '-phred33',
                          os.path.join(Input_dir_path, file_name), os.path.join(Output_dir_path, file_new_name),
                          f'HEADCROP:{cut_nucl}']
                print(command_line)
                subprocess.run(command_line)
            else:
                command_line_no = ['cp', os.path.join(Input_dir_path, file_name), os.path.join(Output_dir_path, file_new_name)]
                subprocess.run(command_line_no)
        else:
            print(f'File {file_name} not found!')
