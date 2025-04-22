import os
import subprocess
import sys

# Set paths

index_path = "/mnt/c/users/thorb/Documents/OfvirkiFiskurinn/RNASequencinganalysis/020_vik/kallisto_index/index.idx"
output_dir ="/mnt/c/users/thorb/Documents/OfvirkiFiskurinn/RNASequencinganalysis/020_vik/output_alt/"
fastq_base_dir = "/mnt/c/users/thorb/Documents/OfvirkiFiskurinn/RNASequencinganalysis/020_vik/all_fastq_thor_boe/"
kallisto_exec = "/mnt/c/users/thorb/Documents/OfvirkiFiskurinn/RNASequencinganalysis/020_vik/kallisto/kallisto"

# Loop through the 8 folders

mysamples = os.listdir(fastq_base_dir)
mysamples.sort()
for folder in mysamples:
    folder_path = os.path.join(fastq_base_dir, folder)
    
    if os.path.isdir(folder_path):
        print(f"Processing folder: {folder_path}")
        
        # Get all .fq files in the folder
        all_files = os.listdir(folder_path)
        all_files_cleaned = []
        for element in all_files:
            if element[0] != ".":
                print(element)
                full_path = folder_path+"/"+element
                all_files_cleaned.append(full_path)


        all_files_cleaned.sort()
        print(all_files_cleaned)
        all_files_cleaned_string = " ".join(all_files_cleaned)
        print(all_files_cleaned_string)

        # define output dir
        output_folder = output_dir+folder
        
        # execute kallisto
        cmd = [kallisto_exec, "quant", "-i", index_path, "-o", output_folder, "--threads=8", all_files_cleaned_string]
        cmd = " ".join(cmd)
        print()
        print(cmd)
        print()
        os.system(cmd)