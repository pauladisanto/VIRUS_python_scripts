#python3 Generate_bigger_fasta_file.py

import os
from Bio import SeqIO

def concatenate_fasta_files(input_dir, output_file):
    # Create or open the output file in write mode
    with open(output_file, "w") as output_handle:
        sequence_ids = set()  # Set to track sequence IDs to avoid duplicates
        
        # Iterate through each file in the input directory
        for file_name in os.listdir(input_dir):
            # Check for .fasta and .fa file extensions
            if file_name.endswith((".fasta", ".fa")):
                file_path = os.path.join(input_dir, file_name)
                print(f"Processing file: {file_path}")
                
                # Read each fasta file
                with open(file_path, "r") as input_handle:
                    for record in SeqIO.parse(input_handle, "fasta"):
                        if record.id not in sequence_ids:
                            SeqIO.write(record, output_handle, "fasta")
                            sequence_ids.add(record.id)
                        else:
                            print(f"Duplicate sequence ID {record.id} found in file {file_path}, skipping.")
    
    print(f"All FASTA files have been concatenated into: {output_file}")

# Directory containing your smaller FASTA files
input_dir = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Sequences_Gensam_Ours"


# Path for the output combined FASTA file
output_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/combined.fasta"

# Concatenate all FASTA files
concatenate_fasta_files(input_dir, output_file)


