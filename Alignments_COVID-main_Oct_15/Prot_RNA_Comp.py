#conda activate python_course 

from Bio import SeqIO
from Bio.Seq import Seq

# Function to process the FASTA file and save DNA complement, RNA, and protein sequences to separate files
def process_fasta(fasta_file, complement_file, rna_file, protein_file):
    with open(complement_file, 'w') as comp_f, open(rna_file, 'w') as rna_f, open(protein_file, 'w') as prot_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            dna_seq = record.seq
            
            # Get DNA complement
            dna_complement = dna_seq.complement()
            
            # Transcribe DNA to RNA
            rna_seq = dna_seq.transcribe()
            
            # Translate RNA to protein
            protein_seq = rna_seq.translate()

            # Write to the DNA complement file
            comp_f.write(f">Complement Sequence: {record.id}\n")
            comp_f.write(f"{dna_complement}\n\n")
            
            # Write to the RNA file
            rna_f.write(f">RNA Sequence: {record.id}\n")
            rna_f.write(f"{rna_seq}\n\n")
            
            # Write to the protein file
            prot_f.write(f">Protein Sequence: {record.id}\n")
            prot_f.write(f"{protein_seq}\n\n")

# Replace 'your_sequences.fasta' with your actual input file path
fasta_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/combined.fasta"

# Define output file paths
complement_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/dna_complement_sequences.fasta"
rna_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/rna_sequences.fasta"
protein_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/protein_sequences.fasta"

# Call the function to process the FASTA file and save output to separate files
process_fasta(fasta_file, complement_file, rna_file, protein_file)