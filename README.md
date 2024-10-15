# VIRUS_python_scripts
Just to play with the fasta sequences

My original idea was to write a Snakefile but NO time...
In this repo you will find the following scripts


1_Generate_bigger_fasta_file.py, this script that generates a bigger fasta file using smaller fasta files

You need to change these input_dir and output_file paths 
input_dir = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Sequences_Gensam_Ours"
output_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/combined.fasta"

2_COVID_alignment.py, this script makes an alignment using the fasta sequences located in combined.fasta
It has to make the aligment between pairs.

You need to change this path 
file_path = '/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/combined.fasta'

3_Blast_DNA.py, this script will male the aligments in a database that contains viruses in my case. I had to dowload the database to my computer 
because it was taking so long so you need to do the same 

