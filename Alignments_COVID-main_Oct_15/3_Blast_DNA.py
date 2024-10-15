import subprocess
import os
from Bio import SeqIO
from Bio.Blast import NCBIXML

# Ensure the PATH includes the directory where blastn is located
blast_path = "/usr/local/ncbi/blast/bin"  # Update this path based on the output of `which blastn`
os.environ["PATH"] += os.pathsep + blast_path

# Full path to the BLAST database
blast_db_path = "/Users/xdisga/ref_viruses_rep_genomes"  # Update this path to where the database files are located

# Path to store BLAST results
blast_results_folder = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Blast_results/"
os.makedirs(blast_results_folder, exist_ok=True)  # Create the folder if it doesn't exist

# Function to perform local BLAST search for a sequence
def blast_local(seq_record, output_file):
    print(f"Performing BLAST for sequence ID: {seq_record.id}")
    try:
        # Write the sequence to a temporary FASTA file
        fasta_file = f"temp_{seq_record.id}.fasta"
        with open(fasta_file, "w") as f:
            f.write(f">{seq_record.id}\n{str(seq_record.seq)}\n")

        # Run the BLAST command locally
        result_file = os.path.join(blast_results_folder, f"blast_results_{seq_record.id}.xml")
        blast_command = [
            "blastn",
            "-query", fasta_file,
            "-db", blast_db_path,
            "-out", result_file,
            "-outfmt", "5"  # XML format
        ]
        print(f"Running command: {' '.join(blast_command)}")
        subprocess.run(blast_command, check=True)

        # Parse the BLAST results
        with open(result_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)

            # Write details about the BLAST results to the output file
            with open(output_file, "a") as out_f:
                out_f.write(f"\n#### Results for Sample: {seq_record.id} ####\n")
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            out_f.write(f"****Alignment****\n")
                            out_f.write(f"Sequence: {alignment.title}\n")
                            out_f.write(f"Length: {alignment.length}\n")
                            out_f.write(f"E value: {hsp.expect}\n")
                            out_f.write(f"Query Start: {hsp.query_start} - Query End: {hsp.query_end}\n")
                            out_f.write(f"Subject Start: {hsp.sbjct_start} - Subject End: {hsp.sbjct_end}\n")
                            out_f.write(f"Score: {hsp.score}\n")
                            out_f.write(f"Identities: {hsp.identities}/{hsp.align_length}\n")
                            out_f.write(f"Gaps: {hsp.gaps}\n")

                            # Write the complete alignment (full query, match, and subject sequences)
                            out_f.write(f"Query Sequence: {hsp.query}\n")
                            out_f.write(f"Match:         {hsp.match}\n")
                            out_f.write(f"Subject:       {hsp.sbjct}\n")
                            out_f.write("\n")
        
        # Remove the temporary FASTA file after BLAST is done
        os.remove(fasta_file)
        print(f"Temporary file {fasta_file} removed.")
        
    except Exception as e:
        print(f"Error processing sequence ID {seq_record.id}: {e}")
        if os.path.exists(fasta_file):
            os.remove(fasta_file)  # Ensure temporary files are cleaned up even on error

# Read sequences from a file (FASTA format)
fasta_file_path = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/combined.fasta"
print(f"Reading sequences from {fasta_file_path}")
sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

# Check if sequences are read correctly
if not sequences:
    print("No sequences found in the FASTA file.")
else:
    print(f"Total sequences read: {len(sequences)}")

# Output file for BLAST results
output_file = "blast_results.txt"

# Process each sequence individually
for seq_record in sequences:
    blast_local(seq_record, output_file)
