from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import concurrent.futures
import os

# Function to run BLASTp remotely using NCBI's online service for SARS-CoV-2
def run_blastp_ncbi(record, index, output_dir):
    fasta_data = f">{record.id}\n{str(record.seq)}\n"
    print(f"Submitting BLASTp search to NCBI for sequence {record.id} (length: {len(record.seq)})")

    try:
        # Taxonomy ID for SARS-CoV-2 is 2697049, using the refseq_select database
        result_handle = NCBIWWW.qblast("blastp", "refseq_protein", fasta_data, 
                                       expect=1e-3, hitlist_size=20, 
                                       entrez_query="txid2697049[ORGN]")  # SARS-CoV-2 filter

        # Save the results in XML format
        output_file = os.path.join(output_dir, f"blast_result_{index+1}.xml")
        with open(output_file, "w") as out_file:
            out_file.write(result_handle.read())
        
        print(f"BLAST search for {record.id} completed. Results saved to {output_file}")
    except Exception as e:
        print(f"Error submitting BLASTp search for sequence {record.id}: {e}")

# Function to parse BLAST XML results
def parse_blast_results(xml_file):
    with open(xml_file, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            print(f"\nQuery: {blast_record.query}")
            if not blast_record.alignments:
                print("No alignments found.")
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    print(f"****Alignment****")
                    print(f"Sequence: {alignment.title}")
                    print(f"Length: {alignment.length}")
                    print(f"Score: {hsp.score}")
                    print(f"E-value: {hsp.expect}")
                    print(f"Identities: {hsp.identities}")
                    print(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n")

# Protein FASTA file
query_protein_fasta = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/ORFs_output_filtered.fasta"
output_dir = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/blast_proteins"

# Create the output directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

# Read the protein sequences from the FASTA file
protein_sequences = list(SeqIO.parse(query_protein_fasta, "fasta"))

# Check if sequences are read correctly (just in case)
if not protein_sequences:
    print("No sequences found in the FASTA file.")
else:
    print(f"Total sequences read: {len(protein_sequences)}")

# Run BLASTp for each sequence in parallel (limit the number of threads)
max_workers = 3  # Limit the number of concurrent BLAST searches to avoid overloading NCBI
with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(run_blastp_ncbi, record, i, output_dir) for i, record in enumerate(protein_sequences)]
    concurrent.futures.wait(futures)

# Parse BLAST results (once the results are saved)
for i in range(len(protein_sequences)):
    parse_blast_results(os.path.join(output_dir, f"blast_result_{i+1}.xml"))
    