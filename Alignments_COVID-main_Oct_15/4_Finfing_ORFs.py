from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to find ORFs in a sequence
def find_orfs(sequence, min_protein_length=100):
    orfs = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    seq_len = len(sequence)

    # Search for ORFs in the 3 reading frames of the forward strand
    for frame in range(3):
        for pos in range(frame, seq_len, 3):
            codon = str(sequence[pos:pos + 3])
            if codon == start_codon:
                for stop_pos in range(pos + 3, seq_len, 3):
                    stop_codon = str(sequence[stop_pos:stop_pos + 3])
                    if stop_codon in stop_codons:
                        # Extract ORF and translate
                        orf_seq = sequence[pos:stop_pos + 3]
                        protein_seq = orf_seq.translate(to_stop=True)
                        if len(protein_seq) >= min_protein_length:
                            orfs.append(protein_seq)
                        break
    return orfs

# Function to process a genome FASTA file and find proteins, saving them in FASTA format
def process_genome_fasta(fasta_file, output_fasta_file, min_protein_length=100):
    # Read the genome sequence from the FASTA file
    with open(fasta_file, "r") as fasta_handle:
        genome_record = SeqIO.read(fasta_handle, "fasta")
        genome_seq = genome_record.seq

        print(f"Genome length: {len(genome_seq)}")
        
        # Find ORFs and translate them to proteins
        proteins = find_orfs(genome_seq, min_protein_length)

        # Write the found proteins to a FASTA output file
        with open(output_fasta_file, "w") as out_handle:
            for i, protein in enumerate(proteins):
                # Create a SeqRecord object for each protein
                protein_record = SeqRecord(protein, 
                                           id=f"Protein_{i + 1}", 
                                           description=f"{len(protein)} amino acids")
                # Write the protein in FASTA format
                SeqIO.write(protein_record, out_handle, "fasta")

# Example usage with a genome file
genome_fasta_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/DNA_sequence_sars_2.fasta"
output_fasta_file = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/ORFs_output.fasta"
process_genome_fasta(genome_fasta_file, output_fasta_file, min_protein_length=100)
