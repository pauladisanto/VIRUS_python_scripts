from Bio import SeqIO

def filter_truncated_proteins(input_fasta, output_fasta):
    # Dictionary to store the longest sequence for each unique protein
    longest_sequences = {}

    # Read sequences from the input FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_str = str(record.seq)
        is_truncated = False
        
        # Check if the current sequence is a truncated version of an existing sequence
        for stored_seq in list(longest_sequences.keys()):
            if seq_str in stored_seq:  # Current sequence is a truncated version of an existing sequence
                is_truncated = True
                break
            elif stored_seq in seq_str:  # Stored sequence is a truncated version of the current sequence
                longest_sequences.pop(stored_seq)  # Remove the truncated sequence

        if not is_truncated:
            longest_sequences[seq_str] = record
            print(f"Stored sequence {record.id} with length {len(record.seq)}")

    # Write the longest sequences to the output FASTA file
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(longest_sequences.values(), output_handle, "fasta")

# Usage
input_fasta = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/ORFs_output.fasta"
output_fasta = "/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/ORFs_output_filtered.fasta"

filter_truncated_proteins(input_fasta, output_fasta)
