from Bio import SeqIO
from Bio.Seq import Seq

def translate_dna_to_protein(dna_id, dna_seq):
    """
    Translate DNA sequence into six possible reading frames (3 forward, 3 reverse).
    Returns a list of protein sequences.
    """
    print(f"\nTranslating DNA sequence (ID: {dna_id}):")
    # Forward frames
    frames = [dna_seq.translate(to_stop=True), 
              dna_seq[1:].translate(to_stop=True), 
              dna_seq[2:].translate(to_stop=True)]
    
    # Reverse complement frames
    reverse_seq = dna_seq.reverse_complement()
    frames += [reverse_seq.translate(to_stop=True), 
               reverse_seq[1:].translate(to_stop=True), 
               reverse_seq[2:].translate(to_stop=True)]
    
    for i, frame in enumerate(frames, start=1):
        print(f"Frame {i}: {frame[:30]}... (truncated for readability)")

    return frames

def load_fasta_sequences(fasta_file):
    """
    Load sequences from a FASTA file using SeqIO from Biopython.
    """
    print(f"\nLoading sequences from FASTA file: {fasta_file}")
    sequences = {}
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = record.seq
            print(f"Loaded sequence {record.id} (length: {len(record.seq)})")
    return sequences

def find_matching_proteins(dna_seqs, protein_seqs):
    """
    For each DNA sequence, translate it into protein in six frames,
    and compare it to the provided protein sequences.
    """
    matches = []
    for dna_id, dna_seq in dna_seqs.items():
        print(f"\nProcessing DNA sequence {dna_id}...")
        translated_frames = translate_dna_to_protein(dna_id, dna_seq)
        
        for protein_id, protein_seq in protein_seqs.items():
            print(f"Comparing to protein sequence {protein_id} (length: {len(protein_seq)})...")
            if protein_seq in translated_frames:
                matches.append((dna_id, protein_id))
                print(f"Match found: DNA {dna_id} -> Protein {protein_id}")
            else:
                print(f"No match found for protein {protein_id}")
    
    return matches

# File paths
dna_fasta_file = '/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/combined.fasta'
protein_fasta_file = '/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/ORFs_output_filtered.fasta'

# Load DNA and protein sequences
dna_sequences = load_fasta_sequences(dna_fasta_file)
protein_sequences = load_fasta_sequences(protein_fasta_file)

# Find and print matching protein sequences
print("\nFinding matches between DNA and protein sequences...")
matches = find_matching_proteins(dna_sequences, protein_sequences)

if matches:
    print("\nMatches found:")
    for match in matches:
        print(f"DNA {match[0]} -> Protein {match[1]}")
else:
    print("\nNo matches found.")
