#python3 /Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/COVID_alignment.py

from Bio import Align
from collections import defaultdict

def parse_fasta(file):
    sequences = {}
    current_sequence_id = None
    current_sequence = ''
    
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # If a sequence ID is already being processed, save the sequence
                if current_sequence_id is not None:
                    sequences[current_sequence_id] = current_sequence
                # Get the new sequence ID
                current_sequence_id = line.strip()[1:]
                current_sequence = ''
            else:
                current_sequence += line.strip()
        
        # Add the last sequence
        if current_sequence_id is not None:
            sequences[current_sequence_id] = current_sequence
    
    return sequences

# Path to the combined FASTA file
file_path = '/Users/xdisga/Documents/Python_for_biologists/Alignments_COVID-main/Merged_fasta_files/combined.fasta'  # Make sure this path is correct

sequences = parse_fasta(file_path)
#print(sequences)

# Getting the keys
keys = sequences.keys()
#print(keys)  # Output: dict_keys(['seq1', 'seq2', 'seq3'])

# If you need the keys as a list
keys = list(keys)
print(keys)
################################################################

def find_pairs(keys):
    # Dictionary to hold the pairs
    pair_dict = defaultdict(list)
    
    # Populate the dictionary with key prefixes
    for key in keys:
        prefix = key.split('_')[0]
        pair_dict[prefix].append(key)
    
    # Extract pairs
    pairs = {prefix: pair_dict[prefix] for prefix in pair_dict if len(pair_dict[prefix]) > 1}
    
    return pairs

pairs = find_pairs(keys)

for prefix, pair in pairs.items():
    print(f"Prefix: {prefix}")
    for p in pair:
        print(f"  {p}")

print("Pairing process completed.")

print(pairs) #dictionary


for key, pair in pairs.items():
    for i, sequence_id in enumerate(pair, start=1):
        print(f"sequence_id{i} = '{sequence_id}'")
    print()


#####################################################################

# Create the aligner
aligner = Align.PairwiseAligner()

# Set alignment parameters (optional)
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

# Function to calculate percentage identity excluding 'N' nucleotides
def calculate_percentage_identity_excluding_N(alignment):
    matches = 0
    valid_positions = 0
    seq1, seq2 = alignment
    for s1, s2 in zip(seq1, seq2):
        if s1 in 'ACGT' and s2 in 'ACGT':
            valid_positions += 1
            if s1 == s2:
                matches += 1
    if valid_positions == 0:
        return 0.0
    return (matches / valid_positions) * 100

# Open a file in write mode to save the output
with open('alignments_output.txt', 'w') as file:
    # Loop over the specified pairs
    for prefix, pair in pairs.items():
        file.write(f"Performing alignment for prefix: {prefix}\n")
        sequence_id1, sequence_id2 = pair
        
        selected_sequence1 = sequences.get(sequence_id1)
        selected_sequence2 = sequences.get(sequence_id2)

        if selected_sequence1 and selected_sequence2:
            # Perform the alignment
            alignments = aligner.align(selected_sequence1, selected_sequence2)

            if alignments:
                # Take the first alignment from the list
                alignment = alignments[0]

                # Print the alignment details
                file.write("Alignment:\n")
                file.write(str(alignment) + '\n')
                file.write(f"Score: {alignment.score}\n")

                # Calculate percentage identity
                alignment_length = sum(1 for s1, s2 in zip(*alignment) if '-' not in (s1, s2))
                num_matches = sum(1 for s1, s2 in zip(*alignment) if s1 == s2 and s1 != '-')
                percentage_identity = (num_matches / alignment_length) * 100
                file.write(f"Percentage Identity: {percentage_identity:.2f}%\n")

                # Calculate percentage identity excluding 'N' nucleotides
                percentage_identity_excluding_N = calculate_percentage_identity_excluding_N(alignment)
                file.write(f"Percentage Identity (excluding 'N'): {percentage_identity_excluding_N:.2f}%\n")

                # Extract and print more detailed information
                aligned_seq1, aligned_seq2 = alignment.aligned

                file.write("Aligned regions:\n")
                for (start1, end1), (start2, end2) in zip(aligned_seq1, aligned_seq2):
                    file.write(f"seq1[{start1}:{end1}] -> seq2[{start2}:{end2}]\n")
                
                file.write("Detailed alignment:\n")
                file.write(f"Length of the alignment: {alignment.shape[1]}\n")

                # Print the sequences with gaps
                file.write("Aligned sequence 1 with gaps: ")
                file.write("".join([str(c) for c in aligned_seq1]) + '\n')
                
                file.write("Aligned sequence 2 with gaps: ")
                file.write("".join([str(c) for c in aligned_seq2]) + '\n')
            else:
                file.write("No alignment found.\n")
        else:
            if not selected_sequence1:
                file.write(f"Sequence ID {sequence_id1} not found in the file.\n")
            if not selected_sequence2:
                file.write(f"Sequence ID {sequence_id2} not found in the file.\n")

        file.write('\n')  # Add an empty line between each pair
