VIRUS Python Scripts

This repository contains Python scripts designed to manipulate and analyze FASTA sequences, particularly for viral genome analysis. The initial goal was to use Snakemake, but due to time constraints, I’ve opted for Python scripts. These scripts were developed in a conda environment with Python 3 and Biopython installed.
Prerequisites

    Install Conda.
    Create and activate a new environment with the required dependencies:

    conda create -n virus_analysis python=3 biopython
    conda activate virus_analysis

Included Scripts
1. Generate_bigger_fasta_file.py

This script merges multiple smaller FASTA files into one larger FASTA file.

    Modify these paths before running:

    input_dir = "/path/to/your/smaller/fasta/files"
    output_file = "/path/to/output/combined.fasta"

2. COVID_alignment.py

This script performs pairwise alignment on the sequences found in the combined.fasta file.

    Modify this path:

    file_path = '/path/to/combined.fasta'

3. Blast_DNA.py

This script runs local BLAST alignments against a downloaded virus database. Since online BLAST can be slow, it's recommended to download the database locally.

    Modify these paths:

    blast_path = "/path/to/blast/bin"  # Location of your BLAST installation
    blast_results_folder = "/path/to/blast_results/"  # Folder to store BLAST results
    fasta_file_path = "/path/to/combined.fasta"  # Path to your combined FASTA file

4. Finding_ORFs.py

This script identifies open reading frames (ORFs) in the SARS-CoV-2 genome. The genome file DNA_sequence_sars_2.fasta is included in the repository.

    Modify these paths:

    genome_fasta_file = "/path/to/DNA_sequence_sars_2.fasta"
    output_fasta_file = "/path/to/ORFs_output.fasta"  # Output containing identified proteins

5. Removing_duplicate_ORFs.py

This script filters out duplicate or truncated proteins from the ORFs_output.fasta file. Due to the nature of the ORF function, some sequences may include partial proteins, which are removed here.

    Modify these paths:

    input_fasta = "/path/to/ORFs_output.fasta"
    output_fasta = "/path/to/ORFs_output_filtered.fasta"

6. Protein_alignment.py

This script aligns the proteins in the filtered ORFs_output_filtered.fasta file against the NCBI database for SARS-CoV-2 proteins. Note that this alignment may take a long time even with a reduced database. (You can just avoid the script). This part is just to see in the found ORF makes sense.  

    Modify these paths:

    query_protein_fasta = "/path/to/ORFs_output_filtered.fasta"
    output_dir = "/path/to/blast_proteins"

7. Alignment_my_seq_and_ORFs.py

This script aligns your DNA sequences from combined.fasta with the proteins obtained from the ORF function. It first translates the DNA into protein sequences (considering all six possible reading frames) and then compares these to the ORFs in ORFs_output_filtered.fasta.

    Modify these paths:

    dna_fasta_file = '/path/to/combined.fasta'
    protein_fasta_file = '/path/to/ORFs_output_filtered.fasta'

Notes:

    My sequences contain many "NNNN" nucleotides, which may reduce alignment accuracy and prevent finding 100% matches. Additionally, some proteins may lack stop codons, further complicating alignments.
