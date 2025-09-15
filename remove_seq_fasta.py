#!/usr/bin/env python

### Removing a sequence from a fasta file knowing the header

from Bio import SeqIO

# Input and output file paths
input_fasta = "PATH/PlasmoDB-68_PknowlesiH_Genome.fasta"
output_fasta = "PATH/PlasmoDB-68_PknowlesiH_Genome-MIT.fasta"
header_to_remove = "AY722797.1 Plasmodium knowlesi isolate Malayan MRA-487 mitochondrion, complete genome"  # Replace with the exact header (without the '>')

# Read the input FASTA and filter sequences
with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
    records_to_keep = (
        record for record in SeqIO.parse(input_handle, "fasta") 
        if record.id != header_to_remove
    )
    # Write the filtered sequences to the output file
    SeqIO.write(records_to_keep, output_handle, "fasta")

print(f"Sequence with header '{header_to_remove}' has been removed.")
print(f"Filtered FASTA file saved to {output_fasta}.")
