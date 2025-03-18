#!/usr/bin/env python3
"""
GFA to FASTA Converter
==========================
This script extracts sequences from a GFA file and converts them into a FASTA format.

Author: Franck Dumetz
Date: 2025-03-18
Version: 1.0

Usage:
    python gfa2fasta.py -i input.gfa -o output.fasta

Dependencies:
    - Python 3.8

Arguments:
    -i, --input  : Input GFA file (required)
    -o, --output : Output FASTA file (required)

Example:
    python gfa2fasta.py -i assembly.gfa -o assembly.fasta
"""

import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Convert GFA to FASTA")
parser.add_argument("-i", "--input", required=True, help="Input GFA file")
parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
args = parser.parse_args()

# Open input and output files
with open(args.output, "w") as fasta_file:
    with open(args.input, "r") as gfa_file:
        for line in gfa_file:
            cols = line.strip().split("\t")
            if cols[0] == "S":  # Only take 'S' lines
                fasta_file.write(f">{cols[1]}\n{cols[2]}\n")

print(f"✅ GFA to FASTA conversion completed: {args.output}")
