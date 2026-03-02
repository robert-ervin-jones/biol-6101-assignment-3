"""
Aggregate individual featureCounts output files into a single count matrix.

This script is called by Snakemake and has access to:
- snakemake.input: list of input files
- snakemake.output: list of output files  
- snakemake.log: log file path
"""

import pandas as pd
import sys

# Set up logging
sys.stderr = open(snakemake.log[0], "w")

def parse_featurecounts(file_path):
    """Parse a featureCounts output file and extract gene counts."""
    # Read the featureCounts output (skip comment lines starting with #)
    df = pd.read_csv(file_path, sep="\t", comment="#")
    
    # Extract sample name from file path
    sample_name = file_path.split("/")[-1].replace(".counts.txt", "")
    
    # Get the last column (counts) and set gene ID as index
    counts = df.iloc[:, -1]
    counts.index = df.iloc[:, 0]
    counts.name = sample_name
    
    return counts

# Main processing
print(f"Processing {len(snakemake.input)} count files...", file=sys.stderr)

# Read all count files
count_series = []
for count_file in snakemake.input:
    print(f"Reading {count_file}...", file=sys.stderr)
    counts = parse_featurecounts(count_file)
    count_series.append(counts)

# Combine into a single dataframe
count_matrix = pd.concat(count_series, axis=1)
count_matrix = count_matrix.sort_index()

print(f"Count matrix shape: {count_matrix.shape}", file=sys.stderr)
print(f"Genes: {count_matrix.shape[0]}, Samples: {count_matrix.shape[1]}", file=sys.stderr)

# Save to output file
count_matrix.to_csv(snakemake.output[0], sep="\t")
print(f"Count matrix saved to {snakemake.output[0]}", file=sys.stderr)

