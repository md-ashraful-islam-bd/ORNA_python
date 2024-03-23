# FASTA Analysis Tool

## Description
This Python script performs quality filtering and k-mer frequency analysis on sequences from a FASTA file. It filters sequences based on average quality scores and median k-mer frequencies, and outputs the filtered sequences to a new FASTA file.

## Features
- Read FASTA files
- Calculate k-mer frequencies
- Filter sequences based on quality scores
- Output filtered sequences to a new FASTA file

## Requirements
- Python 3.x

## Usage
```bash
python fasta_analysis.py <input_fasta_file> <output_fasta_file> <quality_threshold> <k>

<input_fasta_file>: Path to the input FASTA file.
<output_fasta_file>: Path to the output FASTA file.
<quality_threshold>: Quality threshold for filtering sequences (float).
<k>: Length of k-mer (integer)

Example: python fasta_analysis.py input.fasta output.fasta 30 5

Input File Format
The input FASTA file should follow the standard FASTA format:
>Header1
ATCGATCGATCG...
>Header2
ATCGATCGATCG...
...


Output File Format
The output FASTA file will contain the filtered sequences along with their average quality scores and median k-mer frequencies:
>Header1 avg_quality_score=30.5 median_kmer_freq=2
ATCGATCGATCG...
>Header2 avg_quality_score=32.0 median_kmer_freq=3
ATCGATCGATCG...
...


License
This project is licensed under the Apache License 2.0. See the LICENSE file for details.

Copyright
Copyright (c) 2024 Md. Ashraful Islam. All rights reserved.
