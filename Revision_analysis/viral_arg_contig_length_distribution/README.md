# Viral ARG-carrying contig length distribution

## Purpose

This analysis summarizes the length distribution of viral contigs carrying annotated ARGs. It evaluates whether ARG-carrying viral contigs are concentrated near the 5 kb minimum viral-contig length threshold.

## Input data

- `data/viral_arg_contig_lengths_raw.tsv`: ARG-carrying viral-contig records.
- `data/unique_viral_arg_contig_lengths.csv`: unique viral contigs and lengths.

## Analysis

`scripts/analyze_viral_arg_contig_length_distribution.R` removes duplicated ARG records at the viral-contig level, bins unique viral contigs by length and exports summary statistics.

## Outputs

- `results/viral_arg_contig_length_bins.csv`: binned length distribution.
- `results/viral_arg_contig_length_summary.csv`: summary statistics for unique ARG-carrying viral contigs.
- `figures/viral_arg_contig_length_distribution.pdf`: length-distribution summary.

This analysis was used to support Response Fig. 3.
