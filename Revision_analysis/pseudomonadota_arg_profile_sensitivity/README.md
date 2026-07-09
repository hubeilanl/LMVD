# Pseudomonadota ARG profile and short-contig sensitivity

## Purpose

This analysis compares ARG profiles between Pseudomonadota AR-MAGs and livestock-derived Pseudomonadota isolate genomes, and evaluates the sensitivity of ARG counts to short-contig filtering. It is intended to assess whether elevated ARG counts in Pseudomonadota AR-MAGs are consistent with isolate-genome patterns and robust to removal of short contigs.

## Input data

- `data/pseudomonadota_arg_type_proportions_by_genome.tsv`: per-MAG/genome ARG-type proportions.
- `data/pseudomonadota_short_contig_sensitivity.tsv`: ARG-count sensitivity to contig-length cutoffs.
- `data/pseudomonadota_arg_composition_support.tsv`: summary support table for genus-level and ARG-type composition statements.

## Analysis

`scripts/analyze_pseudomonadota_arg_profile_sensitivity.R` summarizes ARG-type proportions and short-contig filtering results from the included source tables.

## Outputs

- `results/`: ARG-type proportion summaries, short-contig sensitivity tables and text summaries.
- `figures/`: PDF summaries of ARG-type profiles and short-contig sensitivity.

The folder includes source-data tables for this analysis, not the full upstream isolate-genome annotation table.

This analysis was used to support Response Figs. 4 and 5.
