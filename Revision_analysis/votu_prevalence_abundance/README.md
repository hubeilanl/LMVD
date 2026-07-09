# vOTU prevalence-abundance analysis

## Purpose

This analysis examines the relationship between the prevalence of livestock-manure vOTUs and their mean abundance in samples where they were detected. The analysis is performed separately for chicken, swine and cattle manure to assess whether low cross-type sharing is accompanied by consistently low abundance.

## Input data

- `data/votu_prevalence_abundance_source_data.tsv`: vOTUs with prevalence greater than 1 percent within chicken, swine or cattle manure.

## Analysis

`scripts/analyze_votu_prevalence_abundance.R` fits manure-type-specific linear models of mean abundance in present samples against prevalence percentage and exports the fitted statistics.

## Outputs

- `results/votu_prevalence_abundance_lm_statistics.tsv`: regression statistics for each manure type.
- `figures/`: PDF summaries of the prevalence-abundance relationships.

This analysis was used to support Response Fig. 1.
