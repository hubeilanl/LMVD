# Viral Shannon diversity technical-sensitivity analysis

## Purpose

This analysis evaluates whether continental differences in assembly-based viral Shannon diversity remain detectable after accounting for sequencing depth and assembly-output metrics. It provides a technical-sensitivity check for diversity patterns derived from assembled viral contigs.

## Input data

- `data/viral_shannon_technical_metrics.tsv`: sample-level viral Shannon diversity and technical metrics.

## Model

`scripts/analyze_viral_shannon_technical_sensitivity.R` fits:

```text
Shannon index ~ log10(nReads + 1) + log10(sum_len + 1) + log10(num_seqs + 1)
```

Residual viral Shannon diversity is then compared among continents. This is a sensitivity analysis of assembly-derived diversity estimates, not full catalogue-wide read recruitment.

## Outputs

- `results/`: residuals, continent-level summaries, test statistics and model diagnostics.
- `figures/viral_shannon_residual_by_continent.pdf`: residual Shannon diversity by continent.

This analysis was used to support Response Fig. 2.
