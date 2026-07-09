# ARG-MGE non-viral counterpart control

## Purpose

This analysis tests whether selected viral ARG-MGE units have same-sample counterparts on non-viral contigs. The comparison provides a targeted control for distinguishing viral ARG-MGE associations from co-occurring non-viral mobilome signals in the same metagenomes.

## Input-table scope

The input tables are target-pair-filtered tables rather than complete annotation tables for all ARGs, MGEs or contigs. They contain only records relevant to the four ARG-MGE pairs listed below. In the file names, `all_contig` indicates that the records were extracted from the all-contig annotation background before the non-viral counterpart check.

## Input data

- `data/target_pair_viral_arg_annotations.tsv`: viral ARG annotations filtered for the target ARGs.
- `data/target_pair_viral_mge_annotations.tsv`: viral MGE annotations filtered for the target MGEs.
- `data/target_pair_all_contig_arg_annotations.tsv`: target ARG annotations extracted from the all-contig background.
- `data/target_pair_all_contig_mge_annotations.tsv`: target MGE annotations extracted from the all-contig background.

## Target ARG-MGE pairs

- `tnp + lnuC`
- `tndX + AAC(6')-Im`
- `tndX + APH(2'')-IIa`
- `tnpR_2 + ACI-1`

## Analysis

`scripts/analyze_arg_mge_nonviral_counterparts.py` identifies viral target ARG-MGE units, searches the corresponding all-contig background for same-sample non-viral counterparts and summarizes counterpart status by unit and by ARG-MGE pair.

## Outputs

- `results/viral_arg_mge_units_with_nonviral_counterpart_status.tsv`: unit-level viral ARG-MGE counterpart status.
- `results/same_sample_nonviral_candidate_pairs.tsv`: same-sample non-viral candidate pairs for the target combinations.
- `results/nonviral_counterpart_summary_by_pair.tsv`: pair-level summary.
- `results/nonviral_counterpart_summary_overall.tsv`: overall summary.

This analysis was used as a targeted control for the Fig. 7d-f ARG-MGE units.
