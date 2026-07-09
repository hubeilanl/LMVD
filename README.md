# A global livestock virome catalog links viruses with resistome assembly

Livestock manure is a major One Health interface for antimicrobial resistance (AMR), yet the viral component of this ecosystem has remained poorly resolved. This study integrates 4,017 publicly available livestock manure metagenomes to establish the Livestock Manure Virome Database (LMVD), a global catalogue of 1,330,372 species-level viral operational taxonomic units (vOTUs). More than 99% of these vOTUs were not represented in existing viral genome resources, showing that livestock manure contains a large, habitat-specific viral sequence space that has been underrepresented in previous databases.

The analysis links this viral catalogue with microbial hosts, antibiotic resistance genes (ARGs), mobile genetic elements and environmental predictors. LMVD reveals restricted biogeographic distributions of livestock manure viruses and identifies antimicrobial usage and farming scale as important correlates of viral alpha diversity. Although viral genomes contained 78 ARG subtypes, representing 9.1% of the livestock manure resistome diversity, direct viral ARG carriage was uncommon. The study therefore places boundaries on the role of viruses as direct ARG carriers while highlighting a broader route through which viruses may shape AMR risk, namely predicted interactions between temperate viruses and ARG-carrying or pathogenic bacterial hosts.

This repository provides the command records, R scripts and processed analysis tables used for the computational and statistical analyses described in the manuscript. It is organized into three components: bioinformatics command templates, statistical analysis scripts and additional sensitivity or control analyses.

## Repository structure

| Directory | Contents | Scope |
|---|---|---|
| `bioinformatics_command_templates/` | Command records for the main sequence-based analyses. | Covers viral discovery, viral genome curation, vOTU clustering, viral taxonomy and lifestyle prediction, MAG reconstruction, mOTU profiling, virus-host linkage prediction, functional annotation, plasmid prediction, Ka/Ks analysis and abundance profiling. |
| `statistical_analysis/` | R scripts and processed input tables for manuscript-level statistical analyses. | Covers environmental predictor analysis, microbial and viral alpha diversity, Bray-Curtis dissimilarity, PCoA/PERMANOVA, Procrustes analysis, distance decay and latitudinal diversity gradients. |
| `Revision_analysis/` | Source tables, scripts and summaries for additional sensitivity and control analyses. | Includes targeted analyses of vOTU prevalence-abundance patterns, technical sensitivity of viral Shannon diversity, viral ARG-carrying contig lengths, Pseudomonadota ARG-profile sensitivity and non-viral counterparts of selected ARG-MGE units. |

## Bioinformatics command templates

The `bioinformatics_command_templates/` directory records the principal command structures used to generate the sequence-derived resources analysed in the manuscript. The files are grouped by analytical task rather than by sample, so that the major tools, parameters and input-output relationships can be inspected directly.

`viral_identification_taxonomy_lifestyle_commands.md` describes the recovery of putative viral contigs from assembled metagenomes using VirSorter2, CheckV-based trimming and quality assessment, SOP-based contig curation, vOTU clustering at the species level, taxonomic assignment and lifestyle prediction.

`mag_reconstruction_motu_profiling_commands.md` documents MAG reconstruction and refinement, CheckM2-based genome quality assessment, dRep-based mOTU dereplication, GTDB-Tk taxonomic classification, read-based abundance profiling and growth-rate estimation.

`virus_host_linkage_prediction_commands.md` describes the two-stage host-prediction strategy. Candidate virus-host pairs were inferred from CRISPR spacer matching, nucleotide sequence homology and tRNA matching, and then evaluated using VirMatcher, which integrates BLAST, CRISPR, tRNA and WIsH evidence.

`functional_annotation_and_sequence_profiling_commands.md` records the command structures for ARG annotation, MGE annotation, virulence-factor annotation, plasmid prediction, Ka/Ks analysis and contig-level abundance estimation across assembled contigs, viral contigs, vOTUs and MAGs.

These command records are intended to document the analytical implementation used in the study. Paths to local files, reference databases, software environments and computing-cluster settings should be specified according to the user's computing environment.

## Statistical analysis scripts

The `statistical_analysis/` directory contains processed tables and R scripts for the statistical analyses included in the manuscript.

`random_forest_environmental_predictors/` contains the sample-level viral Shannon diversity table, country-level environmental predictor summaries, predictor definitions and an R script for random forest regression with permutation-based importance testing. This module corresponds to the analysis of environmental predictors associated with viral alpha diversity, including antimicrobial usage, farming scale, manure production, meat production, human population density and related socioeconomic variables.

`microbial_viral_community_diversity/` contains compressed vOTU and mOTU abundance matrices, sample metadata and R scripts for microbial and viral community analyses. The scripts estimate richness and Shannon diversity, calculate Bray-Curtis dissimilarity matrices, perform PCoA and PERMANOVA, assess viral-microbial community concordance, and evaluate spatial patterns including distance decay and latitudinal richness gradients.

## Sensitivity and control analyses

The `Revision_analysis/` directory contains analysis-specific modules with source data, scripts and output summaries for targeted sensitivity or control analyses. These modules address:

- the relationship between vOTU prevalence and mean abundance across manure types;
- continental differences in viral Shannon diversity after accounting for sequencing depth and assembly-output metrics;
- the length distribution of ARG-carrying viral contigs;
- ARG-profile sensitivity in Pseudomonadota AR-MAGs using livestock-derived Pseudomonadota isolate genomes and short-contig filtering;
- same-sample non-viral counterpart checks for selected viral ARG-MGE units.

Each subdirectory contains a concise README describing the analysis purpose, input tables, script and output files.

## External data resources

LMVD is the main data resource generated by this study. It provides a species-level viral genome catalogue for livestock manure, together with associated metadata, and is designed to support future analyses of viral diversity, viral read recruitment, virus-host associations and the role of viruses in livestock-associated AMR ecology. Viral genomes and associated metadata from LMVD are freely accessible through Zenodo:

- Livestock Manure Virome Database: [https://doi.org/10.5281/zenodo.17226303](https://doi.org/10.5281/zenodo.17226303)

Candidate ARG analyses in this manuscript use the established livestock manure candidate ARG set from the preceding livestock manure resistome study. That study analysed the same broad livestock manure metagenomic resource at the resistome level and identified latent ARGs that are not captured by conventional ARG databases but may encode antibiotic resistance functions. The corresponding resistome data resources are available through Zenodo:

- Global health risks lurking in livestock resistome: [https://doi.org/10.5281/zenodo.15025586](https://doi.org/10.5281/zenodo.15025586)

Raw metagenomic sequencing data were obtained from public repositories including NCBI, ENA and CNSA, as listed in the manuscript Supplementary Data files. Country-level environmental predictors were derived from public sources described in the manuscript, including FAOSTAT, OECD and Our World in Data. Large raw reads, assemblies and intermediate high-performance-computing outputs are not redistributed in this repository.

## Software notes

R scripts are written for the processed input tables included in the relevant analysis folders. Some abundance matrices are large and may require a high-memory computing environment. Bioinformatics command records refer to third-party software and databases, including VirSorter2, CheckV, geNomad, PhaGCN2, metaWRAP, CheckM2, dRep, GTDB-Tk, CoverM, VirMatcher, RGI/CARD, DIAMOND, mobileOG, VFDB and plasmid-prediction tools.
