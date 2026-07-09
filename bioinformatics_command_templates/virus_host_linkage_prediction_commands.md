# Virus-host linkage prediction for vOTUs

This file documents representative command structures for predicting virus-host linkages between livestock manure vOTUs and in situ microbial hosts. The workflow follows a two-stage strategy:

1. Candidate virus-host pairs were first inferred using CRISPR spacer matching, nucleotide sequence homology and tRNA matching.
2. Candidate pairs were then evaluated using VirMatcher, which integrates BLAST, CRISPR, tRNA and WIsH evidence. Only pairs with an interpretation score >= 3 were retained.

The commands are provided to document the principal tools, filtering criteria and input/output relationships used in the analysis. They are not intended to constitute an end-to-end executable workflow for all 4,017 samples, and local paths, software installations and sample identifiers should be specified according to the computing environment.

## 1. Inputs and directory structure

Representative input files:

```bash
VOTU_FASTA="LMVD_vOTU_clustering/LMVD_vOTUs.min5kb.fa"
MOTU_REP_DIR="mOTU_dereplication/dRep/all_MAGs_dRep/dereplicated_genomes"
MOTU_REP_FASTA="host_prediction/in_situ_hosts/LMVD_mOTU_representatives.fa"
THREADS=32
```

Representative working directories:

```text
host_prediction/
├── in_situ_hosts/              # representative mOTU genomes used as in situ host candidates
├── crispr/                     # CRISPR spacer extraction and spacer-vOTU matching
├── nucleotide_homology/        # vOTU-mOTU nucleotide sequence homology
├── trna/                       # tRNA prediction and tRNA matching
├── candidate_pairs/            # merged candidate virus-host pairs
└── VirMatcher/                 # VirMatcher verification output
```

## 2. In situ host genome set

The host candidate set consisted of representative mOTU genomes obtained from CheckM2-filtered MAGs and dRep dereplication. These representative genomes were used as the in situ microbial host database.

Prepare one combined host FASTA file when a single host database FASTA is required:

```bash
mkdir -p host_prediction/in_situ_hosts

cat ${MOTU_REP_DIR}/*.fa \
  > host_prediction/in_situ_hosts/LMVD_mOTU_representatives.fa
```

VirMatcher requires the bacterial host directory to contain host genomes as FASTA files. A host taxonomy file should also be provided as a tab-delimited table:

```text
host_name<TAB>taxonomy
```

Representative taxonomy records:

```text
mOTU_000001    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria
mOTU_000002    d__Bacteria;p__Bacteroidota;c__Bacteroidia
```

## 3. CRISPR spacer matching

CRISPR spacers were extracted from representative mOTU genomes using the MGV CRISPR spacer scripts, which call CRT and PILER-CR and merge predicted CRISPR arrays.

Representative command structure for one split host FASTA file:

```bash
mkdir -p host_prediction/crispr/data/xaa

cd /path/to/MGV/crispr_spacers

python identify_crispr.py \
  -i /path/to/13k_split_MAGs/xaa \
  -o host_prediction/crispr/data/xaa

python merge_crispr.py \
  host_prediction/crispr/data/xaa/crt \
  host_prediction/crispr/data/xaa/pilercr \
  host_prediction/crispr/data/xaa/merged
```

Format merged CRISPR spacers as FASTA:

```bash
cat host_prediction/crispr/data/xaa/merged.spacers \
  | cut -f 1,2,3,6 \
  | awk -v FS="\t" -v OFS="\t" '{print ">",$1,"crispr",$2,$3,$4}' \
  | sed 1d \
  | sed 's/\t/_/' \
  | sed 's/\t/_/' \
  | sed 's/\t/_/' \
  | sed 's/\t/_/' \
  | sed 's/_//' \
  | awk '{if($2!="")print $0}' \
  | sed 's/\t/\n/' \
  > host_prediction/crispr/data/xaa/xaa.crispr.fa
```

After processing all split files, concatenate all CRISPR spacer FASTA files:

```bash
cat host_prediction/crispr/data/*/*.crispr.fa \
  > host_prediction/crispr/LMVD_mOTU_CRISPR_spacers.fa
```

Build a BLAST database from CRISPR spacers:

```bash
makeblastdb \
  -dbtype nucl \
  -in host_prediction/crispr/LMVD_mOTU_CRISPR_spacers.fa \
  -parse_seqids \
  -out host_prediction/crispr/LMVD_mOTU_CRISPR_spacers
```

Match vOTUs against CRISPR spacers:

```bash
blastn \
  -query ${VOTU_FASTA} \
  -db host_prediction/crispr/LMVD_mOTU_CRISPR_spacers \
  -out host_prediction/crispr/vOTU_vs_CRISPR.tsv \
  -evalue 1e-5 \
  -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  -num_threads ${THREADS}
```

Retain perfect spacer matches:

```bash
awk '$3==100 && $4/$14==1 {print $0}' \
  host_prediction/crispr/vOTU_vs_CRISPR.tsv \
  > host_prediction/crispr/vOTU_vs_CRISPR.perfect_matches.tsv
```

## 4. Nucleotide sequence homology matching

Nucleotide sequence homology between vOTUs and representative mOTU genomes was evaluated using BLASTn.

Build a BLAST database from vOTUs:

```bash
mkdir -p host_prediction/nucleotide_homology

makeblastdb \
  -dbtype nucl \
  -in ${VOTU_FASTA} \
  -out host_prediction/nucleotide_homology/LMVD_vOTUs
```

Align mOTU representative genomes to vOTUs:

```bash
blastn \
  -query ${MOTU_REP_FASTA} \
  -db host_prediction/nucleotide_homology/LMVD_vOTUs \
  -out host_prediction/nucleotide_homology/mOTU_vs_vOTU.blast.tsv \
  -evalue 1e-5 \
  -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  -num_threads ${THREADS}
```

Sort BLAST hits by query and bit score:

```bash
sort -k1,1 -k12,12nr \
  host_prediction/nucleotide_homology/mOTU_vs_vOTU.blast.tsv \
  > host_prediction/nucleotide_homology/mOTU_vs_vOTU.blast.sorted.tsv
```

Select best hits per query, if a best-hit script is used:

```bash
python scripts/get_best_hit_contig.py \
  host_prediction/nucleotide_homology/mOTU_vs_vOTU.blast.sorted.tsv \
  host_prediction/nucleotide_homology/mOTU_vs_vOTU.blast.best.tsv
```

Retain homologous matches with alignment length >= 2,500 bp, sequence identity >= 90% and alignment covering >= 75% of the vOTU sequence:

```bash
awk '$4>=2500 && $3>=90 && $4/$14>=0.75 {print $0}' \
  host_prediction/nucleotide_homology/mOTU_vs_vOTU.blast.best.tsv \
  > host_prediction/nucleotide_homology/mOTU_vs_vOTU.homology_filtered.tsv
```

Note: in this BLAST output format, column 14 is the subject length. Because vOTUs were used as the BLAST database, the subject sequence corresponds to the vOTU.

## 5. tRNA matching

tRNAs were predicted from both representative mOTU genomes and vOTUs using tRNAscan-SE with the `-G -Q` parameters.

Predict tRNAs from mOTU representative genomes:

```bash
mkdir -p host_prediction/trna/mOTUs

tRNAscan-SE \
  -G \
  -Q \
  ${MOTU_REP_FASTA} \
  -o host_prediction/trna/mOTUs/LMVD_mOTUs.trnascan.out \
  -f host_prediction/trna/mOTUs/LMVD_mOTUs.trnascan.ss \
  -s host_prediction/trna/mOTUs/LMVD_mOTUs.trnascan.iso \
  -m host_prediction/trna/mOTUs/LMVD_mOTUs.trnascan.stats \
  -b host_prediction/trna/mOTUs/LMVD_mOTUs.trnascan.bed \
  -a host_prediction/trna/mOTUs/LMVD_mOTUs.tRNA.fa \
  --thread 1
```

Predict tRNAs from vOTUs:

```bash
mkdir -p host_prediction/trna/vOTUs

tRNAscan-SE \
  -G \
  -Q \
  ${VOTU_FASTA} \
  -o host_prediction/trna/vOTUs/LMVD_vOTUs.trnascan.out \
  -f host_prediction/trna/vOTUs/LMVD_vOTUs.trnascan.ss \
  -s host_prediction/trna/vOTUs/LMVD_vOTUs.trnascan.iso \
  -m host_prediction/trna/vOTUs/LMVD_vOTUs.trnascan.stats \
  -b host_prediction/trna/vOTUs/LMVD_vOTUs.trnascan.bed \
  -a host_prediction/trna/vOTUs/LMVD_vOTUs.tRNA.fa \
  --thread 1
```

Build a BLAST database from viral tRNAs:

```bash
makeblastdb \
  -dbtype nucl \
  -in host_prediction/trna/vOTUs/LMVD_vOTUs.tRNA.fa \
  -out host_prediction/trna/vOTUs/LMVD_vOTUs.tRNA
```

Compare mOTU tRNAs against vOTU tRNAs:

```bash
blastn \
  -query host_prediction/trna/mOTUs/LMVD_mOTUs.tRNA.fa \
  -db host_prediction/trna/vOTUs/LMVD_vOTUs.tRNA \
  -out host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.tsv \
  -outfmt '6 std qlen slen' \
  -num_threads 8
```

Sort hits and retain best hit per query:

```bash
sort -k1,1 -k12,12nr \
  host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.tsv \
  > host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.sorted.tsv

python scripts/get_best_hit_contig.py \
  host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.sorted.tsv \
  host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.best.tsv
```

Retain full-length exact tRNA matches between mOTUs and vOTUs:

```bash
awk '$3==100 && $4/$13==1 && $4/$14==1 {print $0}' \
  host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.best.tsv \
  > host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.full_length_exact.tsv
```

## 6. Merge candidate virus-host pairs

Candidate virus-host pairs from CRISPR spacer matching, nucleotide homology and tRNA matching were merged according to the corresponding output-table formats. This step was performed by table parsing and formatting rather than by a single standardized standalone script.

The expected inputs were:

```text
host_prediction/crispr/vOTU_vs_CRISPR.perfect_matches.tsv
host_prediction/nucleotide_homology/mOTU_vs_vOTU.homology_filtered.tsv
host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.full_length_exact.tsv
```

The merged candidate-pair table was formatted with at least the following fields:

```text
vOTU_id    host_id    evidence_type    identity    alignment_length    coverage    source_table
```

Representative output:

```text
host_prediction/candidate_pairs/LMVD_candidate_virus_host_pairs.tsv
```

Representative table-formatting logic:

```bash
mkdir -p host_prediction/candidate_pairs

# Extract candidate pairs from CRISPR spacer matches.
# The exact host_id field depends on the CRISPR spacer header format.
awk 'BEGIN{OFS="\t"} {print $1,$2,"CRISPR",$3,$4,$4/$14,"vOTU_vs_CRISPR.perfect_matches.tsv"}' \
  host_prediction/crispr/vOTU_vs_CRISPR.perfect_matches.tsv \
  > host_prediction/candidate_pairs/LMVD_candidate_pairs.CRISPR.tsv

# Extract candidate pairs from nucleotide homology matches.
# In this BLAST format, query is the host genome and subject is the vOTU.
awk 'BEGIN{OFS="\t"} {print $2,$1,"nucleotide_homology",$3,$4,$4/$14,"mOTU_vs_vOTU.homology_filtered.tsv"}' \
  host_prediction/nucleotide_homology/mOTU_vs_vOTU.homology_filtered.tsv \
  > host_prediction/candidate_pairs/LMVD_candidate_pairs.homology.tsv

# Extract candidate pairs from tRNA matches.
# The exact vOTU_id and host_id fields depend on the tRNA FASTA header format.
awk 'BEGIN{OFS="\t"} {print $2,$1,"tRNA",$3,$4,"full_length","mOTU_tRNA_vs_vOTU_tRNA.full_length_exact.tsv"}' \
  host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.full_length_exact.tsv \
  > host_prediction/candidate_pairs/LMVD_candidate_pairs.tRNA.tsv

# Combine all candidate-pair tables.
echo -e "vOTU_id\thost_id\tevidence_type\tidentity\talignment_length\tcoverage\tsource_table" \
  > host_prediction/candidate_pairs/LMVD_candidate_virus_host_pairs.tsv

cat \
  host_prediction/candidate_pairs/LMVD_candidate_pairs.CRISPR.tsv \
  host_prediction/candidate_pairs/LMVD_candidate_pairs.homology.tsv \
  host_prediction/candidate_pairs/LMVD_candidate_pairs.tRNA.tsv \
  >> host_prediction/candidate_pairs/LMVD_candidate_virus_host_pairs.tsv
```

The exact `awk` fields depend on the sequence-header format used for CRISPR spacers, MAGs, vOTUs and tRNAs and should be verified before large-scale execution.

## 7. Prepare VirMatcher input files

VirMatcher requires:

1. a viral FASTA file;
2. a directory containing bacterial or archaeal host genomes in FASTA format;
3. a host taxonomy file.

Prepare host genome directory:

```bash
mkdir -p host_prediction/VirMatcher/bac-dir

cp ${MOTU_REP_DIR}/*.fa \
  host_prediction/VirMatcher/bac-dir/
```

If required by the local VirMatcher installation, rename FASTA files from `.fa` to `.fasta`:

```bash
for f in host_prediction/VirMatcher/bac-dir/*.fa
do
  mv "$f" "${f%.fa}.fasta"
done
```

Prepare viral FASTA:

```bash
cp ${VOTU_FASTA} \
  host_prediction/VirMatcher/LMVD_vOTUs.fasta
```

Prepare host taxonomy table:

```bash
# Required format:
# host_name<TAB>taxonomy

cp metadata/LMVD_mOTU_taxonomy_for_VirMatcher.tsv \
  host_prediction/VirMatcher/LMVD_mOTU_taxonomy.tsv
```

Representative taxonomy records:

```text
mOTU_000001    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria
mOTU_000002    d__Bacteria;p__Bacteroidota;c__Bacteroidia
```

## 8. VirMatcher verification of candidate virus-host links

VirMatcher was used to evaluate the candidate virus-host pairs by integrating BLAST, CRISPR, tRNA and WIsH evidence.

```bash
mkdir -p host_prediction/VirMatcher/out

conda activate VirMatcher

export PATH=$PATH:/path/to/WIsH
export PATH=$PATH:/path/to/virmatcher/bin

VirMatcher \
  --virus-fp host_prediction/VirMatcher/LMVD_vOTUs.fasta \
  --bacteria-host-dir host_prediction/VirMatcher/bac-dir \
  --bacteria-taxonomy host_prediction/VirMatcher/LMVD_mOTU_taxonomy.tsv \
  --threads ${THREADS} \
  -o host_prediction/VirMatcher/out
```

Some VirMatcher installations support `--python-aggregator`. Availability and behaviour of this option are version-dependent, and the default aggregator can be used when it completes successfully:

```bash
VirMatcher \
  --virus-fp host_prediction/VirMatcher/LMVD_vOTUs.fasta \
  --bacteria-host-dir host_prediction/VirMatcher/bac-dir \
  --bacteria-taxonomy host_prediction/VirMatcher/LMVD_mOTU_taxonomy.tsv \
  --threads ${THREADS} \
  -o host_prediction/VirMatcher/out \
  --python-aggregator
```

Virus-host links with interpretation score >= 3 were retained. This filtering step was performed from the VirMatcher result table using the column containing the final interpretation score.

Representative filtering logic:

```bash
# Set the result-table name and score-column number according to the local VirMatcher output.
# This structure assumes that the final interpretation score is stored in column 5.
VIRMATCHER_RESULT="host_prediction/VirMatcher/out/VirMatcher_results.tsv"
SCORE_COLUMN=5

awk -v col=${SCORE_COLUMN} 'BEGIN{FS=OFS="\t"} NR==1 || $col>=3 {print $0}' \
  ${VIRMATCHER_RESULT} \
  > host_prediction/VirMatcher/LMVD_vOTU_mOTU_links.score_ge3.tsv
```

The VirMatcher output file name and score-column position may differ among installations. This step should therefore be set after inspecting the local VirMatcher result table.

## 9. Optional standalone WIsH commands

In the main workflow, WIsH evidence was incorporated through VirMatcher. The commands below document the standalone WIsH logic for transparency.

WIsH expects one directory containing host genome FASTA files and one directory containing viral FASTA files. Host models are first built from bacterial or archaeal genomes:

```bash
mkdir -p host_prediction/WIsH/models

WIsH \
  -c build \
  -g host_prediction/VirMatcher/bac-dir \
  -m host_prediction/WIsH/models \
  -t ${THREADS}
```

Prediction is then performed for viral genomes or contigs:

```bash
mkdir -p host_prediction/WIsH/prediction

WIsH \
  -c predict \
  -g host_prediction/WIsH/virus_fasta_dir \
  -m host_prediction/WIsH/models \
  -r host_prediction/WIsH/prediction \
  -b 1 \
  -t ${THREADS}
```

Main WIsH outputs include:

```text
host_prediction/WIsH/prediction/llikelihood.matrix
host_prediction/WIsH/prediction/prediction.list
```

## 10. Expected key outputs

```text
host_prediction/crispr/vOTU_vs_CRISPR.perfect_matches.tsv
host_prediction/nucleotide_homology/mOTU_vs_vOTU.homology_filtered.tsv
host_prediction/trna/mOTU_tRNA_vs_vOTU_tRNA.full_length_exact.tsv
host_prediction/candidate_pairs/LMVD_candidate_virus_host_pairs.tsv
host_prediction/VirMatcher/LMVD_vOTU_mOTU_links.score_ge3.tsv
```

## 11. Scope notes

1. This file documents representative commands and filtering criteria rather than a complete executable workflow for all 4,017 samples.
2. The in situ host database consisted of representative mOTU genomes obtained after MAG quality assessment and dRep dereplication.
3. Candidate virus-host pairs were initially inferred using CRISPR spacer matching, nucleotide sequence homology and tRNA matching.
4. CRISPR spacer matches were retained only when the spacer showed 100% sequence identity and full-length coverage.
5. Nucleotide homology matches were retained only when the alignment length was >= 2,500 bp, identity was >= 90% and the aligned region covered >= 75% of the vOTU sequence.
6. tRNA matches were retained only when the tRNA sequences showed 100% identity and full-length coverage for both query and subject tRNAs.
7. VirMatcher was used as the final verification step, and only virus-host links with interpretation score >= 3 were retained.
8. VirMatcher output file names may differ depending on installation and version, so downstream filtering scripts should be configured after inspecting the generated result table.
