# Functional annotation, plasmid prediction, Ka/Ks analysis and contig-level abundance profiling

This file documents representative command structures for functional annotation and sequence profiling of assembled contigs, viral contigs (VCs), vOTUs and MAGs. It covers ARG annotation, MGE/VF annotation, plasmid prediction, Ka/Ks calculation and contig-level abundance estimation.

The commands are provided to describe the principal tools, parameters and file relationships used in the analysis. They are not intended to constitute an end-to-end executable workflow for all 4,017 samples, and local paths, database locations and sample identifiers should be specified according to the computing environment.

The same annotation logic can be applied to different sequence sets, including total assembled contigs, curated viral contigs, vOTU representative sequences and MAG sequences.

## 1. Representative input variables

```bash
SAMPLE="Sample_001"
THREADS=24

# Representative sequence inputs
CONTIG_FASTA="all_contigs/${SAMPLE}.fa"
VC_FASTA="viral_prediction/${SAMPLE}/curated/${SAMPLE}.VC.keep1_keep2.min5kb.fa"
MAG_FASTA_DIR="MAGs/${SAMPLE}"

# Representative paired-end reads for abundance profiling
READ1="1-Trim/PE/${SAMPLE}_1.fastq.gz"
READ2="1-Trim/PE/${SAMPLE}_2.fastq.gz"
```

Representative output organization:

```text
functional_annotation/
├── ARG_RGI/
├── MGE_mobileOG/
├── VF_VFDB/
├── plasmid_prediction/
├── KaKs/
└── contig_abundance/
```

## 2. ARG annotation using RGI/CARD

ARGs were annotated using RGI against the CARD database. The same command structure can be applied to assembled contigs, VCs or MAG sequences.

Representative command for assembled contigs:

```bash
mkdir -p functional_annotation/ARG_RGI/${SAMPLE}

conda activate /path/to/rgi_env

rgi main \
  --input_sequence ${CONTIG_FASTA} \
  --output_file functional_annotation/ARG_RGI/${SAMPLE}/${SAMPLE}.contigs.rgi \
  --clean \
  -a DIAMOND \
  -n ${THREADS}
```

Representative command for curated viral contigs:

```bash
rgi main \
  --input_sequence ${VC_FASTA} \
  --output_file functional_annotation/ARG_RGI/${SAMPLE}/${SAMPLE}.VCs.rgi \
  --clean \
  -a DIAMOND \
  -n ${THREADS}
```

Representative command for MAGs:

```bash
mkdir -p functional_annotation/ARG_RGI/${SAMPLE}/MAGs

for MAG in ${MAG_FASTA_DIR}/*.fa
do
  MAG_ID=$(basename ${MAG} .fa)

  rgi main \
    --input_sequence ${MAG} \
    --output_file functional_annotation/ARG_RGI/${SAMPLE}/MAGs/${MAG_ID}.rgi \
    --clean \
    -a DIAMOND \
    -n ${THREADS}
done
```

Main outputs usually include:

```text
*.rgi.txt
*.rgi.json
```

## 3. MGE annotation using DIAMOND and mobileOG

MGEs were annotated by aligning nucleotide sequences against the mobileOG protein database using DIAMOND blastx. The best hit for each query was retained and then filtered by sequence identity and coverage.

The input FASTA can correspond to total contigs, ARG-carrying contigs, ARG-neighbourhood sequences, VCs or MAG sequences. `INPUT_FASTA` should be assigned according to the sequence set under analysis.

```bash
mkdir -p functional_annotation/MGE_mobileOG/${SAMPLE}

INPUT_FASTA="all_contigs/${SAMPLE}.fa"
MOBILEOG_DB="/path/to/mobileOG.faa.dmnd"

diamond blastx \
  --db ${MOBILEOG_DB} \
  -q ${INPUT_FASTA} \
  -o functional_annotation/MGE_mobileOG/${SAMPLE}/${SAMPLE}.mobileOG.tsv \
  --sensitive \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
  -e 1e-5 \
  -p 4
```

Sort hits by query and bit score:

```bash
sort -k1,1 -k12,12nr \
  functional_annotation/MGE_mobileOG/${SAMPLE}/${SAMPLE}.mobileOG.tsv \
  > functional_annotation/MGE_mobileOG/${SAMPLE}/${SAMPLE}.mobileOG.sorted.tsv
```

Retain the best hit per query using the best-hit script:

```bash
python /path/to/get_BestHit_contig.py \
  functional_annotation/MGE_mobileOG/${SAMPLE}/${SAMPLE}.mobileOG.sorted.tsv \
  functional_annotation/MGE_mobileOG/${SAMPLE}/${SAMPLE}.mobileOG.best.tsv
```

Filter best hits using identity >= 80% and subject coverage >= 80%:

```bash
awk '$3>=80 && $4/$14>=0.8 {print $0}' \
  functional_annotation/MGE_mobileOG/${SAMPLE}/${SAMPLE}.mobileOG.best.tsv \
  > functional_annotation/MGE_mobileOG/${SAMPLE}/${SAMPLE}.mobileOG.filtered.tsv
```

## 4. VF annotation using DIAMOND and VFDB

Virulence factors were annotated using the same DIAMOND blastx workflow as the MGE annotation, replacing the database with VFDB.

```bash
mkdir -p functional_annotation/VF_VFDB/${SAMPLE}

INPUT_FASTA="all_contigs/${SAMPLE}.fa"
VFDB_DB="/path/to/VFDB.dmnd"

diamond blastx \
  --db ${VFDB_DB} \
  -q ${INPUT_FASTA} \
  -o functional_annotation/VF_VFDB/${SAMPLE}/${SAMPLE}.VFDB.tsv \
  --sensitive \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
  -e 1e-5 \
  -p 4

sort -k1,1 -k12,12nr \
  functional_annotation/VF_VFDB/${SAMPLE}/${SAMPLE}.VFDB.tsv \
  > functional_annotation/VF_VFDB/${SAMPLE}/${SAMPLE}.VFDB.sorted.tsv

python /path/to/get_BestHit_contig.py \
  functional_annotation/VF_VFDB/${SAMPLE}/${SAMPLE}.VFDB.sorted.tsv \
  functional_annotation/VF_VFDB/${SAMPLE}/${SAMPLE}.VFDB.best.tsv

awk '$3>=80 && $4/$14>=0.8 {print $0}' \
  functional_annotation/VF_VFDB/${SAMPLE}/${SAMPLE}.VFDB.best.tsv \
  > functional_annotation/VF_VFDB/${SAMPLE}/${SAMPLE}.VFDB.filtered.tsv
```

The VFDB database path should be specified according to the local DIAMOND database configuration.

## 5. Plasmid prediction for contigs

Plasmid-associated contigs were predicted using multiple plasmid prediction tools. The same input FASTA can be total assembled contigs or a subset of contigs relevant to downstream analysis.

### 5.1 PlasFlow

```bash
mkdir -p functional_annotation/plasmid_prediction/PlasFlow/${SAMPLE}

conda activate plasflow

PlasFlow.py \
  --input ${CONTIG_FASTA} \
  --output functional_annotation/plasmid_prediction/PlasFlow/${SAMPLE}/${SAMPLE}.PlasFlow.tsv
```

Extract contig-level PlasFlow labels:

```bash
cut -f3,6 \
  functional_annotation/plasmid_prediction/PlasFlow/${SAMPLE}/${SAMPLE}.PlasFlow.tsv \
  | sed '/label/d' \
  > functional_annotation/plasmid_prediction/PlasFlow/${SAMPLE}/${SAMPLE}.PlasFlow.map.tsv
```

### 5.2 PlasClass

```bash
mkdir -p functional_annotation/plasmid_prediction/PlasClass/${SAMPLE}

conda activate plasclass

classify_fasta.py \
  -f ${CONTIG_FASTA} \
  -o functional_annotation/plasmid_prediction/PlasClass/${SAMPLE}/${SAMPLE}.PlasClass.tsv \
  -p 16
```

### 5.3 PlasmidHunter

```bash
mkdir -p functional_annotation/plasmid_prediction/PlasmidHunter/${SAMPLE}

conda activate plasmidhunter

plasmidhunter \
  -i ${CONTIG_FASTA} \
  -o functional_annotation/plasmid_prediction/PlasmidHunter/${SAMPLE}/${SAMPLE}.PlasmidHunter.tsv \
  -c 16
```

### 5.4 geNomad

```bash
mkdir -p functional_annotation/plasmid_prediction/geNomad/${SAMPLE}

conda activate /path/to/genomad_env

genomad end-to-end \
  --cleanup \
  -t 16 \
  ${CONTIG_FASTA} \
  functional_annotation/plasmid_prediction/geNomad/${SAMPLE}/${SAMPLE}.geNomad \
  /path/to/genomad_db
```

### 5.5 PLASMe

```bash
mkdir -p functional_annotation/plasmid_prediction/PLASMe/${SAMPLE}

conda activate plasme

cd /path/to/PLASMe

python PLASMe.py \
  /absolute/path/to/${CONTIG_FASTA} \
  /absolute/path/to/functional_annotation/plasmid_prediction/PLASMe/${SAMPLE}/${SAMPLE}.PLASMe.fa \
  -t 16
```

### 5.6 Merging plasmid prediction outputs

Outputs from PlasFlow, PlasClass, PlasmidHunter, geNomad and PLASMe were parsed according to their tool-specific formats and merged into a contig-level plasmid prediction table.

Representative merged table format:

```text
contig_id    PlasFlow_label    PlasClass_score    PlasmidHunter_label    geNomad_label    PLASMe_label    final_plasmid_assignment
```

This merging step was performed by table parsing and formatting of tool-specific outputs rather than by a single standardized standalone script.

## 6. Ka/Ks analysis of ARG sequences

For each target ARG subtype, ARG sequences were compared with the corresponding CARD reference sequence. Protein sequences were aligned using MAFFT, codon-based nucleotide alignments were generated using PAL2NAL v14, and Ka/Ks ratios were calculated using KaKs_Calculator.

### 6.1 Split multi-FASTA files

Representative command structure for CDS sequences:

```bash
mkdir -p functional_annotation/KaKs/split_CDS
cd functional_annotation/KaKs/split_CDS

cp /path/to/vARG.fna ./

csplit vARG.fna '/>/' -n2 -s {*} -f gene -b "%1d.fa"
rm -f gene0.fa

ls gene*.fa | sed 's/.fa//' | while read i
do
  NAME=$(head -n 1 ${i}.fa | sed 's/>//')
  mv ${i}.fa ${NAME}.fna
done
```

Representative command structure for protein sequences:

```bash
mkdir -p functional_annotation/KaKs/split_protein
cd functional_annotation/KaKs/split_protein

cp /path/to/vARG.faa ./

csplit vARG.faa '/>/' -n2 -s {*} -f gene -b "%1d.fa"
rm -f gene0.fa

ls gene*.fa | sed 's/.fa//' | while read i
do
  NAME=$(head -n 1 ${i}.fa | sed 's/>//')
  mv ${i}.fa ${NAME}.faa
done
```

### 6.2 Pair each ARG sequence with its CARD reference

For each ARG sequence, concatenate the corresponding CARD reference protein sequence and the query protein sequence:

```bash
ARG_ID="test"
CARD_ID="3002837-342"

mkdir -p functional_annotation/KaKs/${ARG_ID}

cat CARD-faa/${CARD_ID}.faa ${ARG_ID}.faa \
  > functional_annotation/KaKs/${ARG_ID}/${ARG_ID}.paired.faa

cat CARD-fna/${CARD_ID}.fna ${ARG_ID}.fna \
  > functional_annotation/KaKs/${ARG_ID}/${ARG_ID}.paired.fna
```

### 6.3 Protein alignment using MAFFT

```bash
conda activate ribotree

mafft --auto \
  functional_annotation/KaKs/${ARG_ID}/${ARG_ID}.paired.faa \
  > functional_annotation/KaKs/${ARG_ID}/${ARG_ID}.aligned.faa
```

### 6.4 Codon alignment using PAL2NAL

```bash
perl /path/to/pal2nal.v14/pal2nal.pl \
  functional_annotation/KaKs/${ARG_ID}/${ARG_ID}.aligned.faa \
  functional_annotation/KaKs/${ARG_ID}/${ARG_ID}.paired.fna \
  -output fasta \
  > functional_annotation/KaKs/${ARG_ID}/${ARG_ID}.codon_alignment.fasta
```

### 6.5 Convert FASTA alignment to AXT and calculate Ka/Ks

```bash
cd functional_annotation/KaKs/${ARG_ID}

perl /path/to/convert_fasta_to_axt.pl \
  ${ARG_ID}.codon_alignment.fasta

KaKs_Calculator \
  -i ${ARG_ID}.codon_alignment.axt \
  -o ${ARG_ID}.kaks.tsv
```

The AXT file name depends on the local FASTA-to-AXT conversion script; the `KaKs_Calculator` input should be set accordingly.

## 7. Contig and VC abundance profiling using CoverM

Coverage of assembled contigs and VCs was calculated by mapping reads back to the sequences assembled from the corresponding metagenome. The `trimmed_mean` metric and the following read-alignment filters were used:

```text
--min-read-aligned-percent 75
--min-read-aligned-length 75
--min-read-percent-identity 80
```

Representative command for all assembled contigs:

```bash
mkdir -p functional_annotation/contig_abundance/all_contigs

coverm contig \
  -m trimmed_mean \
  -t 16 \
  --coupled ${READ1} ${READ2} \
  --reference ${CONTIG_FASTA} \
  -o functional_annotation/contig_abundance/all_contigs/${SAMPLE}.contig_coverage.tsv \
  --min-read-aligned-percent 75 \
  --min-read-aligned-length 75 \
  --min-read-percent-identity 80
```

Representative command for curated viral contigs from the same sample:

```bash
mkdir -p functional_annotation/contig_abundance/VCs

coverm contig \
  -m trimmed_mean \
  -t 16 \
  --coupled ${READ1} ${READ2} \
  --reference ${VC_FASTA} \
  -o functional_annotation/contig_abundance/VCs/${SAMPLE}.VC_coverage.tsv \
  --min-read-aligned-percent 75 \
  --min-read-aligned-length 75 \
  --min-read-percent-identity 80
```

The cell-normalized abundance of each contig or VC was calculated by dividing its coverage by the estimated microbial cell number in the corresponding sample. For vOTUs, abundance was calculated as the mean cell-normalized abundance of all VCs assigned to the same vOTU cluster.

## 8. Scope notes

1. This file documents representative commands and filtering criteria rather than a complete executable workflow for all 4,017 samples or all sequence sets.
2. ARG annotation was performed using RGI/CARD.
3. MGE and VF annotations used the same DIAMOND blastx best-hit workflow, with mobileOG and VFDB used as the corresponding databases.
4. The `get_BestHit_contig.py` script is used to retain the best hit per query after sorting BLAST/DIAMOND output by query ID and bit score.
5. Plasmid prediction was performed using multiple tools, and the final contig-level table was generated by parsing and merging tool-specific outputs.
6. Ka/Ks analysis was performed by comparing each target ARG sequence with the corresponding CARD reference sequence, rather than by all-vs-all comparison among viral, plasmid and chromosomal ARGs.
7. Contig and VC abundance were estimated only by mapping reads back to sequences assembled from the corresponding metagenome.
