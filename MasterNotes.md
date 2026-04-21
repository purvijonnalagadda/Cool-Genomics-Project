# Gut Phageome Analysis: MDD vs Control (PRJEB76994)

## Project Goal

Identify and compare bacteriophage communities in gut metagenomic samples from individuals with major depressive disorder (MDD) and controls using 6 samples.

## Project Directory Structure

```bash
/home/NETID/project_fastq/
│
├── raw/                   # original FASTQ files from SRA
├── trimmed/               # Trimmomatic output
├── fastqc_raw/            # FastQC reports for raw reads
├── fastqc_trimmed/        # FastQC reports for trimmed reads
├── assembly/              # MEGAHIT outputs
├── virsorter/             # VirSorter2 outputs
├── checkv/                # CheckV outputs
├── metadata/              # sample metadata
├── counts/                # summary tables
├── logs/                  # SLURM .out / .err files
├── slurm_scripts/         # SLURM scripts
├── db/                    # VirSorter2 database
├── checkv_db/             # CheckV database
└── sample_names.txt       # list of sample accessions
```
## Full Workflow


# STEP 1: LOAD ENVIRONMENT

```bash
cd ~
module load anaconda3
conda activate phage-env
```

### Loads required bioinformatics tools


# STEP 2: DOWNLOAD FASTQ FILES

```bash
module load sra-toolkit
```

```bash
## Controls
fasterq-dump ERR13348292 --split-files
fasterq-dump ERR13348320 --split-files

## MDD
fasterq-dump ERR13348298 --split-files
fasterq-dump ERR13348304 --split-files

gzip *.fastq
```

### Converts SRA data into FASTQ files


# STEP 3: ORGANIZE PROJECT

```bash
mkdir project_fastq
cd project_fastq
mkdir raw trimmed fastqc_raw fastqc_trimmed assembly virsorter checkv metadata counts logs slurm_scripts db checkv_db

mv ../*.fastq.gz raw/

echo -e "ERR13348292
ERR13348320
ERR13348298
ERR13348304" > sample_names.txt
```

### Organizes files into directories and creates a sample list


# STEP 4: LOCATE TRIMMOMATIC ADAPTER FILE

```bash
find $CONDA_PREFIX -name "TruSeq3-PE.fa"

ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"
```

### Set adapter path and define standard adapter file for trimming


# STEP 5: FASTQC (RAW READS)

```bash
fastqc raw/*.fastq.gz -o fastqc_raw
```
### Checks sequencing quality and adapter content before trimming

# STEP 6: TRIM READS WITH TRIMMOMATIC

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  trimmomatic PE -threads 8 -phred33 \
    raw/${sample}_1.fastq.gz raw/${sample}_2.fastq.gz \
    trimmed/${sample}_R1_paired.fastq.gz trimmed/${sample}_R1_unpaired.fastq.gz \
    trimmed/${sample}_R2_paired.fastq.gz trimmed/${sample}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    SLIDINGWINDOW:4:20 \
    MINLEN:50
```

### Removes adapter sequences and trims low-quality bases. Only paired reads will be used for assembly.


# STEP 7: FASTQC (TRIMMED READS)

```bash
fastqc trimmed/*.fastq.gz -o fastqc_trimmed
```

### Confirms improved read quality


# STEP 8: ASSEMBLY (MEGAHIT)

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  megahit \
    -1 trimmed/${sample}_R1.fastq.gz \
    -2 trimmed/${sample}_R2.fastq.gz \
    -t 8 \
    -o assembly/${sample}
```

### Assembles reads into contigs


# STEP 9: ASSEMBLY STATS

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  seqkit stats -a assembly/${sample}/final.contigs.fa \
    > assembly/${sample}_stats.txt
```

### Generates assembly metrics (N50, GC, contig count)


# STEP 10: SETUP VIRSORTER2

```bash
virsorter setup -d ~/db -j 4 --conda-frontend conda
```

### Downloads VirSorter2 database


# STEP 11: RUN VIRSORTER2

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  virsorter run \
    -w virsorter/${sample} \
    -i assembly/${sample}/final.contigs.fa \
    --keep-original-seq \
    --include-groups dsDNAphage,ssDNA \
    --min-length 5000
```

### Identifies viral contigs


# STEP 12: COUNT VIRAL CONTIGS

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  count=$(grep -c ">" virsorter/${sample}/final-viral-combined.fa)
  echo -e "${sample}\t${count}"
done > counts/viral_counts.tsv
```

### Counts viral contigs per sample


# STEP 13: FILTER CONTIGS ≥ 5 KB

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  seqkit seq -m 5000 \
    virsorter/${sample}/final-viral-combined.fa \
    > virsorter/${sample}/viral_5kb.fa
```

### Keeps high-confidence viral contigs


# STEP 14: SETUP CHECKV

```bash
cd checkv_db
checkv download_database ./
cd ..
```

### Downloads CheckV database


# STEP 15: RUN CHECKV

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  checkv end_to_end \
    virsorter/${sample}/viral_5kb.fa \
    checkv/${sample} \
    -d checkv_db \
    -t 8
```

### Assesses viral genome quality


# STEP 16: SUMMARIZE CHECKV

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  echo "=== ${sample} ==="
  tail -n +2 checkv/${sample}/quality_summary.tsv | cut -f8 | sort | uniq -c
```

### Summarizes viral quality categories


# STEP 17: CREATE METADATA

```bash
echo -e "Sample\tGroup
ERR13348292\tControl
ERR13348320\tControl
ERR13348298\tMDD
ERR13348304\tMDD" > metadata/metadata.tsv
```

### Defines sample groups


# STEP 18: FINAL ANALYSIS

Compare:
Viral contig counts (counts/viral_counts.tsv)
High-quality viral contigs (CheckV/)

### Goal:
Identify differences in phage abundance between MDD and controls
