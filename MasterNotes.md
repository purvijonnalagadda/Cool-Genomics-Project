# Gut Phageome Analysis: MDD vs Control (PRJEB76994)

## Project Goal

Identify and compare bacteriophage communities in gut metagenomic samples from individuals with major depressive disorder (MDD) and controls using 6 samples.

---

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

## Controls
fasterq-dump ERR13348288 --split-files
fasterq-dump ERR13348290 --split-files
fasterq-dump ERR13348291 --split-files

## MDD
fasterq-dump ERR13348393 --split-files
fasterq-dump ERR13348396 --split-files
fasterq-dump ERR13348397 --split-files

gzip *.fastq

### Converts SRA data into FASTQ files


# STEP 3: ORGANIZE PROJECT

```bash
mkdir project_fastq
cd project_fastq
mkdir raw trimmed assembly virsorter checkv logs fastqc_out

mv ../*.fastq.gz raw/
```

### Organizes files into directories


# STEP 4: FASTQC (RAW READS)

```bash
fastqc raw/*.fastq.gz -o fastqc_out
```

### Checks sequencing quality and adapter content


# STEP 5: TRIM READS

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  fastp \
    -i raw/${sample}_1.fastq.gz \
    -I raw/${sample}_2.fastq.gz \
    -o trimmed/${sample}_R1.fastq.gz \
    -O trimmed/${sample}_R2.fastq.gz \
    -h logs/${sample}_fastp.html \
    -j logs/${sample}_fastp.json
```

### Removes adapters and low-quality bases


# STEP 6: FASTQC (TRIMMED READS)

```bash
fastqc trimmed/*.fastq.gz -o fastqc_out
```

### Confirms improved read quality


# STEP 7: ASSEMBLY (MEGAHIT)

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  megahit \
    -1 trimmed/${sample}_R1.fastq.gz \
    -2 trimmed/${sample}_R2.fastq.gz \
    -t 8 \
    -o assembly/${sample}
```

### Assembles reads into contigs


# STEP 8: ASSEMBLY STATS

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  seqkit stats -a assembly/${sample}/final.contigs.fa \
    > assembly/${sample}_stats.txt
```

### Generates assembly metrics (N50, GC, contig count)


# STEP 9: SETUP VIRSORTER2

virsorter setup -d ~/db -j 4 --conda-frontend conda

### Downloads viral database


# STEP 10: RUN VIRSORTER2

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  virsorter run \
    -w virsorter/${sample} \
    -i assembly/${sample}/final.contigs.fa \
    --keep-original-seq \
    --include-groups dsDNAphage,ssDNA \
    --min-length 5000
```

### Identifies viral contigs (bacteriophages)


# STEP 11: COUNT VIRAL CONTIGS

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  count=$(grep -c ">" virsorter/${sample}/final-viral-combined.fa)
  echo -e "${sample}\t${count}"
done > viral_counts.tsv
```

### Counts viral contigs per sample


# STEP 12: FILTER CONTIGS ≥ 5 KB

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  seqkit seq -m 5000 \
    virsorter/${sample}/final-viral-combined.fa \
    > virsorter/${sample}/viral_5kb.fa
```

### Keeps high-confidence viral contigs


# STEP 13: SETUP CHECKV

```bash
mkdir checkv_db
cd checkv_db
checkv download_database ./
cd ..
```

### Downloads CheckV database


# STEP 14: RUN CHECKV

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  checkv end_to_end \
    virsorter/${sample}/viral_5kb.fa \
    checkv/${sample} \
    -d checkv_db \
    -t 8
```

### Assesses viral genome quality


# STEP 15: SUMMARIZE CHECKV

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397

```bash
  echo "=== ${sample} ==="
  tail -n +2 checkv/${sample}/quality_summary.tsv | cut -f8 | sort | uniq -c
```

### Summarizes quality categories


# STEP 16: CREATE METADATA

```bash
echo -e "Sample\tGroup
ERR13348288\tControl
ERR13348290\tControl
ERR13348291\tControl
ERR13348393\tMDD
ERR13348396\tMDD
ERR13348397\tMDD" > metadata.tsv
```

### Defines sample groups


# STEP 17: FINAL ANALYSIS

Compare:
Viral contig counts (viral_counts.tsv)
High-quality viral contigs (CheckV)

### Goal:
Identify differences in phage abundance between MDD and controls
