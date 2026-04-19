# Gut Phageome Analysis: MDD vs Control (PRJEB76994)

## Project Goal

Identify and compare bacteriophage communities in gut metagenomic samples from individuals with major depressive disorder (MDD) and controls using 6 samples.

---

## Full Workflow (All Steps)

```bash
##############################
# STEP 1: LOAD ENVIRONMENT
##############################

cd ~

module load anaconda3
conda activate phage-env

# Purpose:
# Loads all required bioinformatics tools (fastqc, fastp, megahit, virsorter, checkv)


##############################
# STEP 2: DOWNLOAD FASTQ FILES
##############################

module load sra-toolkit

# Controls
fasterq-dump ERR13348288 --split-files
fasterq-dump ERR13348290 --split-files
fasterq-dump ERR13348291 --split-files

# MDD
fasterq-dump ERR13348393 --split-files
fasterq-dump ERR13348396 --split-files
fasterq-dump ERR13348397 --split-files

gzip *.fastq

# Purpose:
# Converts SRA data into FASTQ files for analysis


##############################
# STEP 3: ORGANIZE PROJECT
##############################

mkdir project_fastq
cd project_fastq

mkdir raw trimmed assembly virsorter checkv logs fastqc_out

mv ../*.fastq.gz raw/

# Purpose:
# Organizes files into structured directories


##############################
# STEP 4: FASTQC (RAW READS)
##############################

fastqc raw/*.fastq.gz -o fastqc_out

# Purpose:
# Evaluates sequencing quality, GC content, and adapter contamination


##############################
# STEP 5: TRIM READS (FASTP)
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  fastp \
    -i raw/${sample}_1.fastq.gz \
    -I raw/${sample}_2.fastq.gz \
    -o trimmed/${sample}_R1.fastq.gz \
    -O trimmed/${sample}_R2.fastq.gz \
    -h logs/${sample}_fastp.html \
    -j logs/${sample}_fastp.json
done

# Purpose:
# Removes adapters and low-quality bases to improve downstream accuracy


##############################
# STEP 6: FASTQC (TRIMMED READS)
##############################

fastqc trimmed/*.fastq.gz -o fastqc_out

# Purpose:
# Confirms improved read quality after trimming


##############################
# STEP 7: ASSEMBLY (MEGAHIT)
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  megahit \
    -1 trimmed/${sample}_R1.fastq.gz \
    -2 trimmed/${sample}_R2.fastq.gz \
    -t 8 \
    -o assembly/${sample}
done

# Purpose:
# Assembles short reads into contigs for viral identification


##############################
# STEP 8: ASSEMBLY STATS (SEQKIT)
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  seqkit stats -a assembly/${sample}/final.contigs.fa \
    > assembly/${sample}_stats.txt
done

# Purpose:
# Generates assembly metrics (N50, contig count, GC content)


##############################
# STEP 9: SETUP VIRSORTER2
##############################

virsorter setup -d ~/db -j 4 --conda-frontend conda

# Purpose:
# Downloads viral reference database for detection


##############################
# STEP 10: RUN VIRSORTER2
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  virsorter run \
    -w virsorter/${sample} \
    -i assembly/${sample}/final.contigs.fa \
    --keep-original-seq \
    --include-groups dsDNAphage,ssDNA \
    --min-length 5000
done

# Purpose:
# Identifies viral (bacteriophage) contigs from assemblies


##############################
# STEP 11: COUNT VIRAL CONTIGS
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  count=$(grep -c ">" virsorter/${sample}/final-viral-combined.fa)
  echo -e "${sample}\t${count}"
done > viral_counts.tsv

# Purpose:
# Counts number of predicted viral contigs per sample


##############################
# STEP 12: FILTER ≥5KB CONTIGS
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  seqkit seq -m 5000 \
    virsorter/${sample}/final-viral-combined.fa \
    > virsorter/${sample}/viral_5kb.fa
done

# Purpose:
# Removes short contigs to improve confidence in viral predictions


##############################
# STEP 13: SETUP CHECKV
##############################

mkdir checkv_db
cd checkv_db
checkv download_database ./
cd ..

# Purpose:
# Downloads database for viral genome quality assessment


##############################
# STEP 14: RUN CHECKV
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  checkv end_to_end \
    virsorter/${sample}/viral_5kb.fa \
    checkv/${sample} \
    -d checkv_db \
    -t 8
done

# Purpose:
# Evaluates completeness and quality of viral genomes


##############################
# STEP 15: SUMMARIZE CHECKV RESULTS
##############################

for sample in ERR13348288 ERR13348290 ERR13348291 ERR13348393 ERR13348396 ERR13348397
do
  echo "=== ${sample} ==="
  tail -n +2 checkv/${sample}/quality_summary.tsv | cut -f8 | sort | uniq -c
done

# Purpose:
# Summarizes viral quality categories (complete, high, medium, low)


##############################
# STEP 16: CREATE METADATA FILE
##############################

echo -e "Sample\tGroup
ERR13348288\tControl
ERR13348290\tControl
ERR13348291\tControl
ERR13348393\tMDD
ERR13348396\tMDD
ERR13348397\tMDD" > metadata.tsv

# Purpose:
# Defines sample grouping for comparison


##############################
# STEP 17: FINAL ANALYSIS GOAL
##############################

# Compare:
# - Total viral contigs (viral_counts.tsv)
# - High-quality viral contigs (CheckV output)

# Expected Output:
# Table comparing MDD vs Control viral abundance

# Purpose:
# Enables biological interpretation of phage differences between groups
