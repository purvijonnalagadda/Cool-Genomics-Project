# Gut Phageome Analysis: MDD vs Control (PRJEB76994)

## Project Goal

Identify and compare bacteriophage communities in gut metagenomic samples from individuals with major depressive disorder (MDD) and controls using 6 samples.

## Project Directory Structure

```bash
/home/NETID/FinalProject/
│
├── raw/                   # original FASTQ files from SRA
├── trimmed/               # Trimmomatic output
├── fastqc_raw/            # FastQC reports for raw reads
├── fastqc_trimmed/        # FastQC reports for trimmed reads
├── assembly/              # MEGAHIT outputs
├── virsorter/             # VirSorter2 outputs
├── checkv/                # CheckV outputs
├── Bowtie-Samtooms        # Bowtie-Samtools outputs
└── CoverM                 # CoverM outputs

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
nano slurm_scripts/02_qc_trim.slurm
```
```bash
#!/bin/bash
#SBATCH --job-name=projecttrims.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=user@email.com
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task1 --time=3:00:00

shopt -s expand_aliases
module load trimmomatic

adapters=/home/user/fastq_files_project/TruSeq-PE.fa
input_R1=/home/user/fastq_files_project/raw/<SAMPLE>_1.fastqc.gz
input_R2=/home/user/fastq_files_project/raw/<SAMPLE>.fastqc.gz
output_R1_PE=/home/user/fastq_files_project/trimmed/<SAMPLE>_1_trmPE.fq.gz
output_R1_SE=/home/user/fastq_files_project/trimmed/<SAMPLE>_1_trmPE.fq.gz
output_R2_PE=/home/user/fastq_files_project/trimmed/<SAMPLE_2_trmPE.fq.gz
output_R2_SE=/home/user/fastq_files_project/trimmed/<SAMPLE>_2_trmPE.fq.gz

trimmomatic PE \
$input_R1 \ 
$input_R2 \
$output_R1_PE $output_R1_SE \
$output_R2_PE $output_R2_SE \

LEADING:10
TRAILING:10
ILLUMINACLIP:$adapters:2:30:10 \
SLIDINGWINDOW:4:20 \
MINLEN:50 \

fastqc \
  "${WORKDIR}/trimmed/${SAMPLE}_R1_trmPE.fq.gz" \
  "${WORKDIR}/trimmed/${SAMPLE}_R2_trmPE.fq.gz" \
  -o ${WORKDIR}/fastqc_trimmed
```
Run script:
```bash
sbatch slurm_scripts/02_qc_trim.slurm
```

### Removes adapter sequences and trims low-quality bases. Only paired reads will be used for assembly. Runs FastQC (raw), Trimmomatic trimming, and FastQC (trimmed) via SLURM


# STEP 7: ASSEMBLY (MEGAHIT)

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash

#!/bin/bash
#SBATCH --job-name=megahit_sample298
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/home/sja111/FinalProject/Assembly/logs/megahit_%j.out
#SBATCH --error=/home/sja111/FinalProject/Assembly/megahit_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sja111@georgetown.edu

module load mamba/
source $(mamba info --base)/etc/profile.d/conda.sh
conda activate megahit-env

READ1=/home/sja111/FinalProject/Readcleaning/trimmed/<SAMPLE>_1_trmPE.fq.gz
READ2=/home/sja111/FinalProject/Readcleaning/trimmed/<SAMPLE>_2_trmPE.fq.gz
OUTDIR=/home/sja111/FinalProject/Assembly/Depression/<SAMPLE>_megahit_out

megahit \
  -1 ${READ1} \
  -2 ${READ2} \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${OUTDIR}

echo "Done. Contigs should be in ${OUTDIR}/final.contigs.fa"
```

### Assembles reads into contigs


# STEP 8: ASSEMBLY STATS

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  seqkit stats -a assembly/${sample}/final.contigs.fa \
    > assembly/${sample}_stats.txt
```

### Generates assembly metrics (N50, GC, contig count)


# STEP 9: SETUP VIRSORTER2

```bash
virsorter setup -d ~/db -j 4 --conda-frontend conda
```

### Downloads VirSorter2 database


# STEP 10: RUN VIRSORTER2

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash

#!/bin/bash
#SBATCH --job-name=virsorter292
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sja111@georgetown.edu
#SBATCH --output=/home/sja111/FinalProject/Virsorter/logs/virsorter.j%.out
#SBATCH --error=/home/sja111/FinalProject/Virsorter/logs/virsorter.j%.err

# ==== Load mamba ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where you had VirSorter2 installed
mamba activate vs2-env

# ==== Set paths and filenames ====
#set up directories
INDIR=/home/sja111/FinalProject/Assembly/contigs         #directory where input will come from
OUTROOT=/home/sja111/FinalProject/Virsorter/     #directory output will go
mkdir -p "${OUTROOT}"                                    #new directory to be created for output files

SAMPLE_ID=292                                            #just the basic sample name (sample2 ?)
INPUT="${INDIR}/final.contigs<SAMPLE#>.fa"                     #contig file name/location
OUTDIR="${OUTROOT}/vs2-${SAMPLE_ID}"                     #where you’ll find the output files
mkdir -p "${OUTDIR}"


# ==== Run virsorter2 with >5kb cutoff and DNA virus categories 
echo "Running VirSorter2 on ${INPUT}"
virsorter run \
  -w "${OUTDIR}" \
  -i "${INPUT}" \
  --keep-original-seq \
  --include-groups dsDNAphage,NCLDV,ssDNA \
  --min-length 5000
```

### Identifies viral contigs


# STEP 11: COUNT VIRAL CONTIGS

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  count=$(grep -c ">" virsorter/${sample}/final-viral-combined.fa)
  echo -e "${sample}\t${count}"
done > counts/viral_counts.tsv
```

### Counts viral contigs per sample


# STEP 12: FILTER CONTIGS ≥ 5 KB

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  seqkit seq -m 5000 \
    virsorter/${sample}/final-viral-combined.fa \
    > virsorter/${sample}/viral_5kb.fa
```

### Keeps high-confidence viral contigs


# STEP 13: COMBINE VIRAL CONTIGS ACROSS SAMPLES

```bash
mkdir -p combined_votus

> combined_votus/all_viral_contigs.fa

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304
do
  awk -v sample="$sample" '
    /^>/ {print ">" sample "_" substr($0,2)}
    !/^>/ {print}
  ' virsorter/${sample}/viral_5kb.fa >> combined_votus/all_viral_contigs.fa
done

grep -c ">" combined_votus/all_viral_contigs.fa
```

# STEP 14: STEP 14: CLUSTER VIRAL CONTIGS INTO vOTUs WITH vCLUST

```bash
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

mamba create -y -n votu-env -c bioconda -c conda-forge vclust
mamba activate votu-env
```

```bash
cd combined_votus

vclust prefilter \
  -i all_viral_contigs.fa \
  -o fltr.txt

vclust align \
  -i all_viral_contigs.fa \
  -o ani.tsv \
  --filter fltr.txt

vclust cluster \
  -i ani.tsv \
  -o clusters.tsv \
  --ids ani.ids.tsv \
  --metric ani \
  --ani 0.95 \
  --out-repr
```

```bash
tail -n +2 clusters.tsv | awk '{print $2}' | sort -u > votu_seeds.txt

mamba activate megahit-env

seqkit grep -f votu_seeds.txt all_viral_contigs.fa > votus_final.fna

grep -c ">" votus_final.fna
```

# STEP 15: BUILD BOWTIE2 INDEX FROM vOTUs

```bash
cd /home/NETID/project_fastq

mkdir -p bowtie2

module load bowtie2

bowtie2-build combined_votus/votus_final.fna bowtie2/votu_index
```

# STEP 16: MAP READS BACK TO vOTUs WITH BOWTIE2 AND SAMTOOLS

```bash
#!/bin/bash
#SBATCH --job name=bowtie2_<SAMPLE>
#SBATCH --output=/home/sja111/FinalProject/Bowtie_Samtools/logs/bowtie-%j.out
#SBATCH --error=/home/sja111/FinalProject/Bowtie_Samtools/logs/bowtie-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sja111@georgetown.edu
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

# ---------SET UP----------
SAMPLE="<SAMPLE>"
INDEX="/home/sja111/FinalProject/Bowtie_Samtools/votu_index"
OUTPUTDIR="/home/sja111/FinalProject/Bowtie_Samtools/${SAMPLE}"

# --------- LOAD MODULES ----------
module purge
module load bowtie2/2.5.4

# --------- RUN BOWTIE2 AND PIPE TO SAMTOOLS ----------

# First make output and log directories; move into OUTPUTDIR
mkdir -p "${OUTPUTDIR}"
cd "${OUTPUTDIR}"
mkdir -p logs

echo "Running bowtie2 on sample ${SAMPLE}"

bowtie2 -p 8 -x "${INDEX}" -1 "/home/sja111/FinalProject/Readcleaning/trimmed/<SAMPLE>_1_trmPE.fq.gz" -2 "/home/sja111/FinalProject/Readcleaning/trimmed/<SAMPLE>_2_trmPE.fq.gz" \
| samtools view -bS - > "${SAMPLE}.bam"

echo "Finished running bowtie2 and performing compression"

#---------sort and index files
echo "Sorting"
samtools sort "${SAMPLE}.bam" > "${SAMPLE}_sorted.bam"

echo "Indexing"
samtools index "${SAMPLE}_sorted.bam"

echo "Finished ${SAMPLE}"
```

# STEP 17: ESTIMATE vOTU ABUNDANCE WITH COVERM

```bash
#!/bin/bash
#SBATCH --job-name=coverm_votus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sja111@georgetown.edu
#SBATCH --output=/home/sja111/FinalProject/CoverM/logs/coverm-j%.out
#SBATCH --error=/home/sja111/FinalProject/CoverM/logs/coverm-%j.err

module purge
module load coverm


coverm contig \
  --bam-files \
    /home/sja111/FinalProject/CoverM/sample292_sorted.bam \
    /home/sja111/FinalProject/CoverM/sample320_sorted.bam \
    /home/sja111/FinalProject/CoverM/sample298_sorted.bam \
    /home/sja111/FinalProject/CoverM/sample304_sorted.bam \
  --methods mean rpkm tpm \
  --min-covered-fraction 0.1 \
  --output-file /home/sja111/FinalProject/CoverM/coverm_output/votus_coverm.tsv \
  --output-format dense \
  --threads 8
```

# STEP 18: SETUP CHECKV

```bash
cd checkv_db
checkv download_database ./
cd ..
```

### Downloads CheckV database


# STEP 19: RUN CHECKV

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
#!/bin/bash
#SBATCH --job-name=checkv
#SBATCH --output=/home/sja111/FinalProject/CheckV/logs/checkv-%j.out
#SBATCH --error=/home/sja111/FinalProject/CheckV/logs/checkv-%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sja111@georgetown.edu


# ==== Load checkv program module (students: no need to change) ====

module load checkv


# ==== Set variables, paths, and filenames (students: edit this block!) ====

CHECKVDB="/home/sja111/FinalProject/CheckV/checkv-db-v1.5"

SAMPLE_ID="vOTUs"
INPUT="/home/sja111/FinalProject/Virsorter/votus_final.fna"
OUTDIR="/home/sja111/FinalProject/CheckV/${SAMPLE_ID}"

mkdir -p "${OUTDIR}"


# ==== run checkv (students: no need to change. The second line is the command) ====
echo "Running CheckV on ${INPUT}"
checkv end_to_end "${INPUT}" "${OUTDIR}" -d "${CHECKVDB}" -t ${SLURM_CPUS_PER_TASK}

```

### Assesses viral genome quality


# STEP 20: SUMMARIZE CHECKV

for sample in ERR13348292 ERR13348320 ERR13348298 ERR13348304

```bash
  echo "=== ${sample} ==="
  tail -n +2 checkv/${sample}/quality_summary.tsv | cut -f8 | sort | uniq -c
```

### Summarizes viral quality categories


# STEP 21: CREATE METADATA

Convert .tsv to .csv 

Create vOTU Heatmap
```R
# ==== EDIT THESE THREE LINES ====
filename <- "sja111_votus.csv"  # Excel file name
tpm_threshold <- 10                                  # keep vOTUs with max TPM > this
heatmap_colors <- c("#440154", "#31688e", "#35b779", "#fde725")

# =================================

library(readxl)
library(pheatmap)

# 1. Read the Excel file
cov <- read.csv(filename)

# 2. Keep Contig and TPM columns only
tpm_cols <- grepl("TPM$", names(cov))
cov_tpm <- cov[ , c("Contig", names(cov)[tpm_cols])]

# Remove S1k141_26921_full
cov_tpm <- subset(cov_tpm, Contig != "S1k141_26921||full")

# 3. Optional filter: drop very low-abundance vOTUs
cov_tpm$max_tpm <- apply(cov_tpm[ , -1], 1, max, na.rm = TRUE)
cov_tpm <- subset(cov_tpm, max_tpm > tpm_threshold)
cov_tpm$max_tpm <- NULL

# If everything got filtered out (threshold too high), warn and stop
if (nrow(cov_tpm) == 0) {
  stop("No vOTUs passed the TPM threshold. Try lowering tpm_threshold.")
}

# 4. Make matrix for heatmap (rows = vOTUs, cols = samples)
mat <- as.matrix(cov_tpm[ , -1])
rownames(mat) <- cov_tpm$Contig

# 5. Log-transform for nicer color scaling
mat_log <- log10(mat + 1)

# 6. Draw heat map
pheatmap(mat_log,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(heatmap_colors)(100),
         fontsize_row = 8,
         fontsize_col = 10,
         main = "vOTU relative abundance (log10 TPM + 1)")
```

Determine sample richness and Shannon Diversity

```R
## Install once if needed:
install.packages(c("vegan", "pheatmap", "RColorBrewer"))

## Each session:
library(vegan)
library(pheatmap)
library(RColorBrewer)

## 1. Read your exact file
abund_raw <- read.csv(
  "sja111_votus.csv",
  header = TRUE,
  check.names = FALSE
)

## 2. First column is vOTU / Contig ID. Make that the rownames.
rownames(abund_raw) <- abund_raw$Contig

## 3. Drop the Contig column so we have a numeric matrix only
abund <- abund_raw[, -1]

## 4. Quick sanity check
dim(abund)
head(rownames(abund))   # vOTU IDs
colnames(abund)         # sample1_jbs172_sorted TPM, etc.

##per sample richness and diversity
## Transpose: now rows = samples, columns = vOTUs
abund_t <- t(abund)

## Optional: drop samples with total abundance = 0
abund_t <- abund_t[rowSums(abund_t) > 0, , drop = FALSE]

## 1. Species richness per sample (number of vOTUs with >0)
richness <- specnumber(abund_t)   # named vector, names are sample IDs

## 2. Shannon diversity per sample
shannon <- diversity(abund_t, index = "shannon")

## 3. Combine into a data.frame
alpha_div <- data.frame(
  Sample   = rownames(abund_t),
  Richness = richness,
  Shannon  = shannon
)

alpha_div


## Barplot of richness PER SAMPLE
barplot(
  alpha_div$Richness,
  names.arg = alpha_div$Sample,   # force sample IDs on x-axis
  las = 2,
  ylab = "vOTU richness (count)",
  main = "Per-sample viral richness"
)

## Barplot of Shannon PER SAMPLE
barplot(
  alpha_div$Shannon,
  names.arg = alpha_div$Sample,
  las = 2,
  ylab = "Shannon diversity",
  main = "Per-sample viral diversity"
)
```
### Defines sample groups


# STEP 22: FINAL ANALYSIS

Compare:
Viral contig counts (counts/viral_counts.tsv)
High-quality viral contigs (CheckV/)

### Goal:
Identify differences in phage abundance between MDD and controls

# Final Directory Structure

```bash
/home/NETID/FinalProject/
│
├── raw/                        # original FASTQ files from SRA
│   ├── ERR13348292_1.fastq.gz
│   ├── ERR13348292_2.fastq.gz
│   ├── ERR13348320_1.fastq.gz
│   ├── ERR13348320_2.fastq.gz
│   ├── ERR13348298_1.fastq.gz
│   ├── ERR13348298_2.fastq.gz
│   ├── ERR13348304_1.fastq.gz
│   └── ERR13348304_2.fastq.gz
│
├── trimmed/                    # Trimmomatic outputs
│   ├── ERR13348292_R1_trmPE.fq.gz
│   ├── ERR13348292_R1_trmSE.fq.gz
│   ├── ERR13348292_R2_trmPE.fq.gz
│   ├── ERR13348292_R2_trmSE.fq.gz
│   ├── ERR13348320_R1_trmPE.fq.gz
│   ├── ERR13348320_R1_trmSE.fq.gz
│   ├── ERR13348320_R2_trmPE.fq.gz
│   ├── ERR13348320_R2_trmSE.fq.gz
│   ├── ERR13348298_R1_trmPE.fq.gz
│   ├── ERR13348298_R1_trmSE.fq.gz
│   ├── ERR13348298_R2_trmPE.fq.gz
│   ├── ERR13348298_R2_trmSE.fq.gz
│   ├── ERR13348304_R1_trmPE.fq.gz
│   ├── ERR13348304_R1_trmSE.fq.gz
│   ├── ERR13348304_R2_trmPE.fq.gz
│   └── ERR13348304_R2_trmSE.fq.gz
│
├── fastqc_raw/                 # FastQC reports for raw reads
│   ├── ERR13348292_1_fastqc.html
│   ├── ERR13348292_2_fastqc.html
│   ├── ERR13348320_1_fastqc.html
│   ├── ERR13348320_2_fastqc.html
│   ├── ERR13348298_1_fastqc.html
│   ├── ERR13348298_2_fastqc.html
│   ├── ERR13348304_1_fastqc.html
│   └── ERR13348304_2_fastqc.html
│
├── fastqc_trimmed/             # FastQC reports for trimmed paired reads
│   ├── ERR13348292_R1_paired_fastqc.html
│   ├── ERR13348292_R2_paired_fastqc.html
│   ├── ERR13348320_R1_paired_fastqc.html
│   ├── ERR13348320_R2_paired_fastqc.html
│   ├── ERR13348298_R1_paired_fastqc.html
│   ├── ERR13348298_R2_paired_fastqc.html
│   ├── ERR13348304_R1_paired_fastqc.html
│   └── ERR13348304_R2_paired_fastqc.html
│
├── assembly/                   # MEGAHIT outputs
│   ├── slurm
│   ├── logs
│   ├── contigs
│   ├── ERR13348292/
│   │   └── final.contigs.fa
│   ├── ERR13348320/
│   ├── ERR13348298/
│   └── ERR13348304/
│
├── virsorter/                  # VirSorter2 outputs
│   ├── slurm
│   ├── logs
│   ├── db
│   ├── ERR13348292/
│   │   ├── final-viral-combined.fa
│   │   ├── final-viral-score.tsv
│   │   ├── final-viral-boundary.tsv
│   │   └── viral_10kb.fa
│   ├── ERR13348320/
│   ├── ERR13348298/
│   └── ERR13348304/
│
├── checkv/                     # CheckV outputs
│   ├──slurm
│   ├── db
│   ├── ERR13348292/
│   │   ├── quality_summary.tsv
│   │   ├── completeness.tsv
│   │   ├── contamination.tsv
│   │   └── viruses.fna
│   ├── ERR13348320/
│   ├── ERR13348298/
│   └── ERR13348304/
├── Bowtie-Samtools
│   ├── slurm
│   ├── logs
│   ├── votu.index*
│   └── contigs
└── CoverM
    ├── slurm
    ├── logs
    ├── bam_files
    └── coverm_output

```
