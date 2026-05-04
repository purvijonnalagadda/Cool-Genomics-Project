# Liv, Purvi, Maeve & Sal: Depression vs. Control vOTU Analysis
#### Mora-Martinez et al.
## Details
- Our paper focused on the gut microbiomes of depression and obesity, and they ran a workflow to look at alpha and beta diversity of microbiomes between communities. Since this paper was published relatively recently, we wanted to look at it from a new angle. Specifically, we ran a workflow to look at vOTUs between depression and control groups on a subset of the data. 
- For the control group, we used sample alias' **102** & **128** and for the depression group we used sample alias **108** & **113** from this [https://www.ebi.ac.uk/ena/browser/view/PRJEB76994](url) page where their data is published
## Workflow Summary 
For this project, we analyzed viromics sequencing data between depression and control groups. Our workflow included downloading raw FASTQ files, assessing read quality with FastQC, trimming low quality bases and adaptors with trimmomatic, checking cleaned reads with FastQC, assembling reads into contigs with MEGAHIT, identifying viral contigs with VirSorter2, generating vOTUS with vClust, mapping reads with Bowtie-Samtools, generating .tsv files of vOTU abundance using CoverM, and creating graphs and a heatmap describing vOTU richness, Shannon diversity, and relative vOTU abundance in R. This repository contains notes, scripts and workflow documentation used to make the analysis reproducible. 

## Conclusion
This workflow allowed us to identify possible differences in viral DNA content of the microbiomes of depression vs control patients. In the future, we would like to continue this work on a larger subset of the data. Futhermore, we want to control for more factors (medication, diet, etc.) to see if certain vOTUs can be attributed to depression on a causation level, because at the moment our findings are purely correlational. 
