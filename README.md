# MSc_DSB_Dissertation_2024
This is the official repository for the dissertation.

## Started maintaining the Git Repository (20-06-24)
- All work prior to this was done locally
- Most of the work involved working with the PCA plots and embedding

## Added preprocessing pipeline (27-06-24)
- Tested on lab scale on two samples
- Pipeline takes in fastq file and outputs counts matrix

## Added new code for different combination of PCA plots (28-06-24)
- In total five combination of PCA plots were added
- From the non subsetted no features defined PCA plot the outlying cells were identified and these cells were tracked in all the other plots.

## Completed PCA pipeline on GSE121708 (other mouse dataset) (29-06-24)
- Used the old established pipeline created earlier and adapted the code so that it works with this new data set.

## Seperated Seurat and Bioconductor Workflow in two different rmd files and generated 10 PCA plots with different combination of parameters and normalization (01-07-24)
- Started working on my specific part of the PCA plot, assigned by the supervisor

## Started the ranked SEG pipeline - Yet to finish it (01-07-24)
- Completed until SCE creation

## Removed the outlying cells and re-plot the PCA plot (02-07-24) 
- Completed my parts of the PCA plots for reproducing the old dissertation.

## Switched from bioinfmsc6 to bioinfmsc9 (03-07-24)
- Changed the primary server

## Created/Adjusted the Pre-processing pipeline to process the two human datasets (06-07-24)
- Downloaded the two human data sets on bioinfmsc9 (GSE36552, GSE71318)
- Ran them through the pre-processing pipeline (not yet complete)
- Improved the organization of files in the repository.

## Completed metadata creation script (07-07-24)
- Extracted metadata information from GEO
- Mapped Run Accession to GEO Accession
