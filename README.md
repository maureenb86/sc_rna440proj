# sc_rna440proj
https://academic.oup.com/neuro-oncology/article/24/2/273/6291352?login=false

*Project Title: Building a medulloblastoma tumor subtype classifier through
analysis of immune cell expression patterns in the four known tumor subtypes*

## Overview:
Currently:
This repo contains preliminary analysis of the scRNA-seq data from the cited
Riemondy et al. paper.

Eventually:
This repo will contain scripts to traditionally analyze the scRNA-seq data
from the medulloblastoma study detailed below. Figures from the cited paper
will be recreated including tsne clustering of the immune cells into different
cell types. Novel analysis of the data will reveal the immune cell type ratio
of each of the four different subtypes of medulloblastoma. In a writeup that
will be included in the repo, we will report the statistical significance of
the expression of each immune cell type across the different tumors. Upon proof
of this hypothesis, we will develop a classifier for the various tumor types
which takes the immune cell expression data of a patient as input.

Citations for reused scripts will follow as needed throughout the project.

## Data:
Our project will be using data generated in Riemondy et al using
single-cell capture, RNA library preparation, and sequencing. ScRNA-seq was
performed on 2000 cells per sample with transcripts converted to cDNA, barcoded,
and libraries sequencing on Illumina NovaSeq6000, obtaining 50,000 reads/cell.
From this raw sequencing, tumor cells were identified via observation of
discrete clustering patterns and presence of copy number varients (CNVs) using
inferCNV. Tumor cells were then reanalyzed using Harmony alignment.

Data size is such that it is feasible to upload within code stored in repo.

**citation**: Kent A Riemondy, Sujatha Venkataraman, Nicholas Willard,
Anandani Nellan, Bridget Sanford, Andrea M Griesinger, Vladimir Amani,
Siddhartha Mitra, Todd C Hankinson, Michael H Handler, Martin Sill,
Jennifer Ocasio, Seth J Weir, Daniel S Malawsky, Timothy R Gershon,
Alexandra Garancher, Robert J Wechsler-Reya, Jay R Hesselberth,
Nicholas K Foreman, Andrew M Donson, Rajeev Vibhakar, Neoplastic and
immune single-cell transcriptomics define subgroup-specific intra-tumoral
heterogeneity of childhood medulloblastoma, Neuro-Oncology, Volume 24,
Issue 2, February 2022, Pages 273–286, https://doi.org/10.1093/neuonc/noab135

## Folder Structure:

This repo is divided into three main folders: figures, scripts, and data.

The scripts folder contains the .py files required to generate each figure
relevant to the report. Each script will make reference to the corresponding
figure generated.

The figures folder contains the .png files generated from the data throughout
the report.

The data folder is tracked with `git-lfs`. It contains two files: one with the
scRNA-seq data and one with the relevant metadata.

Other relevant files include the requirements.txt which will contain the
package and version requirements to run the project on your local machine.

## Installation:
Required packages:<br>
`Python 3.9.7`<br>
`numpy 1.21.2`<br>
`pandas 1.3.5`<br>
`scanpy 1.8.2`<br>
