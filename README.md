This repository is for single cell RNA sequencing analysis related to Lawlor, Vanslambrouck, Higgins, et al, Nature Materials, 2020.
Imaging analysis is in the repository KidneyRegeneration/BioprintedOrganoids-imaging.

All data needed to replicate the single cell RNA seq analysis is at :

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152014

Each of the conditions, Ratio0 ('dots'), Ratio40 ('line'), Man ('hand') has an associated set of <a href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger">CellRanger</a> output matrices containing the mRNA data as well as a set of <a href="https://github.com/Hoohm/CITE-seq-Count">Cite-seq-Count</a> output matrices containing the hashing data.

The folder R-analysis contains the R notebooks that document the Seurat based analysis, including QC, clustering and differential expression testing.

R-analysis/Seurat-analysis-QC.Rmd describes the initial QC and incorporation of hash tag data.
R-analysis/Seurat-analysis.Rmd describes all other analysis.

Scrublet/run_scrublet.py is the script used to perform computational doublet prediction using <a href="https://github.com/AllonKleinLab/scrublet">Scrublet</a>.
