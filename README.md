# RNAseq

An analysis of bulk RNAseq data utilizing a prostate cancer dataset available at https://www.ncbi.nlm.nih.gov/sra?term=SRP060235.

The analysis covers: 
- differential expression
- geneset enrichment
- pathway analysis

## Install conda environment

```
conda env create -f environment.yml
```

## Run Jupyter Notebook 

Once you have your conda environment installed, one can interactively run through the Jupyter notebook and run the analysis. 

## Results 

Results can be found interactively throughout the notebook and for larger files, they are outputted into the `results/` directory. 

## Outcome 

With this analysis, we have came up with a list of genes that are differentially expressed within this prostate cancer dataset, done geneset enrichment analysis, and have idenfitifed the different KEGG pathways that our dataset has identified to be involved in the biology of this disease. 

Next steps, is to validate these potential biomarkers using the literature, screening for potential drug molecules, and continuing the drug discovery process. 
