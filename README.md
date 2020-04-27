
# dAsthma source codes for "Latent-space embedding of expression data identifies gene signatures from sputum samples of asthmatic patients"




# Code for training and analysis of the dAsthma framework
## Model implementation
dA.cpp: implementation of the denoising autoencoder
Compile: g++ -o spuda dA.cpp
./spuda -h

## Downstream analysis
analysis: folder containing raw data, training results, code for analysis.
### Data and results
data: the data folder
These data will be loaded automatically in the relevant R scripts

Figures: figures created by the codes will be stored in this folder
genes_hidden_units: selected genes will be written into this folder for GO or network analysis, etc.

### Major codes
0_load_data.R: load and process data (processed data saved in data.RData)
1_gsea.R: run GSEA on the hidden units
2_predict_severity.R: prediction of severity (results saved in data.RData)
3_predict_clinical.R: prediction of FEV1/FVC ratio

### Spinoff codes for visualization 
plot_clinical.R: plot correlation between clinical traits
plot_GO.R: plot p-value for gene ontology results (downloaded from https://david.ncifcrf.gov/)
plot_network: plot gene interaction network (downloaded from https://string-db.org/)
