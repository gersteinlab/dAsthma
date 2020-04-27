### Load and process data
# Reads raw data and autoencoder hidden layer values from the data folder, process them for downstream analysis
# Processed data stored in data.RData
Rscript 0_load_data.R

### Run GSEA on the hidden units
# Perform GSEA on the weights of the hidden units
# Load data from data.RData
Rscript 1_gsea.R

### Prediction of severity
# Load data from data.RData
# Perform random forest classificaiton on gene expression and hidden layer values to predict asthma severity
# Feature selection for most predictive genes
# Selected feature list stored in data_pred.RData
Rscript 2_predict_severity.R

### Prediction of FEV1/FVC ratio
# Load data from data_pred.RData
# Calculate correlation between hidden layer values and various clinical traits
# Perform svm regression on gene expression to predict FEV1/FVC ratio
Rscript 3_predict_clinical.R
