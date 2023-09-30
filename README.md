# imputed_scRNA-seq
#Single-cell data analysis and biomarker discovery using data imputation, feature selection and cell clustering  

Input Dataset: 
All datasets are available in online. We also added our all processed datasets in matrix format in the Data folder. 

Instructions for framework execution:
Open the Code folder.

Run the file scimpute.R  with  raw dataset for computing imputed matrix by scImpute Algorithm. 
Then take Imputed data matrix as input data and run the file MRMRE1.R to fetch top 100 features with respect to their feature scores. 
After that run the file shrinkage_cluster_SS_Final_Version.R  with the input data matrix  100 × #cells . As well as run the same file with the input data matrix 500 ×  #cells. Choose the best one. 
Cluster Validation code is also given in the same file. 

Instructions for Comparative Study:
Run the file shrinkage_100genes_cluster_Ven_SS.R to compare the output biomarkers of this framework with our previous one framework using Venn diagram. Both biomarker lists are provided in Data folder. 
Run the file UpRegulated Markers_Ven.R to compare the output upregulated markers from Imputed data and Raw data. 
Run the file Magic_Impute.py to impute dataset using MAGIC imputation algorithm. Then run the file MRMRE1.R 
 and shrinkage_cluster_SS_Final_Version.R with magic imputed dataset. 

Instructions for Ablation Study:
Case 1: Run the file MRMRE1.R  and shrinkage_cluster_SS_Final_Version.R with raw dataset
Case 2: Run the file shrinkage_cluster_SS_Final_Version.R with raw dataset and all features (without imputation and feature selection)

Pre-requisites: 
R version 4.2.3
R packages: scImpute, mRMRe. shrinkageClust, cluster, edgeR, limma
Python package: magic
