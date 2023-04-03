# imputed_scRNA-seq
#Single-cell data analysis and biomarker discovery using data imputation, feature selection and cell clustering  
Input Dataset: 
All datasets are available in online. We also added our all processed datasets in matrix format in the Data folder. 
Instructions for framework execution:
Run the file scimpute.R  with  raw dataset for computing imputed matrix. 
Then take Imputed data matrix as input data and run the file MRMRE1.R to fetch top 100 features with respect to their feature scores. 
After that run the file shrinkage_100genes_cluster_SS_Final_Version.R  with the input data matrix  100 × #cells . Cluster Validation code is also given in the same file. 
Instructions for Comparative Study:
Run the file shrinkage_100genes_cluster_Ven_SS.R to compare the output biomarkers of this framework with our previous one framework using Venn diagram. Both biomarker lists are provided in Data folder. 
Pre-requisites: 
R version 4.2.3
R packages: scImpute, mRMRe. shrinkageClust, cluster, edgeR, limma
