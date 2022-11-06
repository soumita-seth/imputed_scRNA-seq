setwd("D:/PhD Work/Atikul work/ongoingwork")
#install.packages("usethis")
library(usethis)
#install.packages("devtools")
library(devtools)
#install_github("Vivianstats/scImpute")
library(scImpute)
#install.packages("survival")
library(survival)
#install.packages("dplyr")
library(plyr)
library(dplyr)
BiocManager::install(c("ggplot2"))
install.packages("ggbiplot")
library(ggplot2)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install_github("vqv/ggbiplot",force=TRUE)
BiocManager::install(c("ggbiplot"))
library(ggbiplot)
# install shrinkageClust pakckage and test it on marge dataset

path_to_file = 'D:/PhD Work/Atikul work/ongoingwork/Shrinkage-Clustering-master/shrinkageClust_0.1.0.tar.gz'
library(shrinkageClust)

#load marge dataset of Top 100 features

data1=read.csv("singlecellproj_atikul_sm_mrmrtop100genes (1).csv",header = TRUE,sep=",",check.names=TRUE)
dim(data1)#100 1873
#fix(data1)
View(data1)
data1$label
data1.cls<-c(rep("IOLEN2q0",97),rep("Lgr5",96),rep("Lgr5SC",192),rep("PSANNAq0",96),rep("Reg4_positive_cells_Replicate_1",96),rep("Reg4_positive_cells_Replicate_2",192),rep("Whole_Organoid_Replicate_1",288),rep("Whole_Organoid_Replicate_2",384),rep("YFPpos",432))
View(data1.cls)
#data1 transpose 
dd_t<-t(data1)
#data sheet view 

View(dd_t)
data2<-dd_t
Samples<-rownames(dd_t[2:1873,])
View(Samples)
class(data2)
data2.df<-as.data.frame(data2)
#data3<-data2.df
View(data2.df)
data2.df['label']<-data1.cls
dim(data2.df)
data2.df<-data2.df[2:1873,]
rownames(data2.df)<-NULL

label = data2.df$label
data2.df = data2.df[, c(1:100)]
#data2.df<-lapply(data2.df,as.numeric)
#data3 = data[, c(1:30)]
write.csv(data2.df,file = "finaldata2.csv", row.names = FALSE)
data3=read.csv("finaldata2.csv",header = TRUE)
View(data3)
# compute similarity matrix


S = simiMatrix(data3)
View(S)
dim(S)
# view distribution of similarity scores in S
hist(S)

# run shrinkage clustering
set.seed(10)
clust = SuperCluster(s=S,w=48,k=20,iter=1000,random=10)
 #clust = SuperCluster(s=S,w=100,k=20,iter=500,random=1)
clust_membership = clust$c[,10]
View(clust)
View(clust$c)
View(clust$a)
#clust_record<-clust$a
#clust_record['label']<-data2.df$label
#class(clust_record)
#View(clust_record)
rownames(clust$a)<-Samples
rownames(clust$c)<-Samples
write.csv(clust$c[,10],file = "Sample_Cluster_id2.csv", row.names = TRUE)
# use a confusion tablecompare the clustering solution with true membership
table(label,clust_membership)
View(clust$c)
View(clust$a)
# quantitatively assess the accuracy of the clustering result
eval_scores = evaluation(clust_membership, label)
#BiocManager::install(c("NMI"))
#library(NMI)
eval_scores$NMI # Normalized mutual information
eval_scores$RI # Rand index
eval_scores$F1 #F1 score
# generate a PCA plot to visualize the clustering result
scplot(data3, clust_membership)

library("cluster")

#Compute Silhoutte Index
sil <- silhouette(clust$c[,5],dist(data2.df))
si.sum <- summary(sil)
si.sum
si.sum$clus.avg.widths
si.sum$avg.width
si.sum$clus.sizes
plot(sil, main ="Silhouette plot - Shrinkage")
#BiocManager::install(c("factoextra"))
library(factoextra)
BiocManager::install(c("ggpubr"))
library(ggpubr)
BiocManager::install(c("car"))
library(car)
BiocManager::install(c("carData"))
library(carData)
fviz_silhouette(sil)
View(si.sum)
head(sil[, 1:3], 10)
head(sil,500)
sil
View(data1)
dim(data1)
data4<-data1[,2:1873]
rownames(data4)<-data1$label
View(data4)
data5<-read.csv("Sample_Cluster_id2.csv",header = TRUE)
View(data5)
class(data5)
#Choose best cluster as a diseased sample and rest as normal sample
diseased_samples<-data5[data5$ClusterId==6,]
View(diseased_samples)
dim(diseased_samples)
diseased_samples_Name<-diseased_samples$Sample
class(diseased_samples_Name)
All_samples<-colnames(data4)
View(All_samples)
Expr_samples<-All_samples[which(All_samples %in% diseased_samples_Name)]
View(Expr_samples)
Exper_data<- data4[,which(colnames(data4)%in%Expr_samples)]
View(Exper_data)
dim(Exper_data)
control_samples<-setdiff(All_samples,Expr_samples)
View(control_samples)
control_data<-data4[,which(colnames(data4)%in%control_samples)]
dim(control_data)
control.vs.diseased_data<-cbind(control_data,Exper_data)
x1<-control.vs.diseased_data
group <- as.factor(c(rep("control", 504),rep("diseased", 1368)))
colnames(x1)<-group
x1.ctr<- x1[,c(1:504)]
x1.exp<- x1[,c(505:1872)]
x1.ord<- cbind(x1.ctr,x1.exp)
control.vs.diseased1 <- x1.ord
View(control.vs.diseased1)
dim(control.vs.diseased1)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
control.vs.diseased1.Norm<-matrix(0, nrow = nrow(control.vs.diseased1), ncol = ncol(control.vs.diseased1))
dim(control.vs.diseased1.Norm)
for(i in 1:nrow(control.vs.diseased1)){
  print(paste("i:",i))
  #x.list[[i]] <- control.vs.diseased[i,]
  
  x_r <- control.vs.diseased1[i,]
  dim(x_r)###[1, 215]
  #View(x.list)
  #class(x.list)
  #min(x_r)
  #control.vs.diseased.Norm<- as.data.frame(lapply(x.list, normalize))
  #####control.vs.diseased.Norm[i,]<- as.data.frame(normalize(x_r))
  temp<-as.matrix(normalize(x_r))
  control.vs.diseased1.Norm[i,]<-temp 
}
rownames(control.vs.diseased1.Norm)<-rownames(control.vs.diseased1)
View(control.vs.diseased1.Norm)  
dim(control.vs.diseased1.Norm)
c1<-as.matrix(control.vs.diseased1.Norm)
View(c1)
library(edgeR)
library(limma)
RNA_Data<- c1
y_r <- DGEList(RNA_Data)
names(y_r)### See what slots are stored in y
y_r$samples
nrow(y_r$samples)
ctr_size<-ncol(x1.ctr)
sample_clslabel_ctrexp_r <- rep(c("Control","Diseased"),c(ctr_size,ncol(RNA_Data)-ctr_size))
sample_clslabel_ctrexp_r 
design_controlvsdiseased_r<-model.matrix(~ sample_clslabel_ctrexp_r)
#head(design_controlvsdiseased_r,n=15)
#par(mar = rep(2, 4))
#Apply Voom and Limma
v_ctrvsexp_r <- voom(y_r,design_controlvsdiseased_r,plot = TRUE)
v_ctrvsexp_r
names(v_ctrvsexp_r)
v_ctrvsexp_r$E
############Testing for differential expression for rna###
fit_ctrvsexp_r <- lmFit(v_ctrvsexp_r,design_controlvsdiseased_r)
names(fit_ctrvsexp_r)
fit_ebayes_ctrvsexp_r <- eBayes(fit_ctrvsexp_r)
dim(fit_ebayes_ctrvsexp_r)
options(digits=3)
top2_ctrvsexp_r <- topTable(fit_ebayes_ctrvsexp_r,number=Inf,adjust="BH",sort.by="P") ##use this if use "makeContrasts" function
#sum(top2_ctrvsexp_r$adj.P.Val<0.001) #??
View(top2_ctrvsexp_r)
DEGlist_ctrvsexp_r<-top2_ctrvsexp_r
DEGlist_srt_ctrvsexp_r<-DEGlist_ctrvsexp_r[order(DEGlist_ctrvsexp_r$adj.P.Val, decreasing = FALSE),]
nrow(DEGlist_srt_ctrvsexp_r)
View(DEGlist_srt_ctrvsexp_r)
#Got Differentially Expressed Markers 
write.table(DEGlist_srt_ctrvsexp_r,file="scRNAseq_Biomarkers_List(gene=100)_global.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)
#Comparitive study
Present_Markers <- read.csv("scRNAseq_Biomarkers_List(gene=100)_global.csv",header = TRUE)
Previous_Markers <- read.csv("previous_work_markers.csv",header = TRUE)
View(Present_Markers)
dim(Present_Markers)
View(Previous_Markers)
dim(Previous_Markers)
Common_Markers<-intersect(Present_Markers$Markers,Previous_Markers$gene)
View(Common_Markers)
length(Common_Markers)
write.csv(Common_Markers,file = "Common_Markers_new.csv", row.names = FALSE)
