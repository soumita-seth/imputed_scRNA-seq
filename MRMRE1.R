##setwd("C:/Users/akash/Desktop/project")
##setwd("C:/Users/saura/saurav backup lvo/R Script files")
setwd("D:/PhD Work/Atikul work/ongoingwork")
BiocManager::install(c("mRMRe"))
library(mRMRe)
library(survival)
library(devtools)
##install_github("Vivianstats/scImpute")
library(scImpute)
##install.packages("igraph")
library(igraph)

####increasing the memory limit#############
memory.limit()###16180
memory.limit(size=43000)###41000
##memory.limit(size=30000)


####single cell recorrected data file read##
#set.thread.count(2)
#data(cgps)
#Read Imputed Matrix
dd<-read.csv("scimpute_count.csv",header = TRUE,sep=",",check.names=TRUE,row.names = 1)
class(dd)##"data.frame"
dim(dd)##[23630, 1872]
View(dd)


########mrmr feature selection########
##https://stackoverflow.com/questions/48937143/selecting-features-from-a-feature-set-using-mrmre-package##
dd_t<-t(dd)
dim(dd_t)##[1872, 23630]
gr<-c(rep(1,96),rep(2,96),rep(3,192),rep(4,96),rep(5,96),rep(6,192),rep(7,288),rep(8,384),rep(-1,432))##change accordingly
length(gr)##1872
##str(dd)
dd_t_gr<-cbind(dd_t,gr)
dim(dd_t_gr)##[1872, 23631]

##dd_t_gr11<-dd_t_gr[,23531:23631]
##dim(dd_t_gr11)##[1872, 101]

##df[[7]] <- as.numeric(df[[7]])
myf_data <- mRMR.data(data = data.frame(dd_t_gr))
class(myf_data)##"mRMRe"
myresults <- mRMR.classic("mRMRe.Filter", data = myf_data, target_indices = ncol(dd_t_gr),
                        feature_count = 100)##change "feature_count" denoting how many top features are required
 solutions(myresults)
class(solutions(myresults))##"list"

####list to numeric conversion####
result_id<-as.numeric(as.character(unlist(solutions(myresults))))
class(result_id)##"numeric"
####

####final mrmr top features##
my_mrmr_topgene_dat<-t(dd_t_gr[,result_id])
dim(my_mrmr_topgene_dat)##[100,1872]/[500, 1872]
write.table(cbind(rownames(my_mrmr_topgene_dat),my_mrmr_topgene_dat),file="singlecellproj_atikul_sm_mrmrtop100genes.csv",sep=",",quote=F,row.names=F,col.names=T)
##write.table(cbind(rownames(my_mrmr_topgene_dat),my_mrmr_topgene_dat),file="singlecellproj_atikul_sm_mrmrtop500genes.csv",sep=",",quote=F,row.names=F,col.names=T)

####

