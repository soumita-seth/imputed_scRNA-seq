BiocManager::install("VennDiagram")
BiocManager::install("futile.logger")
BiocManager::install("futile.options")
BiocManager::install("lambda.r")
library(lambda.r)
library(futile.options)
library(futile.logger)
library(VennDiagram)

# your data
#foo <- c('a','b','c','d')
#baa <- c('a','e','f','g')

Present_Markers <- read.csv("scRNAseq_Biomarkers_List(gene=100)_global.csv",header = TRUE)
Previous_Markers <- read.csv("previous_work_markers.csv",header = TRUE)
m1<-Present_Markers$Markers
m2<-Previous_Markers$gene

# Generate plot
v <- venn.diagram(list(Present_Markers=m1, Previous_Markers=m2),height=5000,width=1500, 
                  resolution = 1000,
                  fill = c("red", "blue"),
                  alpha = c(0.5, 0.5), cat.cex = 0.7, cex=0.5,
                  filename=NULL)

# have a look at the default plot
grid.newpage()
grid.draw(v)

# have a look at the names in the plot object v
lapply(v,  names)
# We are interested in the labels
lapply(v, function(i) i$label)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in foo only
v[[5]]$label  <- paste(setdiff(m1, m2), collapse="\n")  
# in baa only
v[[6]]$label <- paste(setdiff(m2, m1)  , collapse="\n")  
# intesection
v[[7]]$label <- paste(intersect(m1, m2), collapse="\n") 
# plot  
grid.newpage()
grid.draw(v)
