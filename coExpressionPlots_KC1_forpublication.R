###Katie Ching
###last updated 2020-08-07
###MAKING PLOTS OF GENE EXPRESSION FROM SCRNASEQ DATA,
###PLOTTING EXPRESSION OF 2 GENES AGAINST EACH OTHER,
###AND COLORING BY CELL TYPE

###The goal is to generate plots to see if cells are 
###expressing different genes at the same time or not.
###Dots for individual cells will be plotted and color
###coded to see if there is a difference in 
###co-expression between cell types.

###---###---###---###---###---###---###---###---###---###---

###SET-UP AND DEFINITIONS

##load the .Rda file in RStudio or use the following:
##example code for loading data below; if you have RStudio installed, then you can just drag the .Rda file into RStudio:
load("C:/Users/Katie/Desktop/oeHBCdiff_lineageData.Rda")

##Definitions:
#mlm: matrix of the microvillous lineage (genes X cells)
#nlm: matrix of the neuronal lineage
#slm: matrix of the sustentacular lineage
#mclus.labels, nclus.labels, sclus.labels: vectors of the cluster assigments for the cells in each respective lineage
#colpal, colpalM,N,&S are vectors of colors for plotting gene expression by cluster along pseudotime/developmental order for each lineage.

###---###---###---###---###---###---###---###---###---###---

#CREATE A FUNCTION THAT PLOTS EXPRESSION OF TWO DIFFERENT GENES
#WITHIN THE NEURONAL LINEAGE,
#ONE ON THE X AXIS, ONE ON THE Y AXIS

###Example of how to use this funtion:
###After running, type: coExpressionPlot("Plk4", "Stil")
###Plk4 will be on the horizontal axis, and Stil will be on the vertical axis.
coExpressionPlot <- function(gene1, gene2) {
  #Keep all RNA levels for gene1 and gene2.
  gene1RNAlevels <- as.vector(nlm[gene1,])
  gene2RNAlevels <- as.vector(nlm[gene2,])
  #This returns an error "subscript out of bounds" if the gene is not in the dataset
  
  #Create a matrix to hold gene information.
  if (length(gene1RNAlevels) == length(gene2RNAlevels)){
    coExpressionMatrix <- matrix(nrow = length(gene1RNAlevels), ncol = 3)
  } else {print("number of cells does not match for your two genes")}
  
  #Make a matrix with expression of gene1, 
  #gene2, and cell type.
  coExpressionMatrix[,1] <- gene1RNAlevels
  coExpressionMatrix[,2] <- gene2RNAlevels
  coExpressionMatrix[,3] <- as.vector(nclus.labels)
  
  #Make a scatterplot with gene1 on x, gene2 on y,
  #and color dots according to cell type.
  #Give the plot labels for x, y, and title.
  #set the axis lengths.
  plot(x = coExpressionMatrix[,1], y = coExpressionMatrix[,2],
       xlab = paste0(gene1, " log2 normalized RNA level"), 
       ylab = paste0(gene2, " log2 normalized RNA level"),
       col=colpalN[nclus.labels], pch=19,
       main = paste0("Co-Expression of ", gene1, " and ", gene2))
}

###---###---###---###---###---###---###---###---###---###---

