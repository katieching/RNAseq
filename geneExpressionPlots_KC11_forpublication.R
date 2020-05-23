###Katie Ching
###started 2018-11-30, last updated 2020-05-21

###MAKING PLOTS OF GENE EXPRESSION FROM SCRNASEQ DATA,
###POOLING CELL TYPE DATA TO AVOID OVER-INTERPRETING WITHIN-TYPE TRENDS
###DERIVED FROM FLETCHER ET AL., 2017 
###https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(17)30127-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1934590917301273%3Fshowall%3Dtrue

###The goal is to generate plots for analyzing expression of specified genes
###across different cell types in the olfactory epithelium.
###Given that the exact position of cells within these cell types is not clear, 
###the following plots are intended to allow comparison between cell types 
###without trying to make comparisons within cell types.

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

###CREATE A FUNCTION THAT DRAWS A DOTPLOT FOR A SINGLE GENE
###WITH INDIVIDUAL CELLS AS DOTS IN COLUMNS
###BASED ON CELL TYPE

###Note: Cell types have now been re-ordered from Fletcher et al., 2017
###Note: This is all focused on the neuronal cell lineage.
###To do this for another cell lineage, change "nlm" and cell groups.

###Make a function with the mouse gene name as the input. 
###For example, type: geneDotPlot("Plk4")
geneDotPlot <- function(gene) {
  #Keep all RNA levels for gene.
  geneRNAlevels <- as.vector(nlm[gene,])
  #This returns an error "subscript out of bounds" if the gene is not in the dataset
  
  #Create a matrix as a reference, based on number code.
  ###NOTE: cell types do not match up with ascending order of cell type number
  ###(just corresponds to colors)
  ###Order seems to be: 1, 8, 5, 3, 2, 14, 10, 9, 12.
  CellTypeNumbers <- c(1, 8, 5, 3, 2, 14, 10, 9, 12)
  AllCellTypes <- c("HBC", "HBC1", "HBC2", "GBC", "INP1", "INP2", "INP3", "iOSN", "mOSN")
  #Create color palette vector in the correct order.
  orderedcolpalN <- c("#1B9E77", "#E6AB02", "#A6CEE3", "cyan", "antiquewhite2", 
                      "#A6761D", "darkorchid2", "#FFED6F", "#FF7F00")
  #Put all information into a matrix
  CellTypeDecoder <- matrix(nrow = length(AllCellTypes), ncol = 4)
  colnames(CellTypeDecoder) <- c("order", "number code ", "cell type name ", 
                                 "color code")
  CellTypeDecoder[,1] <- c(1:9)
  CellTypeDecoder[,2] <- CellTypeNumbers
  CellTypeDecoder[,3] <- AllCellTypes
  CellTypeDecoder[,4] <- orderedcolpalN
  
  
  #Create a matrix to hold gene information.
  ExpressionMatrix <- matrix(nrow = length(geneRNAlevels), ncol = 3)
  #Fill matrix with gene expression and cell type.
  ExpressionMatrix[,1] <- geneRNAlevels
  ExpressionMatrix[,2] <- as.vector(nclus.labels)
  
  #Add ordering value to the matrix based on cell type code.
  for (i in 1:length(ExpressionMatrix[,1])){
    
    #Find the cell type code in the 2nd column.
    celltype = ExpressionMatrix[,2][i]
    if (celltype == 1){
      ExpressionMatrix[,3][i] <- 1
    }
    if (celltype == 8){
      ExpressionMatrix[,3][i] <- 2
    }
    if (celltype == 5){
      ExpressionMatrix[,3][i] <- 3
    }
    if (celltype == 3){
      ExpressionMatrix[,3][i] <- 4
    }
    if (celltype == 2){
      ExpressionMatrix[,3][i] <- 5
    }
    if (celltype == 14){
      ExpressionMatrix[,3][i] <- 6
    }
    if (celltype == 10){
      ExpressionMatrix[,3][i] <- 7
    }
    if (celltype == 9){
      ExpressionMatrix[,3][i] <- 8
    }
    if (celltype == 12){
      ExpressionMatrix[,3][i] <- 9
    }
    #Put the order number into the third column
    
  }
  
  #Find average expression for each cell type.
  cell1 <- c()
  cell2 <- c()
  cell3 <- c()
  cell5 <- c()
  cell8 <- c()
  cell9 <- c()
  cell10 <- c()
  cell12 <- c()
  cell14 <- c()

  for (i in 1:length(nlm[gene,])){
    cellvector <- nlm[gene,][i]
    readcount = nlm[gene,i]
    cellID = names(cellvector)
    celltype = nclus.labels[cellID]
    if (celltype == 1){
      cell1 <- c(cell1, readcount)
    }
    if (celltype == 2){
      cell2 <- c(cell2, readcount)
    }
    if (celltype == 3){
      cell3 <- c(cell3, readcount)
    }
    if (celltype == 5){
      cell5 <- c(cell5, readcount)
    }
    if (celltype == 8){
      cell8 <- c(cell8, readcount)
    }
    if (celltype == 9){
      cell9 <- c(cell9, readcount)
    }
    if (celltype == 10){
      cell10 <- c(cell10, readcount)
    }
    if (celltype == 12){
      cell12 <- c(cell12, readcount)
    }
    if (celltype == 14){
      cell14 <- c(cell14, readcount)
    }
  }
  
  
  #find average expression of the gene in each cell type
  avg1 = sum(cell1)/length(cell1)
  avg2 = sum(cell2)/length(cell2)
  avg3 = sum(cell3)/length(cell3)
  avg5 = sum(cell5)/length(cell5)
  avg8 = sum(cell8)/length(cell8)
  avg9 = sum(cell9)/length(cell9)
  avg10 = sum(cell10)/length(cell10)
  avg12 = sum(cell12)/length(cell12)
  avg14 = sum(cell14)/length(cell14)
  AvgAllCells <- c(avg1, avg8, avg5, avg3, avg2, avg14, avg10, avg9, avg12)
  
  #Find standard deviation for each cell type.
  CellTypeSD <- c(sd(cell1), sd(cell8), sd(cell5), sd(cell3), sd(cell2),
                  sd(cell14), sd(cell10), sd(cell9), sd(cell12))
  
  #Put all names and values into a matrix.
  avg.mat <- matrix(nrow = 9, ncol = 3)
  ##Paste the names of cell types ino the first column.
  avg.mat[,1] <- paste(AllCellTypes)
  ##Paste the average expression for each cell type into the second column.
  avg.mat[,2] <- paste(AvgAllCells)
  ##Paste the standard deviation for each cell type in the third column.
  avg.mat[,3] <- paste(CellTypeSD)
  
  
  #Make a scatterplot with cell type on x axis, RNA levels on y
  #and color dots according to cell type.
  #Give the plot labels for x, y, and title.
  Graph <- plot(x = ExpressionMatrix[,3], y = ExpressionMatrix[,1],
                xlab = "", ylab = paste0(gene, " log2 normalized RNA level"),
                ylim = c(0, 12),
                col=colpalN[nclus.labels], pch=19,
                main = paste0("Expression of ", gene), xaxt = 'n')
  
  ###Add line segments for the error bars and "arrows" for the bar caps.
  ###Thanks to Claire Venard for help making nicer error bars!
  
  ###Arguments are in the form (x0, y0, x1, y1).
  #Draw a line at the average value.
  segments(as.numeric(CellTypeDecoder[,1]), AvgAllCells, as.numeric(CellTypeDecoder[,1])+0.2, 
           AvgAllCells, lwd = 3)
  #Draw "arrows" with flat arrow heads to indicate standard deviation.
  arrows(as.numeric(CellTypeDecoder[,1])+0.1, AvgAllCells - CellTypeSD, 
         as.numeric(CellTypeDecoder[,1])+0.1,
         AvgAllCells + CellTypeSD, lwd = 1.5, angle = 90,
         code = 3, length = 0.05)
  axis(side=1, CellTypeDecoder[,1], CellTypeDecoder[,3], col.axis = "black", par(las=2))
  #Note: par(las) refers to orientation of axis label. See "?par(las)" for details.
  
  
}
###---###---###---###---###---###---###---###---###---###---

#CREATE A FUNCTION THAT PLOTS EXPRESSION OF TWO DIFFERENT GENES
#WITHIN THE NEURONAL LINEAGE,
#ONE ON THE X AXIS, ONE ON THE Y AXIS

###make a function with inputs gene1 and gene2
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