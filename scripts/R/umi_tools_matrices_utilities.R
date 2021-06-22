#!/usr/bin/env Rscript
library(data.table)
library(scater)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)


# **************************************************************************** #
# Read GTF file as a dataframe
# **************************************************************************** #
readGTF <- function(GTFFile) {
    GTFTable <- read.table(GTFFile, sep="\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
    names(GTFTable) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    return(GTFTable)
}

# **************************************************************************** #
# Read UMI-tools sparse matrix as a dataframe
# **************************************************************************** #
readSparseMatrix <- function(filePath, fileName = 'counts.tsv.gz', compressed = TRUE, convertTable) {
    if (compressed) {
        sparseMatrix <- fread(paste0(filePath, fileName), header = TRUE, stringsAsFactors = TRUE)
    }
    else {
        sparseMatrix <- read.table(paste0(filePath, fileName), header = TRUE)
    }
    wideMatrix <- sparseMatrix %>% spread(cell, count, fill = 0)
    #matrixNames <- make.unique(as.character(convertTable[match(wideMatrix$gene, convertTable$transcript_id), 2]))
    #matrixNames <- make.unique(as.character(convertTable[match(wideMatrix$gene, convertTable$gene_id), 2]))
    matrixNames <- make.unique(as.character(convertTable[match(wideMatrix$gene, convertTable$gene_id), 1]))
    row.names(wideMatrix) <- matrixNames
    #print(wideMatrix[1:20, 1:4])
    #wideMatrix <- aggregate(wideMatrix[-which(colnames(wideMatrix) == 'gene')], list(gene = sub("\\.[^\\.]*$", "", row.names(wideMatrix))), sum)
    #print(wideMatrix[1:20, 1:4])
    #row.names(wideMatrix) <- wideMatrix$gene
    wideMatrix$gene <- NULL
    sce <- SingleCellExperiment(
        assays = list(counts = as.matrix(wideMatrix)),
        rowData = rownames(wideMatrix),
        colData = colnames(wideMatrix)
    )
    return(sce)
}


# **************************************************************************** #
# Build geneID convert table
# **************************************************************************** #
buildConvertTable <- function(GTFFile){
    convertTable <- GTFFile$attribute
    convertTable <- strsplit(convertTable, "\\s+|;")
    convertTable <- lapply(convertTable, function(x) x <- x[c(2,5,11,14,17)])
    convertTable <- data.frame(matrix(unlist(convertTable), nrow=length(convertTable), byrow = TRUE))
    colnames(convertTable) <- c('gene_id', 'transcript_id', 'gene_name', 'oId', 'nearest_ref')
    return(convertTable)
}

# **************************************************************************** #
### MAIN TEST ZONE
# **************************************************************************** #
#fileName = 'sc_merged.bam.featureCounts.counts.tsv.gz'
#filePath = '/import/kg_csbws02/lehmann/projects/phd/symasym_sc_analysis/data/output/matrices/ref_NCBI/'
#convertAnnot <- readRDS("/import/kg_csbws02/lehmann/projects/phd/symasym_sc_analysis/data/output/rds/ref_NCBI/convertAnnot.rds")
#sce <- readSparseMatrix(filePath,  fileName = 'sc_merged.bam.featureCounts.counts.tsv.gz', compressed = TRUE, convertAnnot)

