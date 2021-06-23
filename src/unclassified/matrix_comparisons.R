library(MASS)
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
# ---- TEST ----
# --------27 10 2020---
DATA_RDS <- '/Volumes/kingdoms-workspace02/lehmann/projects/phd/symasym_sc_analysis/data/output/rds/SR50LR30/'
DATA_RDS <- '/kingdoms/csb/workspace2/lehmann/projects/phd/symasym_sc_analysis/data/output/rds/SR50LR30/'
sce <- readRDS(file = paste0(DATA_RDS, "sce_all.rds"))
mat <- as.matrix(GetAssayData(object = sce, slot = "counts"))
mathead <- head(mat)
mathead2 <- mathead[,1:5]

mathead2 <- as.data.frame(mathead2) %>% mutate(sumcells = rowSums(.))
mathead2 <- as.data.frame(mathead2$sumcells)
colnames(mathead2) <- 'm2'
mathead3 <- mathead2
colnames(mathead3) <- 'm3'
mathead1 <- mathead2
colnames(mathead1) <- 'm1'
mathead3['CD69L',] <- 6
mathead3['CARMIL3',] <- 7
mathead3['CLC2DL5',] <- 10



#matDij <- outer(mathead2, mathead3, Vectorize(euclidean_distance))

#x <- matrix(rnorm(100), nrow = 5)
#dist(t(x))

# --------

mathead <- cbind(mathead1, mathead2, mathead3)
mathead <- cbind(mathead1, mathead3)
mathead
dist(mathead, diag = T)


# --------
mcor <- cor(t(mathead))
mcor

saveRDS(mathead, file = 'tmp_mathead.rds')
# --------
mathead <- readRDS(file = 'scripts/R/tmp_mathead.rds')

ggplot(mathead, aes(x = m1, y = m3)) +
    geom_point() +
    geom_smooth(method = lm, size = 1)



# -------- with real data ----
mRef <- readRDS(file = '/kingdoms/csb/workspace2/lehmann/projects/phd/symasym_sc_analysis/data/output/rds/ref_NCBI/sce_all.rds')
mRef <- as.matrix(GetAssayData(object = mRef, slot = "counts"))
mRef <- as.data.frame(mRef) %>% mutate(Ref = rowSums(.))
mRef <- as.data.frame(mRef$Ref)
colnames(mRef) <- 'Ref'
mRef$genes <- rownames(mRef)

mSR50LR30 <- readRDS(file = '/kingdoms/csb/workspace2/lehmann/projects/phd/symasym_sc_analysis/data/output/rds/SR50LR30/sce_all.rds')
mSR50LR30 <- as.matrix(GetAssayData(object = mSR50LR30, slot = "counts"))
mSR50LR30 <- as.data.frame(mSR50LR30) %>% mutate(SR50LR30 = rowSums(.))
mSR50LR30 <- as.data.frame(mSR50LR30$SR50LR30)
colnames(mSR50LR30) <- 'mSR50LR30'
mSR50LR30$genes <- rownames(mSR50LR30)

mat <- cbind(mRef , mSR50LR30)
saveRDS(mat, file = '/kingdoms/csb/workspace2/lehmann/projects/phd/symasym_sc_analysis/data/output/rds/ref_NCBI/matrice_comp.rds')


# --------28 10 2020----
convertAnnotFull <- readRDS(paste0(DATA_RDS, "convertAnnotFull.rds"))
convertAnnot <- readRDS(paste0(DATA_RDS, "convertAnnot.rds"))

# [24,530] gene names
# gene_id [24,540]
convertAnnotFull <- convertAnnotFull %>% group_by(gene_id) %>% distinct()
test <- convertAnnotFull %>% group_by(gene_name) %>% distinct(gene_id)
nrow(test) - 24530
# result : 2294 duplicated XLOC


