library(dplyr)
library(Seurat)
library (Matrix)
library(ggplot2)
library(sctransform)
library (EnhancedVolcano)
library (DoubletFinder)
library (pheatmap)
library(circlize)
library(plyr)
library(tidyverse) 
library(ggpubr)
library(cluster)
library(umap)
library(ggrepel)
library(Hmisc)
library(DoubletFinder)


#_____________________________Change this every time_____________________________
#write (brainID)_C(region)
Patient<-'288_C01'
#_____________________________Change this every time_____________________________


setwd("~/TREM2/Pre-processing")
#####Lines 31 - 80 should be run individually for each sample

####Load Dataset---- (Select path to the directory where your matrices are, it is the filtered_feature_bc_matrix_SAMPLE_NAME, and it 
####should contain 3 files - matrix, features and barcodes)
data <- Read10X(data.dir = paste('~/new_analysis_APOE/filtered_feature_bc_matrix_C',print(Patient), sep = ""))

####Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = data, project = paste('C', print(Patient), sep = ""), min.cells = 0, min.features = 0)

#### Store mt.percent
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")

###SUBSET low quality cells
plot1 <- FeatureScatter(seurat_obj , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#_____________________________Change this every time_____________________________
#write desired nFeature max cutoff given plot2
nFeature_max <- 9100
#_____________________________Change this every time_____________________________

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < nFeature_max  & percent.mt < 5)

## Pre-process Seurat object (sctransform)
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

####Doublet Finder
## pK Identification (no ground-truth)
sweep.res.list_seurat_obj <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = TRUE)
sweep.stats_seurat_obj <- summarizeSweep(sweep.res.list_seurat_obj, GT = FALSE)
bcmvn_seurat_obj <- find.pK(sweep.stats_seurat_obj)

## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(0.25)          
nExp_poi <- round(0.030*length(seurat_obj@active.ident))  ## Assuming 3% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#set pK value as value which results in max BCmetric
pK_BC_max<-droplevels(bcmvn_seurat_obj$pK[bcmvn_seurat_obj$BCmetric==max(bcmvn_seurat_obj$BCmetric)])
pK_value<-as.numeric(as.character(pK_BC_max))

## Run DoubletFinder with varying classification stringencies 
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

#look at metatdata in seurat_obj and change the string of numbers after pANN
df<-seurat_obj@meta.data %>% select(starts_with('pANN'))
pANN<-colnames(df)
new_pANN<-substring(pANN, 5)


####check your metadata for the pANN value and replace the next 3 lines of the script with the value for the sample
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = paste('pANN', new_pANN, sep = ""), sct = TRUE)
seurat_obj@meta.data[,"DF_hi.lo"] <- seurat_obj@meta.data[[paste('DF.classifications', new_pANN, sep = "")]]
seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet" & seurat_obj@meta.data[[paste('DF.classifications', new_pANN, sep = "")]] == "Singlet")] <- "Doublet_lo"
seurat_obj@meta.data$DF_hi.lo[which(seurat_obj@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"

assign(paste('C', print(Patient),sep = ""), seurat_obj)

Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj) 

assign(paste('S', print(Patient), sep = ""), seurat_obj)


#______________________________End of individual sample processing___________________



#### Merge objects and remove doublets----

seurat_combined = merge(S139_C12, y=c(S240_C12, S262_C12, S264_C12, S288_C12, S291_C12, S327_C12, S342_C12, S350_C12))

dim (seurat_combined)
Idents(seurat_combined) <- "DF_hi.lo"
seurat_combined <- subset(seurat_combined, idents = "Singlet")
dim (seurat_combined)
head (seurat_combined@meta.data)
tail (seurat_combined@meta.data)
.

####seurat_combined - add metadata ----
write.csv(seurat_combined@meta.data, file = "seurat_combined_metadata.csv")

####open csv file using excel and add the metadata (AAO, AAD, Mutation, Region, Diagnosis, etc), using the information at the inventory spreadsheet
#remove the pANN and DF.classifications metadata columns from the csv
seurat_combined.meta <- read.csv("~/TREM2/Pre-processing/seurat_combined_metadata.csv")
rownames(seurat_combined.meta)=colnames(seurat_combined)
seurat_combined <- AddMetaData(object=seurat_combined, metadata=seurat_combined.meta)
head(seurat_combined@meta.data)
tail(seurat_combined@meta.data)

seurat_combined@meta.data <- seurat_combined@meta.data[, -which(colnames(seurat_combined@meta.data) %in% 'X')]

head(seurat_combined@meta.data)
.
.
.
.
.
.
.

#Integrate to Perform Batch Correction----
#Integrate Data (Batch effect correction)
seurat_combined.list <- SplitObject(seurat_combined, split.by = "orig.ident")

#Normalize Data
seurat_combined.list <- lapply(X = seurat_combined.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat_combined.list, nfeatures = 3000)
seurat_combined.list <- PrepSCTIntegration(object.list = seurat_combined.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_combined.list, normalization.method = "SCT", 
                                        anchor.features = features)
seurat_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")


#scaledata and normalization step before PCA
seurat_combined<-ScaleData(seurat_combined)

#PCA Dimensionality Reduction
seurat_combined <- RunPCA(seurat_combined, features=VariableFeatures(seurat_combined))
ElbowPlot(seurat_combined)

#Visualization
seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:20)
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:20)
DimPlot(seurat_combined, reduction = "umap",label = TRUE) + NoLegend()
DimPlot(seurat_combined)


###test different resolutions, use high resulution (3-7) to get many clusters
seurat_combined <- FindClusters(seurat_combined, resolution=3, 
                                     verbose = FALSE,algorithm=1) 

Idents(seurat_combined) <- "seurat_clusters"
DimPlot (seurat_combined, label=T)

####Change default assay to RNA
DefaultAssay (seurat_combined) <- "RNA"

####Find all markers for each cluster
Cluster_markers_seurat_combined <- FindAllMarkers(seurat_combined, only.pos = TRUE, 
                                              logfc.threshold = 1.0)

write.csv(Cluster_markers_seurat_combined, file = "Cluster.markers_TREM2.csv")
####Explore the csv file for checking the gene markers for each cluster. 
#### 1) Check if there is any cluster that is enriched in MT- genes, or ribosomal genes (genes starting with RPS or RPL) 
#### if yes, you may want to remove these clusters.

cell_count <- table(seurat_combined@meta.data[["orig.ident"]])
cell_count<- as.data.frame(cell_count)
write_csv(cell_count, file = "cell_count_seurat_combined.csv")

#ID celltypes using celltype-ID.R script in /home/c_x_he/TREM2/R_Scripts