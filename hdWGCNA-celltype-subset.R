#https://github.com/smorabit/hdWGCNA/blob/dev/vignettes/basic_tutorial.Rmd

BrainTREM2 <- readRDS("~/TREM2/BrainTREM2.rds")

# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)

setwd("~/TREM2/Oli-subset")


seurat_obj<- subset(x = BrainTREM2, subset = All.Cell.type == "Oli")
DefaultAssay(seurat_obj)<- 'RNA'

seurat_obj<-NormalizeData(seurat_obj)
seurat_obj<-FindVariableFeatures(seurat_obj)
seurat_obj<-ScaleData(seurat_obj)
seurat_obj<-RunPCA(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj<-RunUMAP(seurat_obj, dims=1:20)
seurat_obj<-FindNeighbors(seurat_obj)
seurat_obj<-FindClusters(seurat_obj)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "BrainTREM2-Oli" # the name of the hdWGCNA experiment
)


seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Cases", "All.Cell.type"), # specify the columns in seurat_obj@meta.data to group by
  k = 20, # nearest-neighbors parameter, use 15 cells
  ident.group = 'All.Cell.type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

metacell_obj <- GetMetacellObject(seurat_obj)


seurat_obj <- seurat_obj %>%
  NormalizeMetacells() %>%
  ScaleMetacells(features=VariableFeatures(seurat_obj)) %>%
  RunPCAMetacells(features=VariableFeatures(seurat_obj)) %>%
  RunHarmonyMetacells(group.by.vars='Cases') %>%
  RunUMAPMetacells(reduction='harmony', dims=1:15)
p1 <- DimPlotMetacells(seurat_obj, group.by='All.Cell.type') + umap_theme() + ggtitle("Cell Type")
p2 <- DimPlotMetacells(seurat_obj, group.by='Cases') + umap_theme() + ggtitle("Cases")
p1 | p2

saveRDS(seurat_obj, file='hdWGCNA_after-RunUMAPMetacells.rds')
saveRDS(metacell_obj, file='hdWGCNA_metacell_obj.rds')


#seurat_obj <- readRDS("~/new_frontal_official/hdWGCNA/hdWGCNA_after-RunUMAPMetacells.rds")

#Co-expression of specific celltype

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Oli", # the name of the group of interest in the group.by column
  group.by='All.Cell.type' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)


## Select soft-power threshold
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # set this to FALSE since we did this above
)
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

#works now with lower k value


# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=7,
  setDatExpr=FALSE
)

dev.off()

PlotDendrogram(seurat_obj, main='Oli hdWGCNA Dendrogram')


# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="Diagnosis"
)


# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'All.Cell.type', group_name = 'Oli')

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Oli-M"
)

saveRDS(seurat_obj, file='hdWGCNA_Oli_post-eigengenes.rds')

# plot genes ranked by kME for each module
#p <- PlotKMEs(seurat_obj, ncol=5)
#p






#seurat_obj <- readRDS("~/New_Frontal_Pole/hdWGCNA/hdWGCNA_Oli_post-eigengenes.rds")

# get the module assignment table:
modules <- GetModules(seurat_obj)
# show the first 6 columns:
head(modules[,1:6])

write.csv(modules, "Oli_modules.subset-first.csv")

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# plot module correlagram
library(igraph)
library(corrplot)

ModuleCorrelogram(seurat_obj)


#add MEvalues to metadata

seurat_obj@meta.data <- cbind(
  seurat_obj@meta.data,
  GetMEs(seurat_obj, harmonized=TRUE)
)

# Plot Oli hME using Seurat VlnPlot function
p <- VlnPlot(seurat_obj,
             features = c('Oli-M1','Oli-M2','Oli-M3'),
             #'Oli-M1','Oli-M2','Oli-M3' 
             #'Oli-M4', 'Oli-M5', 'Oli-M6'
             #'Oli-M7'
             group.by = 'Diagnosis',
             pt.size = 0)
#p= p+geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()
# plot output
p


#wilcoxon rank sum test
#https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/#prerequisites
library(tidyverse)
library(rstatix)
library(ggpubr)


hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[["Oli-M7"]])
rownames(hME)=rownames(seurat_obj@meta.data)

# Show a sample of the data by group
set.seed(123)
hME %>% sample_n_by(diagnosis, size = 3)

#get summary data
hME %>%
  group_by(diagnosis) %>%
  get_summary_stats(hME, type = "median_iqr")

#signifiance test
stat.test <- hME %>% 
  rstatix::wilcox_test(hME ~ diagnosis) %>%
  add_significance()
stat.test

hME %>% wilcox_effsize(hME ~ diagnosis)

#make boxplot
library(ggplot2)

bxp <-ggboxplot(
  hME, x = "diagnosis", y = "hME", 
  ylab = "hME", xlab = "Diganosis", add = "jitter"
)
bxp

#add to graph
stat.test <- stat.test %>% add_xy_position(x = "diagnosis")
bxp + 
  stat_pvalue_manual(stat.test, tip.length = 0) #+
  #labs(subtitle = get_test_label(stat.test, detailed = TRUE))


saveRDS(seurat_obj, "hdWGCNA-seurat_obj-basic-end.rds")



seurat_obj <- readRDS("~/new_frontal_official/hdWGCNA-test/Oli-subset/hdWGCNA-seurat_obj-basic-end.rds")
