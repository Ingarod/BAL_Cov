
library('Seurat')
library('Matrix')
library('dplyr')
library('viridis')
library('reshape2')
library('data.table')
library('tidyverse')

setwd('~/BAL_cov/')

dat <- readRDS('nCoV.rds')

#Annotate
anno <- read.delim('all.cell.annotation.meta.txt')

anno$celltype <- as.character(anno$celltype)

idx <- which(dat@meta.data$ID %in% anno$ID)

dat@meta.data$celltype <- 'unknown'

dat@meta.data$celltype[idx] <- anno$celltype 

DimPlot(dat, group.by = 'celltype')


DefaultAssay(dat) <- 'RNA'

cyt <- c('CXCL9', 'CXCL10', 'CXCL11', 'CXCL16', 'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CCL7', 'CCL8', 'CCL13')
Idents(dat) <- 'group'

levels(dat) <- c('S/C', 'O', 'HC') 

DotPlot(dat, features = cyt, cols = 'RdYlBu', dot.scale = 25)+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('Total_cytokines_BalCov.pdf', w = 15, h = 10)



#===================================================================================================

#T and NK

#===================================================================================================

Idents(dat) <- 'celltype'
dat_NKT <- subset(dat, idents = c('T', 'NK'))
dim(dat_NKT)

#Annotate based on author's NK and T cell annotation
anno <- read.delim('NKT.cell.annotation.meta.txt')

anno$celltype <- as.character(anno$celltype)

idx <- which(dat_NKT@meta.data$ID %in% anno$ID)

dat_NKT@meta.data$celltype2 <- 'unknown'

dat_NKT@meta.data$celltype2[idx] <- anno$celltype


dat_NKT <- SCTransform(dat_NKT, vars.to.regress = 'percent.mito', verbose = FALSE)
DefaultAssay(dat_NKT) <- 'SCT'
dat_NKT <- RunPCA(dat_NKT, verbose = FALSE)
ElbowPlot(dat_NKT, ndims = 25)
dat_NKT <- RunUMAP(dat_NKT, dims = 1:20, verbose = FALSE) %>%
	FindNeighbors(dims = 1:20, verbose = FALSE) %>%
	FindClusters(verbose = FALSE)

#Verify cell identification
DimPlot(dat_NKT, group.by = 'celltype2', pt.size = 0.8)
DimPlot(dat_NKT, group.by = 'group')

FeaturePlot(dat_NKT, features = c('CD8A', 'CD4', 'IL7R',  'CCR7', 'FCGR3A', 'NCR1', 'GNLY'))

DimPlot(dat_NKT, group.by = 'SCT_snn_res.0.8') 
ggsave('NKT_UMAP_0.8_BalCov.pdf', h = 7, w = 7)



#===================================================================================================

#NK cells

#===================================================================================================

Idents(dat_NKT) <- 'celltype2'
dat_NK <- subset(dat_NKT, idents = 'NK')

dat_NK <- SCTransform(dat_NK, vars.to.regress = 'percent.mito', verbose = FALSE)
DefaultAssay(dat_NK) <- 'SCT'
dat_NK <- RunPCA(dat_NK, verbose = FALSE)
ElbowPlot(dat_NK, ndims = 25)
dat_NK <- RunUMAP(dat_NK, dims = 1:20, verbose = FALSE) %>%
	FindNeighbors(dims = 1:20, verbose = FALSE) %>%
	FindClusters(verbose = FALSE, resolution = 1.5)
DimPlot(dat_NK, label = TRUE)
DimPlot(dat_NK, group.by = 'group')
DimPlot(dat_NK, group.by = 'sample_new')


#Normalise and scale for gene expression visualisation
DefaultAssay(dat_NK) <- 'RNA'
dat_NK <- NormalizeData(dat_NK)
all.genes <- rownames(dat_NK)
dat_NK <- ScaleData(dat_NK, features = all.genes, vars.to.regress = 'percent.mito')

FeaturePlot(dat_NK, features = c('FCGR3A', 'CXCR3'))


Idents(dat_NK) <- 'group'
levels(dat_NK) <- c('S/C', 'O', 'HC') 
chem <- c('CXCR3', 'CXCR6', 'CCR2', 'CCR5')

DotPlot(dat_NK, features = c('FCGR3A', chem), cols = 'RdYlBu', dot.scale = 35)+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1))+	
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('NK_chemr_BalCov.pdf', w = 12, h = 10)

ef <- c('GZMA', 'GZMB', 'PRF1')

Idents(dat_NK) <- 'group'

levels(dat_NK) <- c('HC', 'O', 'S/C') 

VlnPlot(dat_NK, features = ef, cols = c("#999999", '#F49419', '#E85304'))
ggsave('NK_ef_Vln_BalCov.pdf', w = 16, h = 8)


#Additional chemokine receptors
chem_2 <- c('CCR1', 'CCR2','CCR3', 'CCR4', 'CCR5','CCR6', 'CCR7', 'CCR8', 'CCR9', 'CCR10', 'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CXCR5', 'CXCR6', 'GPR35', 'XCR1', 'CX3CR1', 'ACKR2', 'ACKR3', 'ACKR4')
levels(dat_NK) <- c('S/C', 'O', 'HC') 
DotPlot(dat_NK, features = chem_2, cols = 'RdYlBu',  dot.scale = 15) + 
	theme(text = element_text(size = 30))+
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('NK_AllChemR_BalCov.pdf',  w = 15, h = 10)
 


#===================================================================================================

#T cells

#===================================================================================================

Idents(dat_NKT) <- 'celltype2'
table(dat_NKT@meta.data$celltype2)
CD8T <- subset(dat_NKT, idents = 'CD8 T') 
DefaultAssay(CD8T) <- 'RNA'
CD8T <- NormalizeData(CD8T)
all.genes <- rownames(CD8T)
CD8T <- ScaleData(CD8T, features = all.genes, vars.to.regress = 'percent.mito')

Idents(CD8T) <- 'group'
levels(CD8T) <- c('S/C', 'O', 'HC') 

DotPlot(CD8T, features = chem, cols = 'RdYlBu', dot.scale = 35)+ 
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1))+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD8T_DotPlot_chem_BalCov.pdf', w = 12, h = 10)


levels(CD8T) <- c('HC', 'O', 'S/C') 

VlnPlot(CD8T, features = ef, cols = c("#999999", '#F49419', '#E85304'))
ggsave('CD8T_ef_Vln_BalCov.pdf', w = 16, h = 8)

#Additional chemokine receptors
DotPlot(CD8T, features = chem_2, cols = 'RdYlBu',  dot.scale = 15) + 
	theme(text = element_text(size = 30))+
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD8T_AllChemR_BalCov.pdf',  w = 15, h = 10)


CD4T <- subset(dat_NKT, idents = c('CCR7+ T', 'Treg'))

DefaultAssay(CD4T) <- 'RNA'
CD4T <- NormalizeData(CD4T)
all.genes <- rownames(CD4T)
CD4T <- ScaleData(CD4T, features = all.genes, vars.to.regress = 'percent.mito')

Idents(CD4T) <- 'group'
levels(CD4T) <- c('S/C', 'O', 'HC') 
DotPlot(CD4T, features = chem, cols = 'RdYlBu', dot.scale = 35)+ 
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1))+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD4T_DotPlot_chem_BalCov.pdf', w = 12, h = 10)


DotPlot(CD4T, features = chem_2, cols = 'RdYlBu',  dot.scale = 15) + 
	theme(text = element_text(size = 30))+
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD4T_AllChemR_BalCov.pdf',  w = 15, h = 10)
