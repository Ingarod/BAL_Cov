
library('Seurat')
library('Matrix')
library('dplyr')
library('viridis')
library('reshape2')
library('data.table')
library('tidyverse')

setwd('~/Documents/CovIAV')

dat <- readRDS('Final_nCoV_0716_upload.RDS')

#Remove COV-5 severe patient

Idents(dat) <- 'batch'

dat@meta.data$group <- 'unknown'
 
dat@meta.data$group[grep(dat@meta.data$batch, pattern = 'Ctrl')] <- 'Healthy'
dat@meta.data$group[grep(dat@meta.data$batch, pattern = 'COV')] <- 'COVID-19'
dat@meta.data$group[grep(dat@meta.data$batch, pattern = 'Flu')] <- 'Influenza'
dat@meta.data$group[grep(dat@meta.data$batch, pattern = 'COV-5')] <- 'COVID-19-severe'

table(dat@meta.data$batch, dat@meta.data$group)
Idents(dat) <- 'group'
dat <- subset(dat, idents = c('Healthy', 'Influenza', 'COVID-19'))

#===============================================================================

#Subset on NK cells 

#===============================================================================

Idents(dat) <- 'cell_type'

dat_NK <- subset(dat, idents = c('NKs', 'XCL+ NKs'))

dat_NK <- SCTransform(dat_NK, vars.to.regress = 'percent.mt', verbose = FALSE)
DefaultAssay(dat_NK) <- 'SCT'
dat_NK <- RunPCA(dat_NK, verbose = FALSE)
ElbowPlot(dat_NK, ndims = 25)
dat_NK <- RunUMAP(dat_NK, dims = 1:20, verbose = FALSE) %>%
	FindNeighbors(dims = 1:20, verbose = FALSE) %>%
	FindClusters(resolution = 0.8, verbose = FALSE)

str(dat_NK@meta.data)

DimPlot(dat_NK, pt.size = 1, group.by = 'SCT_snn_res.0.8')

DimPlot(dat_NK, pt.size = 1, group.by = 'group')


#RNA assay 
DefaultAssay(dat_NK) <- 'RNA'
dat_NK <- NormalizeData(dat_NK)
all.genes <- rownames(dat_NK)
dat_NK <- ScaleData(dat_NK, features = all.genes, vars.to.regress = 'percent.mt')

FeaturePlot(dat_NK, features = 'FCGR3A', pt.size = 1) 

NK_markers <- FindAllMarkers(dat_NK, min.pct = 0.25, logfc.threshold = 0.25)
NK_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(dat_NK, features = top10$gene, group.by = 'SCT_snn_res.0.8') + theme(text = element_text(size = 8)) 


#CD16 negative
CD16neg <- subset(dat_NK, idents = c('7', '4', '6'))
table(CD16neg@meta.data$group)
# COVID-19   Healthy Influenza
#      829       114       271


DefaultAssay(CD16neg) <- 'RNA'
CD16neg <- NormalizeData(CD16neg)
all.genes <- rownames(CD16neg)
CD16neg <- ScaleData(CD16neg, features = all.genes, vars.to.regress = 'percent.mt')

chem <- c('CXCR3', 'CXCR6', 'CCR2', 'CCR5')

Idents(CD16neg) <- 'group'
levels(CD16neg) <- c('COVID-19', 'Influenza', 'Healthy')
DotPlot(CD16neg, features = chem, cols = 'RdYlBu', dot.scale = 35)+ 
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1))+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD16neg_DotPlot_chem_CovIAV.pdf', w = 12, h = 8)



#Save average expression per donor 
Idents(CD16neg) <- 'batch'

table(CD16neg@meta.data$batch)
CD16neg_avg <- AverageExpression(CD16neg)

CD16neg_chem <- CD16neg_avg$RNA[chem, ]
head(CD16neg_chem)

#Add group 
CD16neg_chem <- data.frame(t(CD16neg_chem))
CD16neg_chem$group <- 'unknown'
CD16neg_chem[grep('Ctrl', rownames(CD16neg_chem)), ]$group <- 'Healthy'
CD16neg_chem[grep('Flu', rownames(CD16neg_chem)), ]$group <- 'Influenza'
CD16neg_chem[grep('COV', rownames(CD16neg_chem)), ]$group <- 'COVID-19'

write.csv(CD16neg_chem, 'CD16neg_avg_chem_CovIAV.csv')


#CD16positive 
CD16pos <- subset(dat_NK, idents = c('7', '4', '6'), invert = TRUE)
table(CD16pos@meta.data$group)
# COVID-19   Healthy Influenza
#     2288      1481       612


DefaultAssay(CD16pos) <- 'RNA'
CD16pos <- NormalizeData(CD16pos)
all.genes <- rownames(CD16pos)
CD16pos <- ScaleData(CD16pos, features = all.genes, vars.to.regress = 'percent.mt')

Idents(CD16pos) <- 'group'
levels(CD16pos) <- c('COVID-19', 'Influenza', 'Healthy')

DotPlot(CD16pos, features = chem, cols = 'RdYlBu', dot.scale = 35)+ 
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1))+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD16pos_DotPlot_chem_CovIAV.pdf', w = 12, h = 8)

#Save average expression per donor 
Idents(CD16pos) <- 'batch'

table(CD16pos@meta.data$batch)
CD16pos_avg <- AverageExpression(CD16pos)

CD16pos_chem <- CD16pos_avg$RNA[chem, ]
head(CD16pos_chem)

#Add group 
CD16pos_chem <- data.frame(t(CD16pos_chem))
CD16pos_chem$group <- 'unknown'
CD16pos_chem[grep('Ctrl', rownames(CD16pos_chem)), ]$group <- 'Healthy'
CD16pos_chem[grep('Flu', rownames(CD16pos_chem)), ]$group <- 'Influenza'
CD16pos_chem[grep('COV', rownames(CD16pos_chem)), ]$group <- 'COVID-19'

write.csv(CD16pos_chem, 'CD16pos_avg_chem_CovIAV.csv')

#===============================================================================

#Subset on T cells 

#===============================================================================

table(dat@meta.data$cell_type)

str(dat@meta.data)

Idents(dat) <- 'cell_type'

dat_T <- subset(dat, idents = c('Cytotoxic CD8 T cells', 'Naive T cells', 'Activated CD4 T cells', 'Cycling T cells'))
str(dat_T@meta.data)

#SCTransform for clustering

dat_T <- SCTransform(dat_T, vars.to.regress = 'percent.mt', verbose = FALSE)
DefaultAssay(dat_T) <- 'SCT'
dat_T <- RunPCA(dat_T, verbose = FALSE)
ElbowPlot(dat_T, ndims = 25)
dat_T <- RunUMAP(dat_T, dims = 1:20, verbose = FALSE) %>%
	FindNeighbors(dims = 1:20, verbose = FALSE) %>%
	FindClusters(verbose = FALSE)

DimPlot(dat_T, pt.size = 0.8)

DimPlot(dat_T, pt.size = 0.8, group.by = 'cell_type')

DimPlot(dat_T, pt.size = 0.8, group.by = 'group')


#Normalise and scale RNA assay for visualisation and DEGs
DefaultAssay(dat_T) <- 'RNA'
dat_T <- NormalizeData(dat_T)
all.genes <- rownames(dat_T)
dat_T <- ScaleData(dat_T, features = all.genes, vars.to.regress = 'percent.mt')

FeaturePlot(dat_T, features = c('CD8A', 'LYZ','CD4', 'NKG7', 'IL7R', 'CCR7', 'KLRB1'))


#CD8 T cells
Idents(dat_T) <- 'cell_type'
CD8T <- subset(dat_T, idents = c('Cytotoxic CD8 T cells', 'Cycling T cells'))

DefaultAssay(CD8T) <- 'RNA'
CD8T <- NormalizeData(CD8T)
all.genes <- rownames(CD8T)
CD8T <- ScaleData(CD8T, features = all.genes, vars.to.regress = 'percent.mt')


str(CD8T@meta.data)
Idents(CD8T) <- 'group'
DotPlot(CD8T, features = chem, cols = 'RdYlBu', dot.scale = 35)+ 
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1))+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD8T_DotPlot_chem_CovIAV.pdf', w = 12, h = 8)

#Save average expression per donor 
Idents(CD8T) <- 'batch'

table(CD8T@meta.data$batch)
CD8T_avg <- AverageExpression(CD8T)

CD8T_chem <- CD8T_avg$RNA[chem, ]

#Add group 
CD8T_chem <- data.frame(t(CD8T_chem))
CD8T_chem$group <- 'unknown'
CD8T_chem[grep('Ctrl', rownames(CD8T_chem)), ]$group <- 'Healthy'
CD8T_chem[grep('Flu', rownames(CD8T_chem)), ]$group <- 'Influenza'
CD8T_chem[grep('COV', rownames(CD8T_chem)), ]$group <- 'COVID-19'

write.csv(CD8T_chem, 'Jan_22/CD8T_avg_chem_CovIAV.csv')



#CD4T cells 

CD4T <- subset(dat_T, idents = c('Naive T cells', 'Activated CD4 T cells'))

DefaultAssay(CD4T) <- 'RNA'
CD4T <- NormalizeData(CD4T)
all.genes <- rownames(CD4T)
CD4T <- ScaleData(CD4T, features = all.genes, vars.to.regress = 'percent.mt')
 

Idents(CD4T) <- 'group'
DotPlot(CD4T, features = chem, cols = 'RdYlBu', dot.scale = 35)+ 
	theme(axis.text.x = element_text(size = 20, angle = 45, hjust=1))+ 
	theme(text = element_text(size = 30))+
	theme(axis.text.y = element_text(size = 20))+
	theme(legend.text = element_text(size = 20))+
	theme(legend.key.size = unit(1, 'cm'))
ggsave('CD4T_DotPlot_chem_CovIAV.pdf', w = 12, h = 8)

#Save average expression per donor 
Idents(CD4T) <- 'batch'

CD4T_avg <- AverageExpression(CD4T)

CD4T_chem <- CD4T_avg$RNA[chem, ]

#Add group 
CD4T_chem <- data.frame(t(CD4T_chem))
CD4T_chem$group <- 'unknown'
CD4T_chem[grep('Ctrl', rownames(CD4T_chem)), ]$group <- 'Healthy'
CD4T_chem[grep('Flu', rownames(CD4T_chem)), ]$group <- 'Influenza'
CD4T_chem[grep('COV', rownames(CD4T_chem)), ]$group <- 'COVID-19'

write.csv(CD4T_chem, 'Jan_22/CD4T_avg_chem_CovIAV.csv')










