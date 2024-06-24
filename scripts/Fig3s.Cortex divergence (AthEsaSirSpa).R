#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(tidyverse)
library(readxl)
library(Augur)

source("scRNAseq_funs.R")




data.dir <- "path/to/data/dir"
analysis.dir <- "path/to/analysis/dir"
#####################################################################################################################################################################################

AthEsaSpa.CtrlCortex <- readRDS(file.path(data.dir, "AthEsaSpaCtrlCortex.RDS")) 
Sir.CtrlCortex <- readRDS(file.path(data.dir, "SirCtrlCortex.RDS"))
SirCtrlCortex.Group <-  readRDS("SirCtrlCortex_subpopulations.RDS")
Sir.IntCortex <- readRDS(file.path(analysis.dir, "intSirCtrlCortex.RDS")) %>% RunPCA(features = VariableFeatures(.), npcs = 50, umap.method = "umap-learn", metric = "correlation", verbose = FALSE) %>% RunUMAP(dims = 1:50, reduction = "pca", return.model = TRUE, verbose = FALSE)

if (TRUE) {
	### transfer sub-population annotation and umap metrics to the 1-to-1 ortholog based Sir control Cortex seuobj.
	Sir.CtrlCortex@meta.data <- left_join(Sir.CtrlCortex@meta.data %>% rownames_to_column(var = "cellid"), SirCtrlCortex.Group@meta.data %>% rownames_to_column(var = "cellid") %>% mutate(cellid = ifelse(str_sub(cellid, -1) == 1, paste0(cellid, "_17"), paste0(cellid, "_18"))), by = "cellid") %>% column_to_rownames(var = "cellid")
	Sir.CtrlCortex@reductions$umap@cell.embeddings <- SirCtrlCortex.Group@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column(var = "cellid") %>% mutate(cellid = ifelse(str_sub(cellid, -1) == 1, paste0(cellid, "_17"), paste0(cellid, "_18"))) %>% column_to_rownames(var = "cellid") %>% as.matrix()
	Sir.CtrlCortex@reductions$umap@misc$model <- Sir.IntCortex@reductions$umap@misc$model
	print(table(Sir.CtrlCortex@meta.data$group))

	### project Cortex cells from other species to Sir Cortex
	common_features <- lapply(list(Sir.CtrlCortex, AthEsaSpa.CtrlCortex), row.names) %>% Reduce(intersect, .) 
	int.anchors <- FindTransferAnchors(reference = Sir.CtrlCortex, query = AthEsaSpa.CtrlCortex, reference.assay = "SCT", query.assay = "SCT", dims = 1:50, normalization.method = "SCT", reference.reduction = "pca", features = common_features, recompute.residuals = FALSE)
	MappedCortex <- MapQuery(anchorset = int.anchors, reference = Sir.CtrlCortex, query = AthEsaSpa.CtrlCortex, refdata = list(group = "group"), reference.reduction = "pca", reduction.model = "umap")
}


#####################################################################################################################################################################################
### Identify genes that are differentially expressed between sub-populations
DefaultAssay(SirCtrlCortex) <- "RNA"
Idents(SirCtrlCortex) <- "group"

groupDE <- Reduce(function (x, y) {full_join(x, y, by = "geneID")}, 
    list(FindMarkers(SirCtrlCortex, slot = "data", ident.1 = "middle", ident.2 = "left") %>% dplyr::rename_with(~ paste0("middle_left.", .x)) %>% rownames_to_column(var = "geneID"), 
         FindMarkers(SirCtrlCortex, slot = "data", ident.1 = "middle", ident.2 = "right") %>% dplyr::rename_with(~ paste0("middle_right.", .x)) %>% rownames_to_column(var = "geneID"), 
         FindMarkers(SirCtrlCortex, slot = "data", ident.1 = "left", ident.2 = "right") %>% dplyr::rename_with(~ paste0("left_right.", .x)) %>% rownames_to_column(var = "geneID"))
    ) %>% column_to_rownames(var = "geneID")



#####################################################################################################################################################################################
q("no")

