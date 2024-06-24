#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(tidyverse)
library(readxl)

source("scRNAseq_funs.R")
#####################################################################################################################################################################################
IntegratedAth <- readRDS("integratedAth.RDS")
IntegratedEsa <- readRDS("integratedEsa.RDS")
IntegratedSir <- readRDS("integratedSir.RDS")
IntegratedSpa <- readRDS("integratedSpa.RDS")
IntegratedCsa <- readRDS("integratedCsa.RDS")

for (species in c("Ath", "Esa", "Sir", "Spa", "Csa")) {
	print(sprintf("Currently processing %s ...", species))
	seuobj_tmp <- get(paste0("Integrated", species))
	
	if (species == "Csa") {
		pseudo.findDEG(seuobj_tmp, anncol = "celltype", concol = "treatment", repcol = "replicate", cond.1 = "Ctrl", method = "edgeR", pcut = 1, logFC.cutoff = 0) %>% saveRDS(., file = paste0(species, "_pseudobulkDE.RDS"))
	} else {
		pseudo.findDEG(seuobj_tmp, anncol = "intspace_celltype", concol = "intspace_treatment", repcol = "intspace_replicate", cond.1 = "Ctrl", method = "edgeR", pcut = 1, logFC.cutoff = 0) %>% saveRDS(., file = paste0(species, "_pseudobulkDE.RDS")) 
	}

	print(strrep("-", 100))
}

#####################################################################################################################################################################################

q("no")

