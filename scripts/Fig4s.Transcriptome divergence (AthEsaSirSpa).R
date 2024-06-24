#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(tidyverse)
library(readxl)
library(Augur)

source("scRNAseq_funs.R")


#####################################################################################################################################################################################
#####################################################################################################################################################################################
### Control samples only
#####################################################################################################################################################################################
AthEsaSirSpa_SCTIntegrated <- readRDS("integratedAthEsaSirSpa.RDS")
AthEsaSirSpa_SCTIntegrated_Ctrls <- subset(AthEsaSirSpa_SCTIntegrated, subset = treatment == "Ctrl")

input_AthEsaSirSpa_Ctrls <- GetAssayData(AthEsaSirSpa_SCTIntegrated_Ctrls, assay = "RNA", slot = "data")
meta_AthEsaSirSpa_Ctrls <- AthEsaSirSpa_SCTIntegrated_Ctrls@meta.data
print(sprintf("There are %s cells from control samples (4 species together).", nrow(meta_AthEsaSirSpa_Ctrls)))

#####################################################################################################################################################################################
#### between species for a given cell type
multiclassaugur_AthEsaSirSpa_Ctrls <- calculate_auc(input_AthEsaSirSpa_Ctrls, meta = meta_AthEsaSirSpa_Ctrls, label_col = "species", cell_type_col = "celltype", n_threads = 18)

pairwiseAugurList <- list()
comparisonlist <- list(c("Ath", "Esa"), c("Ath", "Sir"), c("Ath", "Spa"), c("Esa", "Sir"), c("Esa", "Spa"), c("Sir", "Spa"))
for (comparison in comparisonlist) {
	print(sprintf("Currently processing %s and %s ...", comparison[1], comparison[2]))
	seuobj_tmp <- subset(AthEsaSirSpa_SCTIntegrated_Ctrls, subset = species %in% comparison)
	input_tmp <- GetAssayData(seuobj_tmp, assay = "RNA", slot = "data")
	meta_tmp <- seuobj_tmp@meta.data
	print(sprintf("There are %s cells from %s and %s.", nrow(meta_tmp), comparison[1], comparison[2]))
	pairwiseaugur_tmp <- calculate_auc(input_tmp, meta = meta_tmp, label_col = "species", cell_type_col = "celltype", n_threads = 18)
	pairwiseAugurList[[paste(comparison, collapse = "_")]] <- pairwiseaugur_tmp
}

saveRDS(multiclassaugur_AthEsaSirSpa_Ctrls, file = "AthEsaSirSpa_Ctrls_multiclassAugur_species.RDS")
saveRDS(pairwiseAugurList, file = "AthEsaSirSpa_Ctrls_pairwiseAugur_species.RDS")

#####################################################################################################################################################################################
#### between cell types for a given species
multiclassaugur_AthEsaSirSpa_Ctrls <- calculate_auc(input_AthEsaSirSpa_Ctrls, meta = meta_AthEsaSirSpa_Ctrls, label_col = "celltype", cell_type_col = "species", n_threads = 18)

pairwiseAugurList <- list()
comparisonlist <- c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Xylem", "Phloem") %>% combn(., 2) %>% as.data.frame() %>% as.list()
for (comparison in comparisonlist) {
	print(sprintf("Currently processing %s and %s ...", comparison[1], comparison[2]))
	seuobj_tmp <- subset(AthEsaSirSpa_SCTIntegrated_Ctrls, subset = celltype %in% comparison)
	input_tmp <- GetAssayData(seuobj_tmp, assay = "RNA", slot = "data")
	meta_tmp <- seuobj_tmp@meta.data
	print(sprintf("There are %s cells from %s and %s.", nrow(meta_tmp), comparison[1], comparison[2]))
	pairwiseaugur_tmp <- calculate_auc(input_tmp, meta = meta_tmp, label_col = "celltype", cell_type_col = "species", n_threads = 18)
	pairwiseAugurList[[paste(comparison, collapse = "_")]] <- pairwiseaugur_tmp
}

saveRDS(multiclassaugur_AthEsaSirSpa_Ctrls, file = "220716_AthEsaSirSpa_Ctrls_multiclassAugur_celltype.RDS")
saveRDS(pairwiseAugurList, file = "220716_AthEsaSirSpa_Ctrls_pairwiseAugur_celltype.RDS")


#####################################################################################################################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################







### Control and treated samples together
#####################################################################################################################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################
#### between cell types for a given species
IntegratedAth <- readRDS("integratedAth.RDS")
IntegratedEsa <- readRDS("integratedEsa.RDS")
IntegratedSir <- readRDS("integratedSir.RDS")
IntegratedSpa <- readRDS("integratedSpa.RDS")
IntegratedCsa <- readRDS("integratedCsa.RDS")


multiclassAugurList <- list()
pairwiseAugurList <- list()
for (species in c("Ath", "Esa", "Sir", "Spa")) {
	print(sprintf("Currently processing %s ...", species))
	seuobj_tmp <- get(paste0("Integrated", species))
	input_tmp <- GetAssayData(seuobj_tmp, assay = "RNA", slot = "data")
	meta_tmp <- seuobj_tmp@meta.data
	print(sprintf("There are %s cells from %s.", nrow(meta_tmp), species))
	multiclassAugurList[[species]] <- calculate_auc(input_tmp, meta = meta_tmp, label_col = "intspace_treatment", cell_type_col = "intspace_celltype", n_threads = 26)
	for (treatmentinfo in c("5uM_ABA", "100uM_NaCl", "100mM_NaCl")) {
		print(sprintf("Currently processing %s in %s ...", treatmentinfo, species))
		treatedseuobj_tmp <- subset(seuobj_tmp, subset = intspace_treatment %in% c("Ctrl", treatmentinfo))
		treatedinput_tmp <- GetAssayData(treatedseuobj_tmp, assay = "RNA", slot = "data")
		treatedmeta_tmp <- treatedseuobj_tmp@meta.data
		print(sprintf("There are %s cells from Ctrl and %s in %s.", nrow(treatedmeta_tmp), treatmentinfo, species))
		pairwiseAugurList[[paste0(species, "_", treatmentinfo)]] <- calculate_auc(treatedinput_tmp, meta = treatedmeta_tmp, label_col = "intspace_treatment", cell_type_col = "intspace_celltype", n_threads = 26)
	}
}

saveRDS(multiclassAugurList, file = "AthEsaSirSpa_multiclassAugur_species.RDS")
saveRDS(pairwiseAugurList, file = "AthEsaSirSpa_pairwiseAugur_species.RDS")
#####################################################################################################################################################################################
#### between species for a given cell type under different conditions
AthEsaSirSpa_SCTIntegrated <- readRDS("integratedAthEsaSirSpa.RDS")


multiclassAugurList <- list()
pairwiseAugurList <- list()
for (celltypeinfo in c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Xylem", "Phloem")){
	print(sprintf("Currently processing %s ...", celltypeinfo))
	seuobj_tmp <- subset(AthEsaSirSpa_SCTIntegrated, subset = celltype == celltypeinfo)
	input_tmp <- GetAssayData(seuobj_tmp, assay = "RNA", slot = "data")
	meta_tmp <- seuobj_tmp@meta.data
	print(sprintf("There are %s cells from %s (4 species together).", nrow(meta_tmp), celltypeinfo))
	multiclassAugurList[[celltypeinfo]] <- calculate_auc(input_tmp, meta = meta_tmp, label_col = "treatment", cell_type_col = "species", n_threads = 26)
	for (treatmentinfo in c("5uM_ABA", "100uM_NaCl", "100mM_NaCl")) {
		print(sprintf("Currently processing %s for %s ...", treatmentinfo, celltypeinfo))
		treatedseuobj_tmp <- subset(seuobj_tmp, subset = treatment %in% c("Ctrl", treatmentinfo))
		treatedinput_tmp <- GetAssayData(treatedseuobj_tmp, assay = "RNA", slot = "data")
		treatedmeta_tmp <- treatedseuobj_tmp@meta.data
		print(sprintf("There are %s cells from Ctrl and %s in %s (4 species togther).", nrow(treatedmeta_tmp), treatmentinfo, celltypeinfo))
		pairwiseAugurList[[paste0(celltypeinfo, "_", treatmentinfo)]] <- calculate_auc(treatedinput_tmp, meta = treatedmeta_tmp, label_col = "treatment", cell_type_col = "species", n_threads = 26)
	}
}

saveRDS(multiclassAugurList, file = "AthEsaSirSpa_multiclassAugur_celltype.RDS")
saveRDS(pairwiseAugurList, file = "AthEsaSirSpa_pairwiseAugur_celltype.RDS")
#####################################################################################################################################################################################


q("no")

