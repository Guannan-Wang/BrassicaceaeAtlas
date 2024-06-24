#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(tidyverse)
library(readxl)
library(Augur)

source("scRNAseq_funs.R")


IntegratedCsa <- readRDS("IntegratedCsa.RDS")
print(sprintf("There are %s cells from Csa.", nrow(IntegratedCsa@meta.data)))


#####################################################################################################################################################################################
#### for cell types within a given species in response to treatments
multiclassAugurList <- list()
pairwiseAugurList <- list()
for (species in c("Csa")) {
	print(sprintf("Currently processing %s ...", species))
	seuobj_tmp <- get(paste0("Integrated", species))
	input_tmp <- GetAssayData(seuobj_tmp, assay = "RNA", slot = "data")
	meta_tmp <- seuobj_tmp@meta.data
	multiclassAugurList[[species]] <- calculate_auc(input_tmp, meta = meta_tmp, label_col = "treatment", cell_type_col = "celltype", n_threads = 26)
	for (treatmentinfo in c("5uM_ABA", "100uM_NaCl", "100mM_NaCl")) {
		print(sprintf("Currently processing %s in %s ...", treatmentinfo, species))
		treatedseuobj_tmp <- subset(seuobj_tmp, subset = treatment %in% c("Ctrl", treatmentinfo))
		treatedinput_tmp <- GetAssayData(treatedseuobj_tmp, assay = "RNA", slot = "data")
		treatedmeta_tmp <- treatedseuobj_tmp@meta.data
		print(sprintf("There are %s cells from Ctrl and %s in %s.", nrow(treatedmeta_tmp), treatmentinfo, species))
		pairwiseAugurList[[paste0(species, "_", treatmentinfo)]] <- calculate_auc(treatedinput_tmp, meta = treatedmeta_tmp, label_col = "treatment", cell_type_col = "celltype", n_threads = 26)
	}
}

saveRDS(multiclassAugurList, file = "Csa_multiclassAugur_species.RDS")
saveRDS(pairwiseAugurList, file = "Csa_pairwiseAugur_species.RDS")




####################################################################################################################################################################################
#### divergence between cell types under different conditions
multiclassAugurList <- list()
pairwiseAugurList <- list()
comparisonlist <- c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Xylem", "Phloem") %>% combn(., 2) %>% as.data.frame() %>% as.list()
for (treatmentinfo in c("Ctrl", "5uM_ABA", "100mM_NaCl")) {  ##c("Ctrl", "5uM_ABA", "100uM_NaCl", "100mM_NaCl")
	print(sprintf("Currently processing %s ...", treatmentinfo))
	treatedseuobj_tmp <- subset(IntegratedCsa, subset = treatment == treatmentinfo)
	treatedinput_tmp <- GetAssayData(treatedseuobj_tmp, assay = "RNA", slot = "data")
	treatedmeta_tmp <- treatedseuobj_tmp@meta.data
	print(sprintf("There are %s cells from %s.", nrow(treatedmeta_tmp), treatmentinfo))
	tryCatch({
		multiclassAugurList[[treatmentinfo]] <- calculate_auc(treatedinput_tmp, meta = treatedmeta_tmp, label_col = "celltype", cell_type_col = "species", n_threads = 18)
		}, error = function(e){
			print(sprintf("ERROR during Augur multiclass testing for %s: %s", treatmentinfo, conditionMessage(e)))
			multiclassAugurList[[treatmentinfo]] <- NA
		}
	)
	
	for (comparison in comparisonlist) {
		print(sprintf("Currently processing %s and %s ...", comparison[1], comparison[2]))
		ctpairseuobj_tmp <- subset(treatedseuobj_tmp, subset = celltype %in% comparison)
		ctpairinput_tmp <- GetAssayData(ctpairseuobj_tmp, assay = "RNA", slot = "data")
		ctpairmeta_tmp <- ctpairseuobj_tmp@meta.data
		print(sprintf("There are %s cells from %s and %s.", nrow(ctpairmeta_tmp), comparison[1], comparison[2]))
		tryCatch({
			pairwiseaugur_tmp <- calculate_auc(ctpairinput_tmp, meta = ctpairmeta_tmp, label_col = "celltype", cell_type_col = "species", n_threads = 18)
			pairwiseAugurList[[treatmentinfo]][[paste(comparison, collapse = "_")]] <- pairwiseaugur_tmp
			}, error = function(e){
				print(sprintf("ERROR during Augur pairwise testing between %s and %s under %s: %s", comparison[1], comparison[2], treatmentinfo, conditionMessage(e)))
				pairwiseAugurList[[treatmentinfo]][[paste(comparison, collapse = "_")]] <- NA
			}
		)	
	}
	cat(rep("-", 100), sep = "", end = "\n")
}

saveRDS(multiclassAugurList, file = "Csa_multiclassAugur_celltypediv.RDS")
saveRDS(pairwiseAugurList, file = "Csa_pairwiseAugur_celltypediv.RDS")



#####################################################################################################################################################################################
#### for cell types within a given species in response to treatments using each sub-genome

### OLD subgenome information from: 10.1038/ncomms4706
Csa.subG <- list(Csa.subG1 = c("Csa17", "Csa16", "Csa15", "Csa04", "Csa13", "Csa11"), Csa.subG2 = c("Csa14", "Csa07", "Csa19", "Csa06", "Csa08", "Csa10", "Csa18"), Csa.subG3 = c("Csa03", "Csa05", "Csa01", "Csa09", "Csa20", "Csa02", "Csa12"))

### 10.1534/g3.119.400957
redefined.subG <- list(Csa.subG1 = c("Csa14", "Csa07", "Csa19", "Csa04", "Csa08", "Csa11"), Csa.subG2 = c("Csa03", "Csa16", "Csa01", "Csa06", "Csa13", "Csa10", "Csa18"), Csa.subG3 = c("Csa17", "Csa05", "Csa15", "Csa09", "Csa20", "Csa02", "Csa12"))

pairwiseAugurList <- list()
for (species in c("Csa")) {
	print(sprintf("Currently processing %s ...", species))
	seuobj_tmp <- get(paste0("Integrated", species))
	for (treatmentinfo in c("5uM_ABA", "100mM_NaCl")) {
		print(sprintf("Currently processing %s in %s ...", treatmentinfo, species))
		treatedseuobj_tmp <- subset(seuobj_tmp, subset = treatment %in% c("Ctrl", treatmentinfo))
		treatedinput_tmp <- GetAssayData(treatedseuobj_tmp, assay = "RNA", slot = "data")
		treatedmeta_tmp <- treatedseuobj_tmp@meta.data
		
		## run Augur with each sub-genome
		for (subG in c("subG1", "subG2", "subG3")) {
			geneIDs.subG <- grep(paste0(redefined.subG[[paste0("Csa.", subG)]], collapse = "|"), rownames(treatedinput_tmp), value = TRUE)
			subG.treatedinput_tmp <- treatedinput_tmp[geneIDs.subG, ]
		
			print(sprintf("For %s, there are %s features x %s cells (%s features x %s cells from all sub-genomes) from Ctrl and %s in %s.", subG, nrow(subG.treatedinput_tmp), ncol(subG.treatedinput_tmp), nrow(treatedinput_tmp), ncol(treatedinput_tmp), treatmentinfo, species))
			pairwiseAugurList[[paste0(species, subG, "_", treatmentinfo)]] <- calculate_auc(subG.treatedinput_tmp, meta = treatedmeta_tmp, label_col = "treatment", cell_type_col = "celltype", n_threads = 26)
		}
		
	}
}

saveRDS(pairwiseAugurList, file = "Csa.subG_pairwiseAugur_species.RDS")



#####################################################################################################################################################################################
q("no")

