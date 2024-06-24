#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(symphony)
library(irlba)
library(tidyverse)
library(readxl)
library(sctransform)
library(glmGamPoi)
library(future)


#######################################################################################################################################
## 1. annotation transfer using Seurat
#######################################################################################################################################
### annotation transfer with different parameters
Ath_Ctrls_cCellTypes <- readRDS("Ath_Ctrls_cCellTypes.RDS") ### Ath_Ctrls_cCellTypes.RDS is generated in consensusann4Ath.R


AthEsaSirSpa_SCTIntegrated <- readRDS("integratedAthEsaSirSpa.RDS")
AthEsaSirSpa_SCTIntegrated@meta.data <- AthEsaSirSpa_SCTIntegrated@meta.data %>% mutate(spesample = str_sub(orig.ident, 1, -3), sample = gsub("AthEsaSirSpa_", "", orig.ident)) %>% mutate(species = str_sub(sample, 1, 3), treatment = str_sub(sample, 4, nchar(sample)-2), replicate = str_sub(sample, -2, -1)) %>% mutate(treatment = lapply(treatment, function(x) {ifelse(x == "C", "Ctrl", ifelse(x == "ABA", "5uM_ABA", ifelse(x == "uMNa", "100uM_NaCl", "100mM_NaCl")))}) %>% unlist())
print(sprintf("AthEsaSirSpa_SCTIntegrated has %s genes x %s cells.", nrow(AthEsaSirSpa_SCTIntegrated), ncol(AthEsaSirSpa_SCTIntegrated)))


# Extract Arabidopsis data and add the revised annotation information.
AthEsaSirSpa_SCTIntegrated_AthCtrls <- subset(AthEsaSirSpa_SCTIntegrated, subset = sample == "AthCR1" | sample == "AthCR2")
print(sprintf("AthEsaSirSpa_SCTIntegrated_AthCtrls has %s genes x %s cells.", nrow(AthEsaSirSpa_SCTIntegrated_AthCtrls), ncol(AthEsaSirSpa_SCTIntegrated_AthCtrls)))
identlist <- list("1_1" = "1_1_1", "1_2" = "1_2_2")
AthEsaSirSpa_SCTIntegrated_AthCtrls@meta.data <- left_join(AthEsaSirSpa_SCTIntegrated_AthCtrls@meta.data %>% rownames_to_column(var = "cellid"), Ath_Ctrls_cCellTypes %>% rownames_to_column(var = "cellid") %>% mutate(cellid = gsubfn::gsubfn(paste(names(identlist), collapse = "|"), identlist, cellid)), by = "cellid") %>% column_to_rownames(var = "cellid")
print(table(AthEsaSirSpa_SCTIntegrated_AthCtrls@meta.data$cCellType))

# Predicting with different ndim and kweight.
for (ndim in seq(100, 500, 10)){
	for (kscore in c(10, 30, 50)){
		TransferAnchors <- FindTransferAnchors(reference = AthEsaSirSpa_SCTIntegrated_AthCtrls, query = AthEsaSirSpa_SCTIntegrated, normalization.method = "SCT", k.score = kscore, npcs = ndim, dims = 1:ndim)
		for (kweight in c(5, 10, 20, 30, 50, 80, 100)){
		tryCatch({
			TransferredLabels <- TransferData(anchorset = TransferAnchors, refdata = AthEsaSirSpa_SCTIntegrated_AthCtrls$cCellType, k.weight = kweight, dims = 1:ndim)
			saveRDS(TransferredLabels, file = paste0("AthEsaSirSpa_SCTPredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS")) 
			}, error = function(e){print(sprintf("ERROR at ndim %s, kscore %s and kweight %s: %s", ndim, kscore, kweight, conditionMessage(e)))}
			)
		}
	}
	
}


#***************************************************************************************************************************************************#
#***************************************************************************************************************************************************#
## get concensus annotation from annotation transfer
AnnotationList <- list()
for (ndim in c(30, seq(50, 500, 50))){
	for (kscore in c(10, 30, 50)){
		for (kweight in c(10, 20, 30, 50, 80, 100)){
			if (file.exists(paste0("AthEsaSirSpa_", string, "PredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS"))){
				AnnotationList[[paste0("ndim", ndim, "kscore", kscore, "kweight", kweight)]] <- readRDS(paste0("AthEsaSirSpa_", string, "PredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS")) %>% select(predicted.id, prediction.score.max) %>% rename_with( ~ paste0("ndim", ndim, "kscore", kscore, "kweight", kweight, "_", .x)) %>% rownames_to_column(var = "cellid")
			} else {
				print(sprintf("Cannot find %s...", paste0("AthEsaSirSpa_", string, "PredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS")))
			}
		}
	}    
}
Annotations <- Reduce(function(x,y) {dplyr::left_join(x, y, by = "cellid")}, AnnotationList) %>% column_to_rownames(var = "cellid")


Annotations$consensusann <- sapply(1:nrow(Annotations), function (i) {
	anns <- Annotations[i, ] %>% select(contains("_predicted.id")) %>% as.character() %>% table() %>% sort(., decreasing = TRUE)
	if (length(anns) > 1) {
		first <- anns[1]
		second <- anns[2]
		if (as.numeric(first) == as.numeric(second)) {
			### When the top two annotations have equal support, select the one with higher probability
			firstprob <- mean(as.numeric(Annotations[i, ][which(Annotations[i, ] == names(first)) + 1]))
			secondprob <- mean(as.numeric(Annotations[i, ][which(Annotations[i, ] == names(second)) + 1]))
			if (firstprob > secondprob) {
				names(first)
			} else if (firstprob < secondprob) {
				 names(second)
			} else {
				sample(c(names(first), names(second)), 1)
			}
		} else {
			names(first)
		}
	} else {
		names(anns[1])
	}
})

### Get the annotation frequency for each cell.
Annotations$annfrequency <- apply(Annotations, 1, function(ann) {sum(ann[1:(ncol(Annotations) - 1)] == ann[ncol(Annotations)])})
### Calculate the average annotation score for each cell.
Annotations$average.score <- sapply(1:nrow(Annotations), function(i) {mean(as.numeric(Annotations[i, ][which(Annotations[i, 1:(ncol(Annotations) - 2)] == Annotations[i, (ncol(Annotations) - 1)]) + 1]))})
### Check how many cells with the top two annotaions having equal support.
Annotations$ambiguousAnn <- apply(Annotations %>% select(contains("_predicted.id")), 1, function(x) {as.numeric(sort(table(as.character(x)), decreasing = TRUE))[1] == as.numeric(sort(table(as.character(x)), decreasing = TRUE))[2]}) %>% replace(., is.na(.), FALSE)
#######################################################################################################################################
#######################################################################################################################################




















#######################################################################################################################################
## 2. scArches
#######################################################################################################################################
## running scArches with different parameters, see ann4AthEsaSirSpa_scArches.py


#***************************************************************************************************************************************************#
#***************************************************************************************************************************************************#
## get concensus annotation from annotation transfer
AnnotationList <- list()
for (hvg_num in seq(2000, 14000, 2000)) {
	if (file.exists(paste0("AthEsaSirSpa_", hvg_num,"HVGs_scArchesAnnotaion.txt"))){
		print(sprintf("Processing %s...", paste0("AthEsaSirSpa_", hvg_num,"HVGs_scArchesAnnotaion.txt")))
		AnnotationList[[paste0("topn", hvg_num)]] <- read.csv(paste0("AthEsaSirSpa_", hvg_num,"HVGs_scArchesAnnotaion.txt"), sep = "\t", header = TRUE) %>% separate(X, c("cellid", "batchid"), "-") %>% mutate(cellid = gsub("\\.", "-", cellid)) %>% column_to_rownames(var = "cellid") %>% select(predictions) %>% rename_with( ~ paste0("topn", hvg_num, "_", .x)) %>% rownames_to_column(var = "cellid")
	}
}

scArchesAnnotations <- Reduce(function(x,y) {dplyr::left_join(x, y, by = "cellid")}, AnnotationList) %>% column_to_rownames(var = "cellid")


scArchesAnnotations$consensusann <- sapply(1:nrow(scArchesAnnotations), function (i) {
	anns <- scArchesAnnotations[i, ] %>% select(contains("_predictions")) %>% as.character() %>% table() %>% sort(., decreasing = TRUE)
	if (length(anns) > 1) {
		first <- anns[1]
		second <- anns[2]
		if (as.numeric(first) == as.numeric(second)) {
			### When the top two annotations have equal support, select the one with higher probability
				sample(c(names(first), names(second)), 1)
		} else {
			names(first)
		}
	} else {
		names(anns[1])
	}
})

### Get the annotation frequency for each cell.
scArchesAnnotations$annfrequency <- apply(scArchesAnnotations, 1, function(ann) {sum(ann[1:(ncol(scArchesAnnotations) - 1)] == ann[ncol(scArchesAnnotations)])})
### Check how many cells with the top two annotaions having equal support.
scArchesAnnotations$ambiguousAnn <- apply(scArchesAnnotations %>% select(contains("_predictions")), 1, function(x) {as.numeric(sort(table(as.character(x)), decreasing = TRUE))[1] == as.numeric(sort(table(as.character(x)), decreasing = TRUE))[2]}) %>% replace(., is.na(.), FALSE)
#######################################################################################################################################
#######################################################################################################################################



















#######################################################################################################################################
## 3. Symphony
#######################################################################################################################################
### annotation transfer with different parameters
AthEsaSirSpa <- readRDS("integratedAthEsaSirSpa.RDS")
AthEsaSirSpa@meta.data <- AthEsaSirSpa@meta.data %>% mutate(species = str_sub(orig.ident, 14, 16), treatment = str_sub(orig.ident, 17, -3), replicate = str_sub(orig.ident, -2, -1))

AthEsaSirSpa_AthCtrls <- subset(AthEsaSirSpa, subset = orig.ident %in% c("AthEsaSirSpa_AthCR1", "AthEsaSirSpa_AthCR2"))
AthEsaSirSpa_nonAthCtrls <- subset(AthEsaSirSpa, subset = orig.ident %in% c("AthEsaSirSpa_AthCR1", "AthEsaSirSpa_AthCR2"), invert = TRUE)

## Extract Arabidopsis data and add the revised annotation information.
identlist <- list("1_1" = "1_1_1", "1_2" = "1_2_2")
Ath_Ctrls_cCellTypes <- readRDS("Ath_Ctrls_cCellTypes.RDS") ### Ath_Ctrls_cCellTypes.RDS is generated in consensusann4Ath.R
AthEsaSirSpa_AthCtrls@meta.data <- left_join(AthEsaSirSpa_AthCtrls@meta.data %>% rownames_to_column(var = "cellid"), Ath_Ctrls_cCellTypes %>% select(cCellType) %>% rownames_to_column(var = "cellid") %>% mutate(cellid = gsubfn::gsubfn(paste(names(identlist), collapse = "|"), identlist, cellid)), by = "cellid") %>% column_to_rownames(var = "cellid") %>% as.data.frame() %>% replace(., is.na(.), "")
print(table(AthEsaSirSpa_AthCtrls@meta.data$cCellType))


AthEsaSirSpa_AthCtrls_NorExp <- AthEsaSirSpa_AthCtrls@assays$RNA@data
AthEsaSirSpa_AthCtrls_MetaData <- AthEsaSirSpa_AthCtrls@meta.data
AthEsaSirSpa_nonAthCtrls_NorExp <- AthEsaSirSpa_nonAthCtrls@assays$RNA@data
AthEsaSirSpa_nonAthCtrls_MetaData <- AthEsaSirSpa_nonAthCtrls@meta.data

print(sprintf("The reference dataset has %s features x %s cells with %s cells x %s attributes in meta.data; the query dataset has %s features x %s cells with %s cells x %s attributes in meta.data.", nrow(AthEsaSirSpa_AthCtrls_NorExp), ncol(AthEsaSirSpa_AthCtrls_NorExp), nrow(AthEsaSirSpa_AthCtrls_MetaData), ncol(AthEsaSirSpa_AthCtrls_MetaData), nrow(AthEsaSirSpa_nonAthCtrls_NorExp), ncol(AthEsaSirSpa_nonAthCtrls_NorExp), nrow(AthEsaSirSpa_nonAthCtrls_MetaData), ncol(AthEsaSirSpa_nonAthCtrls_MetaData)))


for (hvg_num in seq(2000, 14000, 2000)) {
	for (pc_num in seq(20, 100, 20)) {
	# Build reference
		set.seed(0)
		reference = symphony::buildReference(
			AthEsaSirSpa_AthCtrls_NorExp,
			AthEsaSirSpa_AthCtrls_MetaData,
			vars = c("orig.ident"),         # variables to integrate over
			K = 100,                   # number of Harmony clusters
			verbose = TRUE,            # verbose output
			do_umap = TRUE,            # can set to FALSE if want to run umap separately later
			do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
			vargenes_method = "vst",   # method for variable gene selection ('vst' or 'mvp')
			vargenes_groups = "orig.ident", # metadata column specifying groups for variable gene selection 
			topn = hvg_num,               # number of variable genes to choose per group
			d = pc_num,                    # number of PCs
			save_uwot_path = paste0("220710_AthEsaSirSpa_topn", hvg_num, "d", pc_num, "_ref_model")
		)
		reference$normalization_method = "log(CP10k+1)" # optionally save normalization method in custom slot
		# Save reference (modify with your desired output path)
		saveRDS(reference, paste0("220710_AthEsaSirSpa_topn", hvg_num, "d", pc_num, "_symphony_ref.rds"))
		
		
		# Map query
		query <- symphony::mapQuery(AthEsaSirSpa_nonAthCtrls_NorExp,             # query gene expression (genes x cells)
			AthEsaSirSpa_nonAthCtrls_MetaData,        # query metadata (cells x attributes)
			reference,             # Symphony reference object
			vars = c("orig.ident"), # Query batch variable(s) to integrate over (column names in metadata)
			do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
			do_umap = TRUE)        # project query cells into reference UMAP
		
		
		for (k_num in c(5, 10, 20, 30, 50, 80, 100)) {
			query <- symphony::knnPredict(query, reference, reference$meta_data$cCellType, k = 5) # For k-NN prediction, we would recommend that users alter the k parameter so that it is ideally no larger than the number of cells in the rarest cell type of the reference. For example, if the reference contains only 10 cells of a rare cell type, then we recommend the user set k no higher than 10, to ensure that rare cell types in the reference have the chance of being predicted given a majority vote k-NN classifier.
			
			write.table(query$meta_data, file = paste0("AthEsaSirSpa_topn", hvg_num, "d", pc_num, "k", k_num, "_SymphonyAnn.csv"), sep = ",", row.names = T, col.names = T, quote = F)

		}
	}
}

#***************************************************************************************************************************************************#
#***************************************************************************************************************************************************#
## get concensus annotation from annotation transfer
AnnotationList <- list()
for (hvg_num in seq(2000, 14000, 2000)) {
    for (pc_num in seq(20, 100, 20)) {
        for (k_num in c(5, 10, 20, 30, 50, 80, 100)) {
            if (file.exists(paste0("AthEsaSirSpa_topn", hvg_num, "d", pc_num, "k", k_num, "_SymphonyAnn.csv"))){
                 AnnotationList[[paste0("topn", hvg_num, "d", pc_num, "k", k_num)]] <- read.csv(paste0("AthEsaSirSpa_topn", hvg_num, "d", pc_num, "k", k_num, "_SymphonyAnn.csv"), header = TRUE) %>% select(cell_type_pred_knn, cell_type_pred_knn_prob) %>% rename(c("pred.celltype" = "cell_type_pred_knn", "pred.prob" = "cell_type_pred_knn_prob")) %>% rename_with( ~ paste0("topn", hvg_num, "d", pc_num, "k", k_num, "_", .x)) %>% rownames_to_column(var = "cellid")
            }
        }
    }
}

SymphonyAnnotations <- Reduce(function(x,y) {dplyr::left_join(x, y, by = "cellid")}, AnnotationList) %>% column_to_rownames(var = "cellid")

SymphonyAnnotations$consensusann <- sapply(1:nrow(SymphonyAnnotations), function (i) {
	anns <- SymphonyAnnotations[i, ] %>% select(contains("_pred.celltype")) %>% as.character() %>% table() %>% sort(., decreasing = TRUE)
	if (length(anns) > 1) {
		first <- anns[1]
		second <- anns[2]
		if (as.numeric(first) == as.numeric(second)) {
			### When the top two annotations have equal support, select the one with higher probability
			firstprob <- mean(as.numeric(SymphonyAnnotations[i, ][which(SymphonyAnnotations[i, ] == names(first)) + 1]))
			secondprob <- mean(as.numeric(SymphonyAnnotations[i, ][which(SymphonyAnnotations[i, ] == names(second)) + 1]))
			if (firstprob > secondprob) {
				names(first)
			} else if (firstprob < secondprob) {
				 names(second)
			} else {
				sample(c(names(first), names(second)), 1)
			}
		} else {
			names(first)
		}
	} else {
		names(anns[1])
	}
})

### Get the annotation frequency for each cell.
SymphonyAnnotations$annfrequency <- apply(SymphonyAnnotations, 1, function(ann) {sum(ann[1:(ncol(SymphonyAnnotations) - 1)] == ann[ncol(SymphonyAnnotations)])})
### Calculate the average annotation score for each cell.
SymphonyAnnotations$average.score <- sapply(1:nrow(SymphonyAnnotations), function(i) {mean(as.numeric(SymphonyAnnotations[i, ][which(SymphonyAnnotations[i, 1:(ncol(SymphonyAnnotations) - 2)] == SymphonyAnnotations[i, (ncol(SymphonyAnnotations) - 1)]) + 1]))})
### Check how many cells with the top two annotaions having equal support.
SymphonyAnnotations$ambiguousAnn <- apply(SymphonyAnnotations %>% select(contains("_pred.celltype")), 1, function(x) {as.numeric(sort(table(as.character(x)), decreasing = TRUE))[1] == as.numeric(sort(table(as.character(x)), decreasing = TRUE))[2]}) %>% replace(., is.na(.), FALSE)
#######################################################################################################################################
#######################################################################################################################################












#######################################################################################################################################
## 4. find the consensus annotation
#######################################################################################################################################
### read the integrated data
AthEsaSirSpa_SCTIntegrated <- readRDS("integratedAthEsaSirSpa.RDS")
AthEsaSirSpa_MetaData <- AthEsaSirSpa_SCTIntegrated@meta.data %>% mutate(sample = gsub("AthEsaSirSpa_", "", orig.ident), species = str_sub(orig.ident, 14, 16), treatment = str_sub(orig.ident, 17, -3), replicate = str_sub(orig.ident, -2, -1)) %>% mutate(treatment = lapply(treatment, function(x) {ifelse(x == "C", "Ctrl", ifelse(x == "ABA", "5uM_ABA", ifelse(x == "uMNa", "100uM_NaCl", "100mM_NaCl")))}) %>% unlist())
print(sprintf("AthEsaSirSpa_MetaData has %s cells and %s annotations.", nrow(AthEsaSirSpa_MetaData), ncol(AthEsaSirSpa_MetaData)))


Ath_Ctrls_cCellTypes <- readRDS("Ath_Ctrls_cCellTypes.RDS") ### Ath_Ctrls_cCellTypes.RDS is generated in consensusann4Ath.R

identlist <- list("1_1" = "1_1_1", "1_2" = "1_2_2")
AthEsaSirSpa_MetaData <- left_join(AthEsaSirSpa_MetaData %>% rownames_to_column(var = "cellid"), Ath_Ctrls_cCellTypes %>% select(cCellType) %>% rownames_to_column(var = "cellid") %>% mutate(cellid = gsubfn::gsubfn(paste(names(identlist), collapse = "|"), identlist, cellid)), by = "cellid") %>% column_to_rownames(var = "cellid") 

### read annotations from the three different methods
AthEsaSirSpa_TransferConsensusAnn <- Annotations
AthEsaSirSpa_SymphonyConsensusAnns <- scArchesAnnotations
AthEsaSirSpa_scArchesConsensusAnns <- SymphonyAnnotations

print(sprintf("AthEsaSirSpa_TransferConsensusAnn has %s rows and %s columns; AthEsaSirSpa_SymphonyConsensusAnns has %s rows and %s columns; AthEsaSirSpa_scArchesConsensusAnns has %s rows and %s columns.", nrow(AthEsaSirSpa_TransferConsensusAnn), ncol(AthEsaSirSpa_TransferConsensusAnn), nrow(AthEsaSirSpa_SymphonyConsensusAnns), ncol(AthEsaSirSpa_SymphonyConsensusAnns), nrow(AthEsaSirSpa_scArchesConsensusAnns), ncol(AthEsaSirSpa_scArchesConsensusAnns)))



AthEsaSirSpa_Annotations <- left_join(AthEsaSirSpa_MetaData %>% rownames_to_column(var = "cellid"), AthEsaSirSpa_TransferConsensusAnn %>% select(consensusann, annfrequency, average.score, ambiguousAnn) %>% rename_with(~ paste0("standard_", .x)) %>% rownames_to_column(var = "cellid"), by = "cellid") %>% left_join(., AthEsaSirSpa_SymphonyConsensusAnns %>% select(consensusann, annfrequency, average.score, ambiguousAnn) %>% rename_with(~ paste0("symphony_", .x)) %>% rownames_to_column(var = "cellid"), by = "cellid") %>% left_join(., AthEsaSirSpa_scArchesConsensusAnns %>% select(consensusann, annfrequency, ambiguousAnn) %>% rename_with(~ paste0("scarches_", .x)) %>% rownames_to_column(var = "cellid"), by = "cellid") %>% column_to_rownames(var = "cellid")


AthEsaSirSpa_Annotations$celltype <- sapply(1:nrow(AthEsaSirSpa_Annotations), function (i) {
    if (grepl("AthCR", AthEsaSirSpa_Annotations[i, "sample"], fixed = TRUE)) {
        return(as.character(AthEsaSirSpa_Annotations[i, "cCellType"]))
    } else {
        anns <- AthEsaSirSpa_Annotations[i, ] %>% select(-contains("SCT_")) %>% select(contains("_consensusann")) %>% as.character() %>% table() %>% sort(., decreasing = TRUE)
        if (length(anns) > 2) {
            ### When the top two annotations have equal support, select the one with higher probability
            if (as.numeric(AthEsaSirSpa_Annotations[i, "standard_average.score"]) > as.numeric(AthEsaSirSpa_Annotations[i, "symphony_average.score"])) {
                return(as.character(AthEsaSirSpa_Annotations[i, "standard_consensusann"]))
            } else if (as.numeric(AthEsaSirSpa_Annotations[i, "standard_average.score"]) < as.numeric(AthEsaSirSpa_Annotations[i, "symphony_average.score"])) {
                return(as.character(AthEsaSirSpa_Annotations[i, "symphony_consensusann"]))
            } else {
                sample(c(as.character(AthEsaSirSpa_Annotations[i, "standard_consensusann"]), as.character(AthEsaSirSpa_Annotations[i, "symphony_consensusann"])), 1)
            }
        } else {
            names(anns[1])
        }
    }
})
#######################################################################################################################################
#######################################################################################################################################


print("Done")