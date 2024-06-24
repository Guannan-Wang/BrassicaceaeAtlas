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


AthCsa_SCTIntegrated <- readRDS("integratedAthCsa.RDS")
AthCsa_SCTIntegrated@meta.data <- AthCsa_SCTIntegrated@meta.data %>% mutate(sample = gsub("AthCsa_", "", orig.ident)) %>% mutate(species = str_sub(sample, 1, 3), treatment = str_sub(sample, 4, nchar(sample)-2), replicate = str_sub(sample, -2, -1)) %>% mutate(treatment = lapply(treatment, function(x) {ifelse(x == "C", "Ctrl", ifelse(x == "ABA", "5uM_ABA", ifelse(x == "uMNa", "100uM_NaCl", "100mM_NaCl")))}) %>% unlist())
print(sprintf("AthCsa_SCTIntegrated has %s genes x %s cells.", nrow(AthCsa_SCTIntegrated), ncol(AthCsa_SCTIntegrated)))


# Extract Arabidopsis data and add the revised annotation information.
AthCsa_SCTIntegrated_AthCtrls <- subset(AthCsa_SCTIntegrated, subset = sample == "AthCR1" | sample == "AthCR2")
print(sprintf("AthCsa_SCTIntegrated_AthCtrls has %s genes x %s cells.", nrow(AthCsa_SCTIntegrated_AthCtrls), ncol(AthCsa_SCTIntegrated_AthCtrls)))
identlist <- list("1_1" = "1_1_1", "1_2" = "1_2_2")
AthCsa_SCTIntegrated_AthCtrls@meta.data <- left_join(AthCsa_SCTIntegrated_AthCtrls@meta.data %>% rownames_to_column(var = "cellid"), Ath_Ctrls_cCellTypes %>% rownames_to_column(var = "cellid") %>% mutate(cellid = gsubfn::gsubfn(paste(names(identlist), collapse = "|"), identlist, cellid)), by = "cellid") %>% column_to_rownames(var = "cellid")
print(table(AthCsa_SCTIntegrated_AthCtrls@meta.data$cCellType))


for (ndim in c(30, seq(50, 500, 50))){
	for (kscore in c(10, 30, 50)){
		TransferAnchors <- FindTransferAnchors(reference = AthCsa_SCTIntegrated_AthCtrls, query = AthCsa_SCTIntegrated, reduction = "pcaproject", k.score = kscore, npcs = ndim, dims = 1:ndim)
		for (kweight in c(10, 20, 30, 50, 80, 100)){
		tryCatch({
			TransferredLabels <- TransferData(anchorset = TransferAnchors, refdata = AthCsa_SCTIntegrated_AthCtrls$cCellType, k.weight = kweight, dims = 1:ndim) 
			saveRDS(TransferredLabels, file = paste0("AthCsa_PredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS"))
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
			if (file.exists(paste0("AthCsa_", string, "PredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS"))){
				AnnotationList[[paste0("ndim", ndim, "kscore", kscore, "kweight", kweight)]] <- readRDS(paste0("AthCsa_", string, "PredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS")) %>% select(predicted.id, prediction.score.max) %>% rename_with( ~ paste0("ndim", ndim, "kscore", kscore, "kweight", kweight, "_", .x)) %>% rownames_to_column(var = "cellid")
			} else {
				print(sprintf("Cannot find %s...", paste0("AthCsa_", string, "PredictedID_ndim", ndim, "kscore", kscore, "kweight", kweight, ".RDS")))
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
print(Annotations$ambiguousAnn %>% table())
#######################################################################################################################################
#######################################################################################################################################




















#######################################################################################################################################
## 2. scArches
#######################################################################################################################################
## running scArches with different parameters, see ann4Csa_scArches.py


#***************************************************************************************************************************************************#
#***************************************************************************************************************************************************#
## get concensus annotation from annotation transfer
AnnotationList <- list()
for (hvg_num in seq(2000, 14000, 2000)) {
	if (file.exists(paste0("AthCsa_", hvg_num,"HVGs_scArchesAnnotaion.txt"))){
		print(sprintf("Processing %s...", paste0("AthCsa_", hvg_num,"HVGs_scArchesAnnotaion.txt")))
		AnnotationList[[paste0("topn", hvg_num)]] <- read.csv(paste0("AthCsa_", hvg_num,"HVGs_scArchesAnnotaion.txt"), sep = "\t", header = TRUE) %>% separate(X, c("cellid", "batchid"), "-") %>% mutate(cellid = gsub("\\.", "-", cellid)) %>% column_to_rownames(var = "cellid") %>% select(predictions) %>% rename_with( ~ paste0("topn", hvg_num, "_", .x)) %>% rownames_to_column(var = "cellid")
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
AthCsa <- readRDS("integratedAthCsa.RDS")
AthCsa@meta.data <- AthCsa@meta.data %>% mutate(species = str_sub(orig.ident, 8, 10), treatment = str_sub(orig.ident, 11, -3), replicate = str_sub(orig.ident, -2, -1))

AthCsa_AthCtrls <- subset(AthCsa, subset = orig.ident %in% c("AthCsa_AthCR1", "AthCsa_AthCR2"))
AthCsa_CsaCtrls <- subset(AthCsa, subset = orig.ident %in% c("AthCsa_AthCR1", "AthCsa_AthCR2"), invert = TRUE)

## Extract Arabidopsis data and add the revised annotation information.
identlist <- list("1_1" = "1_1_1", "1_2" = "1_2_2")
Ath_Ctrls_cCellTypes <- readRDS("Ath_Ctrls_cCellTypes.RDS") ### Ath_Ctrls_cCellTypes.RDS is generated in consensusann4Ath.R
AthCsa_AthCtrls@meta.data <- left_join(AthCsa_AthCtrls@meta.data %>% rownames_to_column(var = "cellid"), Ath_Ctrls_cCellTypes %>% select(cCellType) %>% rownames_to_column(var = "cellid") %>% mutate(cellid = gsubfn::gsubfn(paste(names(identlist), collapse = "|"), identlist, cellid)), by = "cellid") %>% column_to_rownames(var = "cellid") %>% as.data.frame() %>% replace(., is.na(.), "")
print(table(AthCsa_AthCtrls@meta.data$cCellType))


AthCsa_AthCtrls_NorExp <- AthCsa_AthCtrls@assays$RNA@data
AthCsa_AthCtrls_MetaData <- AthCsa_AthCtrls@meta.data
AthCsa_CsaCtrls_NorExp <- AthCsa_CsaCtrls@assays$RNA@data
AthCsa_CsaCtrls_MetaData <- AthCsa_CsaCtrls@meta.data

print(sprintf("The reference dataset has %s features x %s cells with %s cells x %s attributes in meta.data; the query dataset has %s features x %s cells with %s cells x %s attributes in meta.data.", nrow(AthCsa_AthCtrls_NorExp), ncol(AthCsa_AthCtrls_NorExp), nrow(AthCsa_AthCtrls_MetaData), ncol(AthCsa_AthCtrls_MetaData), nrow(AthCsa_CsaCtrls_NorExp), ncol(AthCsa_CsaCtrls_NorExp), nrow(AthCsa_CsaCtrls_MetaData), ncol(AthCsa_CsaCtrls_MetaData)))


for (hvg_num in seq(2000, 14000, 2000)) {
	for (pc_num in seq(20, 100, 20)) {
	# Build reference
		set.seed(0)
		reference = symphony::buildReference(
			AthCsa_AthCtrls_NorExp,
			AthCsa_AthCtrls_MetaData,
			vars = c("orig.ident"),         # variables to integrate over
			K = 100,                   # number of Harmony clusters
			verbose = TRUE,            # verbose output
			do_umap = TRUE,            # can set to FALSE if want to run umap separately later
			do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
			vargenes_method = "vst",   # method for variable gene selection ('vst' or 'mvp')
			vargenes_groups = "orig.ident", # metadata column specifying groups for variable gene selection 
			topn = hvg_num,               # number of variable genes to choose per group
			d = pc_num,                    # number of PCs
			save_uwot_path = paste0("220805_AthCsa_topn", hvg_num, "d", pc_num, "_ref_model")
		)
		reference$normalization_method = "log(CP10k+1)" # optionally save normalization method in custom slot
		# Save reference (modify with your desired output path)
		saveRDS(reference, paste0("220805_AthCsa_topn", hvg_num, "d", pc_num, "_symphony_ref.rds"))
		
		
		# Map query
		query <- symphony::mapQuery(AthCsa_CsaCtrls_NorExp,             # query gene expression (genes x cells)
			AthCsa_CsaCtrls_MetaData,        # query metadata (cells x attributes)
			reference,             # Symphony reference object
			vars = c("orig.ident"), # Query batch variable(s) to integrate over (column names in metadata)
			do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
			do_umap = TRUE)        # project query cells into reference UMAP
		
		
		for (k_num in c(5, 10, 20, 30, 50, 80, 100)) {
			query <- symphony::knnPredict(query, reference, reference$meta_data$cCellType, k = 5) # For k-NN prediction, we would recommend that users alter the k parameter so that it is ideally no larger than the number of cells in the rarest cell type of the reference. For example, if the reference contains only 10 cells of a rare cell type, then we recommend the user set k no higher than 10, to ensure that rare cell types in the reference have the chance of being predicted given a majority vote k-NN classifier.
			
			write.table(query$meta_data, file = paste0("AthCsa_topn", hvg_num, "d", pc_num, "k", k_num, "_SymphonyAnn.csv"), sep = ",", row.names = T, col.names = T, quote = F)

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
            if (file.exists(paste0("AthCsa_topn", hvg_num, "d", pc_num, "k", k_num, "_SymphonyAnn.csv"))){
                 AnnotationList[[paste0("topn", hvg_num, "d", pc_num, "k", k_num)]] <- read.csv(paste0("AthCsa_topn", hvg_num, "d", pc_num, "k", k_num, "_SymphonyAnn.csv"), header = TRUE) %>% select(cell_type_pred_knn, cell_type_pred_knn_prob) %>% rename(c("pred.celltype" = "cell_type_pred_knn", "pred.prob" = "cell_type_pred_knn_prob")) %>% rename_with( ~ paste0("topn", hvg_num, "d", pc_num, "k", k_num, "_", .x)) %>% rownames_to_column(var = "cellid")
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
AthCsa_SCTIntegrated <- readRDS("integratedAthCsa.RDS")
AthCsa_MetaData <- AthCsa_SCTIntegrated@meta.data %>% mutate(sample = gsub("AthCsa_", "", orig.ident), species = str_sub(orig.ident, 8, 10), treatment = str_sub(orig.ident, 11, -3), replicate = str_sub(orig.ident, -2, -1)) %>% mutate(treatment = lapply(treatment, function(x) {ifelse(x == "C", "Ctrl", ifelse(x == "ABA", "5uM_ABA", ifelse(x == "uMNa", "100uM_NaCl", "100mM_NaCl")))}) %>% unlist())
print(sprintf("AthCsa_MetaData has %s cells and %s annotations.", nrow(AthCsa_MetaData), ncol(AthCsa_MetaData)))


Ath_Ctrls_cCellTypes <- readRDS("Ath_Ctrls_cCellTypes.RDS") ### Ath_Ctrls_cCellTypes.RDS is generated in consensusann4Ath.R

identlist <- list("1_1" = "1_1_1", "1_2" = "1_2_2")
AthCsa_MetaData <- left_join(AthCsa_MetaData %>% rownames_to_column(var = "cellid"), Ath_Ctrls_cCellTypes %>% select(cCellType) %>% rownames_to_column(var = "cellid") %>% mutate(cellid = gsubfn::gsubfn(paste(names(identlist), collapse = "|"), identlist, cellid)), by = "cellid") %>% column_to_rownames(var = "cellid") 

### read annotations from the three different methods
AthCsa_TransferConsensusAnn <- Annotations
AthCsa_SymphonyConsensusAnns <- scArchesAnnotations
AthCsa_scArchesConsensusAnns <- SymphonyAnnotations

print(sprintf("AthCsa_SeuratConsensusAnn has %s rows and %s columns; AthCsa_SymphonyConsensusAnns has %s rows and %s columns; AthCsa_scArchesConsensusAnns has %s rows and %s columns.", nrow(AthCsa_SeuratConsensusAnn), ncol(AthCsa_SeuratConsensusAnn), nrow(AthCsa_SymphonyConsensusAnns), ncol(AthCsa_SymphonyConsensusAnns), nrow(AthCsa_scArchesConsensusAnns), ncol(AthCsa_scArchesConsensusAnns)))



AthCsa_Annotations <- left_join(AthCsa_MetaData %>% rownames_to_column(var = "cellid"), AthCsa_SeuratConsensusAnn %>% select(consensusann, annfrequency, average.score, ambiguousAnn) %>% rename_with(~ paste0("seurat_", .x)) %>% rownames_to_column(var = "cellid"), by = "cellid") %>% left_join(., AthCsa_SymphonyConsensusAnns %>% select(consensusann, annfrequency, average.score, ambiguousAnn) %>% rename_with(~ paste0("symphony_", .x)) %>% rownames_to_column(var = "cellid"), by = "cellid") %>% left_join(., AthCsa_scArchesConsensusAnns %>% select(consensusann, annfrequency, ambiguousAnn) %>% rename_with(~ paste0("scarches_", .x)) %>% rownames_to_column(var = "cellid"), by = "cellid") %>% column_to_rownames(var = "cellid")



AthCsa_Annotations$celltype <- sapply(1:nrow(AthCsa_Annotations), function (i) {
    if (grepl("AthCR", AthCsa_Annotations[i, "sample"], fixed = TRUE)) {
        return(as.character(AthCsa_Annotations[i, "cCellType"]))
    } else {
        anns <- AthCsa_Annotations[i, ] %>% select(-contains("SCT_")) %>% select(contains("_consensusann")) %>% as.character() %>% table() %>% sort(., decreasing = TRUE)
        if (length(anns) > 2) {
            ### When the top two annotations have equal support, select the one with higher probability
            if (as.numeric(AthCsa_Annotations[i, "seurat_average.score"]) > as.numeric(AthCsa_Annotations[i, "symphony_average.score"])) {
                return(as.character(AthCsa_Annotations[i, "seurat_consensusann"]))
            } else if (as.numeric(AthCsa_Annotations[i, "seurat_average.score"]) < as.numeric(AthCsa_Annotations[i, "symphony_average.score"])) {
                return(as.character(AthCsa_Annotations[i, "symphony_consensusann"]))
            } else {
                sample(c(as.character(AthCsa_Annotations[i, "seurat_consensusann"]), as.character(AthCsa_Annotations[i, "symphony_consensusann"])), 1)
            }
        } else {
            names(anns[1])
        }
    }
})
#######################################################################################################################################
#######################################################################################################################################


print("Done")