#!/usr/bin/env Rscript

load_libraries <- c("Matrix", "Seurat", "tidyverse", "readxl", "Cairo", "SeuratWrappers", "clustree")
invisible(lapply(load_libraries, function(x){suppressMessages(library(x, character.only = TRUE, quietly = TRUE))}))

#######################################################################################################################################
## 1. correlation based method
#######################################################################################################################################
### calculate correlation to reference transcriptomes
for (mergedtable in list.files("/path/to/output/dir", pattern = "ann4Ath_Ath_Ctrl_Combined.*.txt", full.names = TRUE)){
	if (grepl("_Cor.txt", mergedtable, fixed = TRUE)){
		cat(basename(mergedtable), "\t", "Skip this file!", "\n", sep = "")
	} else {
		mergedtable_temp <- read.delim(mergedtable, header = TRUE)
		if (grepl("wBrady", mergedtable, fixed = TRUE)){
			mergedtable_cellcor_temp <- suppressWarnings(sapply(130:ncol(mergedtable_temp), function(i) sapply(2:129, function(j) cor.test(mergedtable_temp[, i], mergedtable_temp[, j], method = "pearson")[c("p.value", "estimate")]))) ##Each column represents a cell; each row represents a reference cell type.
			assign(paste(strsplit(basename(mergedtable), "ann4Ath_|.txt")[[1]][2], "_CorTable",sep = ""), mergedtable_cellcor_temp)
		
			mergedtable_cellcor_temp_cor <- mergedtable_cellcor_temp[seq(2,nrow(mergedtable_cellcor_temp),2),]
			mergedtable_cellcor_temp_maxcor <- sapply(1:ncol(mergedtable_cellcor_temp_cor), function(i) max(as.numeric(mergedtable_cellcor_temp_cor[,i])))
			mergedtable_cellcor_temp_pvalue <- mergedtable_cellcor_temp[seq(1,nrow(mergedtable_cellcor_temp)-1,2),]

			mergedtable_cellcor_temp_cellident <- sapply(1:ncol(mergedtable_cellcor_temp_cor), function(i) colnames(mergedtable_temp)[2:129][which(as.numeric(mergedtable_cellcor_temp_cor[,i]) == max(as.numeric(mergedtable_cellcor_temp_cor[,i])))])

			mergedtable_cellcor_temp_maxp <- sapply(1:(ncol(mergedtable_temp)-129), function(i) as.numeric(mergedtable_cellcor_temp_pvalue[,i])[which(as.numeric(mergedtable_cellcor_temp_cor[,i]) == max(as.numeric(mergedtable_cellcor_temp_cor[,i])))])
		
			mergedtable_cellcortable_temp <- data.frame(row.names = gsub("\\.", "-", colnames(mergedtable_temp)[130:ncol(mergedtable_temp)]), maxcor = as.character(mergedtable_cellcor_temp_maxcor), maxp = as.character(mergedtable_cellcor_temp_maxp), cellident = as.character(mergedtable_cellcor_temp_cellident))

			assign(paste(strsplit(basename(mergedtable), "ann4Ath_|.txt")[[1]][2], "_Cor",sep = ""), mergedtable_cellcortable_temp)
			write.table(mergedtable_cellcortable_temp, file = paste("ann4Ath_", strsplit(basename(mergedtable), "ann4Ath_|.txt")[[1]][2], "_ppmCor.txt",sep = ""), quote = FALSE, sep = "\t")
		
		} else if (grepl("wLi", mergedtable, fixed = TRUE)) {
			mergedtable_cellcor_temp <- suppressWarnings(sapply(17:ncol(mergedtable_temp), function(i) sapply(2:16, function(j) cor.test(mergedtable_temp[, i], mergedtable_temp[, j], method = "pearson")[c("p.value", "estimate")]))) ##Each column represents a cell; each row represents a reference cell type.
			assign(paste(strsplit(basename(mergedtable), "ann4Ath_|.txt")[[1]][2], "_CorTable",sep = ""), mergedtable_cellcor_temp)
			
			mergedtable_cellcor_temp_cor <- mergedtable_cellcor_temp[seq(2,nrow(mergedtable_cellcor_temp),2),]
			mergedtable_cellcor_temp_maxcor <- sapply(1:ncol(mergedtable_cellcor_temp_cor), function(i) max(as.numeric(mergedtable_cellcor_temp_cor[,i])))
			mergedtable_cellcor_temp_pvalue <- mergedtable_cellcor_temp[seq(1,nrow(mergedtable_cellcor_temp)-1,2),]

			mergedtable_cellcor_temp_cellident <- sapply(1:ncol(mergedtable_cellcor_temp_cor), function(i) colnames(mergedtable_temp)[2:16][which(as.numeric(mergedtable_cellcor_temp_cor[,i]) == max(as.numeric(mergedtable_cellcor_temp_cor[,i])))])

			mergedtable_cellcor_temp_maxp <- sapply(1:(ncol(mergedtable_temp)-16), function(i) as.numeric(mergedtable_cellcor_temp_pvalue[,i])[which(as.numeric(mergedtable_cellcor_temp_cor[,i]) == max(as.numeric(mergedtable_cellcor_temp_cor[,i])))])
		
			mergedtable_cellcortable_temp <- data.frame(row.names = gsub("\\.", "-", colnames(mergedtable_temp)[17:ncol(mergedtable_temp)]), maxcor = as.character(mergedtable_cellcor_temp_maxcor), maxp = as.character(mergedtable_cellcor_temp_maxp), cellident = as.character(mergedtable_cellcor_temp_cellident))
		
			assign(paste(strsplit(basename(mergedtable), "ann4Ath_|.txt")[[1]][2], "_Cor",sep = ""), mergedtable_cellcortable_temp)
			write.table(mergedtable_cellcortable_temp, file = paste("ann4Ath_", strsplit(basename(mergedtable), "ann4Ath_|.txt")[[1]][2], "_ppmCor.txt",sep = ""), quote = FALSE, sep = "\t")
			}
	}
}

#***************************************************************************************************************************************************#
#***************************************************************************************************************************************************#
## get concensus annotation from correlation based method
mergeannfiles <- function(dp, reg = NULL, nameext = NULL, sep = "\t"){
## mergeannfiles is intended for merging selected annotation files in a directory, the first column of the each file will be used for merging.
## dp (required), path to the directory
## reg (optional), regular expression to find matching file names, default: NULL
## nameext (optional), a named vector of extensions to be added to the colnames in each of the files that will be merged to avoid identical colnames, default: NULL (setNames(, ))
## sep (optional), the separator in all the files to be merged
    filelist <- list()
    n = 1
    for (file in list.files(path = dp, pattern = reg, full.names = TRUE)){
        if (grepl("wBrady", file, fixed = TRUE)) {
            filelist[[n]] <- read.table(file, sep = "\t", header = TRUE, row.names = 1) %>% mutate(cellident = gsub("_.*", "", cellident)) %>% setNames(paste0(nameext[basename(file)], "_", names(.))) %>% rownames_to_column(var = "matchid") #strsplit(cellident, "_")[[1]][1]
        } else if (grepl("wLi", file, fixed = TRUE)) {
            filelist[[n]] <- read.table(file, sep = "\t", header = TRUE, row.names = 1) %>% setNames(paste0(nameext[basename(file)], "_", names(.))) %>% rownames_to_column(var = "matchid")
        }
        n = n + 1
    }
    mergedfile <- Reduce(function(x, y, ...) {left_join(x, y, by = "matchid", ...)}, filelist)
    return(mergedfile)
}

AthCor <- "/path/to/correlation/files"


Ath_Ctrls_Cor_CellType <- mergeannfiles(AthCor, reg = "ann4Ath_Ath_Ctrl_Combined_.*_top.*_Cor.txt", nameext = sapply(sapply(list.files(AthCor, pattern = "ann4Ath_Ath_Ctrl_Combined_.*_top.*_Cor.txt"), function(e){gsub("ann4Ath_Ath_Ctrl_Combined_", "", e)}), function(e){gsub("_Cor.txt", "", e)})) %>% rowwise() %>% mutate(matchid = paste(strsplit(matchid, "_")[[1]][2], "_", ifelse(grepl("ABA", matchid, fixed = TRUE), 1, 2), sep = "")) %>% column_to_rownames(var = "matchid")


## Consensus cell type using Brady et al.
Ath_Ctrls_Cor_CellType$cBrady <- apply(Ath_Ctrls_Cor_CellType %>% select(matches("wBradyHVG_top.*_cellident")), 1, function(x) {names(sort(table(as.character(x)), decreasing = TRUE)[1])})
## Consensus cell type using Li et al.
Ath_Ctrls_Cor_CellType$cLiwStele <- apply(Ath_Ctrls_Cor_CellType %>% select(matches("wLiHVG_top.*_cellident")) , 1, function(x) {names(sort(table(as.character(x)), decreasing = TRUE)[1])})
Ath_Ctrls_Cor_CellType$cLi <- apply(Ath_Ctrls_Cor_CellType %>% select(matches("wLiHVG_top.*_cellident")) , 1, function(x) {ifelse(names(sort(table(as.character(x)), decreasing = TRUE)[1]) == "Stele", names(sort(table(as.character(x)), decreasing = TRUE)[2]), names(sort(table(as.character(x)), decreasing = TRUE)[1]))})

Ath_Ctrls_Cor_cCellType <- Ath_Ctrls_Cor_CellType %>% rownames_to_column(var = "cell") %>% rowwise() %>% mutate(cBradyCellType = as.character(BradyConversion[cBrady]), cLiwSteleCellType = as.character(LiwSteleConversion[cLiwStele]), cLiCellType = as.character(LiConversion[cLi])) %>% select(cell, cBradyCellType, cLiwSteleCellType, cLiCellType) %>% column_to_rownames(var = "cell")

#######################################################################################################################################
#######################################################################################################################################




















#######################################################################################################################################
## 2. marker based method
#######################################################################################################################################
## running SEMITONES with different parameters, see ann4Ath_SEMITONES.py


#***************************************************************************************************************************************************#
#***************************************************************************************************************************************************#
## get concensus annotation from annotation transfer
Ath_Ctrls_SEMITONES_CellType <- read.table("Ath_Ctrls_RNAscaledata_CellTypes_noStele.txt", sep = "\t", header = TRUE, row.names = 1)

#- Considering adding weights to annotations at each nsds. The larger the nsds is, the heavier the wight will be.
# Ath_Ctrls_SEMITONES_CellType$cSEMITONES <- apply(Ath_Ctrls_SEMITONES_CellType, 1, function(x) {names(sort(table(as.character(x)), decreasing = TRUE)[1])})
Ath_Ctrls_SEMITONES_CellType$cSEMITONES <- apply(Ath_Ctrls_SEMITONES_CellType, 1, function(x) {
    anns <- sort(table(as.character(x)), decreasing = TRUE)
    if (length(anns) > 1) {
        first <- anns[1]
        second <- anns[2]
        if (as.numeric(first) == as.numeric(second)) {
            ## adding weights to annotations at each nsds
            x_rev <- rep(as.character(x), times = c(rep(c(0, 1, 2, 3, 4), each = 2), 5)) %>% .[. %in% c(names(first), names(second))]
            names(sort(table(x_rev), decreasing = TRUE)[1])
        } else {
            names(first)
        }
    } else {
        names(anns[1])
    }
})

Ath_Ctrls_SEMITONES_cCellType <- Ath_Ctrls_SEMITONES_CellType %>% rownames_to_column(var = "cell") %>% rowwise() %>% mutate(cSEMITONESCellType = as.character(SEMITONESConversion[cSEMITONES])) %>% select(cell, cSEMITONESCellType) %>% column_to_rownames(var = "cell")
#######################################################################################################################################
#######################################################################################################################################



















#######################################################################################################################################
## 3. annotation transfer
#######################################################################################################################################
### annotation transfer with different parameters
Shahan_Ath_RootAtlas <- readRDS("Shahan_Ath_RootAtlas.rds") %>% RunPCA(Shahan_Ath_RootAtlas, npcs = 500, verbose = FALSE, approx = FALSE) %>% RunUMAP(Shahan_Ath_RootAtlas, reduction = "pca", dims = 1:500, umap.method = "umap-learn", metric = "correlation", return.model = TRUE)
Ath_Ctrls <- readRDS("integratedAth") %>% subset(., subset = intspace_treatment == "Ctrl")


for (ndim in c(30, seq(50, 500, 50))){
	for (kweight in c(5, 10, 20, 30, 50, 80, 100)){
	tryCatch({
		Ath_TransferAnchors <- FindTransferAnchors(reference = Shahan_Ath_RootAtlas, query = Ath_Ctrls, normalization.method = "LogNormalize", reference.reduction = "pca", npcs = ndim, dims = 1:ndim)
		Ath_TransferredLabels <- TransferData(anchorset = Ath_TransferAnchors, refdata = Shahan_Ath_RootAtlas$celltype.anno.crude, k.weight = kweight, dims = 1:ndim)
		saveRDS(Ath_TransferredLabels, file = paste0("Ath_PredictedID_ndim", ndim, "kweight", kweight, ".RDS"))
		
		}, error = function(e){print(sprintf("ERROR at ndim %s and kweight %s: %s", ndim, kweight, conditionMessage(e)))}
		)
	}
}

#***************************************************************************************************************************************************#
#***************************************************************************************************************************************************#
## get concensus annotation from annotation transfer

AnnotationList <- list()
for (ndim in c(30, seq(50, 500, 50))){
    for (kweight in c(5, 10, 20, 30, 50, 80, 100)){
        if (file.exists(paste0("Ath_PredictedID_ndim", ndim, "kweight", kweight, ".RDS"))){
            AnnotationList[[paste0("ndim", ndim, "kweight", kweight)]] <- readRDS(paste0("Ath_PredictedID_ndim", ndim, "kweight", kweight, ".RDS")) %>% select(predicted.id, prediction.score.max) %>% rename_with( ~ paste0("ndim", ndim, "kweight", kweight, "_", .x)) %>% rownames_to_column(var = "cellid")
        }
    }    
}

Annotations <- Reduce(function(x,y) {dplyr::left_join(x, y, by = "cellid")}, AnnotationList) %>% column_to_rownames(var = "cellid")
ConsensusAnnotations <- Annotations %>% select(contains("_predicted.id")) %>% mutate(consensusann = apply(., 1, function(ann) {names(which.max(table(ann)))}), annfrequency = apply(., 1, function(ann) {sort(table(ann), decreasing = TRUE)[1]})) %>% select(contains("ann")) %>% cbind(Annotations, .)
ConsensusAnnotations$average.score <- sapply(1:nrow(ConsensusAnnotations), function(i) {mean(as.numeric(ConsensusAnnotations[i, ][which(ConsensusAnnotations[i, 1:98] == ConsensusAnnotations[i, 99]) + 1]))})
Ath_Ctrls_ConsensusPredictedCellType <- ConsensusAnnotations %>% rownames_to_column(var = "cell") %>% rowwise() %>% mutate(PredictedCellType = as.character(PredictedIDConversion[consensusann])) %>% select(cell, PredictedCellType) %>% column_to_rownames(var = "cell")
#######################################################################################################################################
#######################################################################################################################################












#######################################################################################################################################
## 4. find the consensus annotation
#######################################################################################################################################
Ath_cAnnotations <- Reduce(function(x, y, ...) {dplyr::full_join(x %>% rownames_to_column(var = "geneid"), y %>% rownames_to_column(var = "geneid"), by = "geneid", all = FALSE, ...) %>% column_to_rownames(var = "geneid")}, list(Ath_Ctrls_Cor_cCellType, Ath_Ctrls_SEMITONES_cCellType, Ath_Ctrls_ConsensusPredictedCellType))


Ath_cAnnotations$cCellType <- apply(Ath_cAnnotations %>% select(cBradyCellType, cLiCellType, cSEMITONESCellType, PredictedCellType), 1, function(x) {
    anns <- sort(table(as.character(x)), decreasing = TRUE)
    if (length(anns) > 1) {
        first <- anns[1]
        second <- anns[2]
        ## when the there are equivalent supports for the top two annotations, use the annotation based on cell type markers.
        if (as.numeric(first) == as.numeric(second)) {
            x["cSEMITONESCellType"]
        } else {
            names(first)
        }
    } else {
        names(anns[1])
    }
})
#######################################################################################################################################
#######################################################################################################################################


print("Done")