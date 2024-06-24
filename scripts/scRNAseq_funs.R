## This file will be used to store common functions that are useful in scRNA-seq anlysis.
## Quote, Unquote and Quasiquote: https://masterr.org/r/quote-n-quasiquote/
## Quasiquotation: https://adv-r.hadley.nz/quasiquotation.html
## Non-standard evaluation: http://adv-r.had.co.nz/Computing-on-the-language.html
## https://sebastiansauer.github.io/prop_fav/#:~:text=See%20the%20!!%20%3B%20they%20mean,out%20of%20a%20data%20frame.
## as.name vs sym: https://stackoverflow.com/questions/55870158/whats-the-difference-between-as-name-and-sym
## https://stackoverflow.com/questions/61380829/what-is-the-difference-among-evalparse-evalstr2lang-evalstr2expressio
## Converting a String to a Variable Name On-The-Fly and Vice-versa in R: https://www.r-bloggers.com/2010/12/converting-a-string-to-a-variable-name-on-the-fly-and-vice-versa-in-r/




#####################################################################################################################################################################################
## util
read_all_sheets <- function(xlsxFile, ...) {
	sheet_list <- list()
	sheet_names <- openxlsx::getSheetNames(xlsxFile)
	for (sn in sheet_names) {
		sheet_list[[sn]] <- openxlsx::read.xlsx(xlsxFile, sheet = sn, ...)
	}
	return(sheet_list)
}



#####################################################################################################################################################################################

CellIndAgg <- function(MetaData, res, identcol, topcutoff = 0.9){
## CellIndAgg is a function to aggregate cell identity information from a meta.data when each cell is assigned to a cell identity (e.g. by matching to reference transcriptomes)
## MetaData, the meta.data slot from seurat object with a column for cell identity 
## MetaData basic format: orig.ident nCount_RNA nFeature_RNA percent.mt integrated_snn_res.1.5 seurat_clusters cellident
## res, resolution for FindClusters
## identcol, the name of the column for cell identity
## topcutoff, cutoff used for finding top cell identities, cumulative percentage of top cell identities should be larger than this cutoff 
    print(sprintf("Currently running cell cluster assignment at resolution of %s...", res))
    rescol <- paste("integrated_snn_res.", res, sep = "")
    res_temp <- MetaData[, c(rescol, identcol)]

    clstrList <- list()
    for (cluster in as.numeric(levels(res_temp[[rescol]]))){
        res_cluster_temp <- subset(res_temp, get(rescol) == cluster)
        if(nrow(res_cluster_temp) > 0) {
            celltype_count_temp <- setNames(aggregate(data = res_cluster_temp, get(rescol) ~ get(identcol), function(x) {length(x)}), c(identcol, "integrated_snn_res")) %>% arrange(desc(integrated_snn_res)) %>% mutate(cumperc = cumsum(integrated_snn_res)/sum(integrated_snn_res))
            #cat(res, "\t", cluster, "\t", sum(celltype_count_temp$cumperc <= topcutoff)+1, "\t", celltype_count_temp[[identcol]][1], "\t", round(celltype_count_temp$cumperc[1], digits = 2), "\t", paste0(celltype_count_temp[0:sum(celltype_count_temp$cumperc <= topcutoff)+1, ][[identcol]], collapse = ","), "\t", paste0(sapply(celltype_count_temp[0:sum(celltype_count_temp$cumperc <= topcutoff)+1, ]$cumperc, function(x) round(x, digits = 2)), collapse = ","), "\n", sep = "")
            clstrList[cluster + 1] <- paste(res, "\t", cluster, "\t", nrow(res_cluster_temp), "\t", celltype_count_temp[[identcol]][1], "\t", round(celltype_count_temp$cumperc[1], digits = 2), "\t", sum(celltype_count_temp$cumperc <= topcutoff)+1, "\t", paste0(celltype_count_temp[0:sum(celltype_count_temp$cumperc <= topcutoff)+1, ][[identcol]], collapse = ","), "\t", paste0(sapply(celltype_count_temp[0:sum(celltype_count_temp$cumperc <= topcutoff)+1, ]$cumperc, function(x) round(x, digits = 2)), collapse = ","), sep = "")
        } else {
            clstrList[cluster + 1] <- paste(res, "\t", cluster, "\t", 0, "\t", NA, "\t", NA, "\t", NA, "\t", NA, "\t", NA, sep = "")
        }
    }
    clstrInfo <- as.data.frame(do.call(rbind, sapply(clstrList, function(x) {strsplit(x, "\t")}))) %>% setNames(., c("res","cluster","numcells","topident","topidentperc","numtopidents","topidents","topidentcumpercs"))
    return(clstrInfo)
}


# AssignCellType <- function(SeuObject, BradyAnnFile, LiAnnFile, res = NULL){
# AssignCellType assigns cell type to Seurat clusters based on correlation to existing reference transcriptomes
# SeuObject, a seurat object, after running FindNeighbors, FindClusters
# res, which resolution to use for cell clusters, default "seurat_clusters"
# BradyAnnFile, annotation file when comparing to Brady et al. 2007 dataset
# LiAnnFile, annotation file when comparing to Li et al. 2016 dataset 
    # ifelse(is.null(res), Clusters <- "seurat_clusters", Clusters <- paste0("integrated_snn_res.", res))
    # SeuObject_Ctrl_Metadata <- subset(SeuObject@meta.data, subset = orig.ident %in% c("Ath_ABA_Ctrl_R1", "Ath_NaCl_Ctrl_R1")) %>% rownames_to_column(var = "matchid")

    # BradyAnn <- read.delim(BradyAnnFile, header = TRUE) %>% rownames_to_column(var = "rowid") %>% rename_with( ~ paste0("wBrady" , .x), .cols = everything()) %>% rowwise() %>% mutate(matchid = paste(strsplit(wBradyrowid, "_")[[1]][2], "_", ifelse(strsplit(wBradyrowid, "_")[[1]][1] == "ABACtrl", 1, 2), sep = ""), wBradycellidentm = strsplit(wBradycellident, "_")[[1]][1]) %>% mutate(wBradycellidentm = unlist(BradyCellType[[wBradycellidentm]]["Tissue"], use.names = FALSE))
    # LiAnn <- read.delim(LiAnnFile, header = TRUE) %>% rownames_to_column(var = "rowid") %>% rename_with( ~ paste0("wLi" , .x), .cols = everything()) %>% rowwise() %>% mutate(matchid = paste(strsplit(wLirowid, "_")[[1]][2], "_", ifelse(strsplit(wLirowid, "_")[[1]][1] == "ABACtrl", 1, 2), sep = "")) %>% mutate(wLicellident = unlist(LiCellType[[wLicellident]]["Tissue"], use.names = FALSE))

    # SeuObject_Ctrl_CellType <- Reduce(function(x, y, ...) {merge(x, y, by = "matchid", all = FALSE, ...)}, list(SeuObject_Ctrl_Metadata, BradyAnn, LiAnn))

    # SeuObject_Ctrl_wBradyCellTypeSec <- CellIndAgg(SeuObject_Ctrl_CellType, res, "wBradycellident", topcutoff = 0.8)
    # SeuObject_Ctrl_wBradyCellType <- CellIndAgg(SeuObject_Ctrl_CellType, res, "wBradycellidentm", topcutoff = 0.8)
    # SeuObject_Ctrl_wLiCellType <- CellIndAgg(SeuObject_Ctrl_CellType, res, "wLicellident", topcutoff = 0.8)

    # SeuObject@meta.data$wBradycellidentSecAgg <- sapply(SeuObject@meta.data[[Clusters]], function(x) {SeuObject_Ctrl_wBradyCellTypeSec[as.numeric(x), c("topidents")]})
    # SeuObject@meta.data$wBradycellidentAgg <- sapply(SeuObject@meta.data[[Clusters]], function(x) {SeuObject_Ctrl_wBradyCellType[as.numeric(x), c("topidents")]})
    # SeuObject@meta.data$wLicellidentAgg <- sapply(SeuObject@meta.data[[Clusters]], function(x) {SeuObject_Ctrl_wLiCellType[as.numeric(x), c("topidents")]})
    
    # return(SeuObject)
# }


DominantCellType <- function(SeuObject, AnnFile, AnnCol, res = NULL, cutoff = 0.8){
## DominantCellType assigns cell type to Seurat clusters based on correlation of each cell to existing reference transcriptomes
## SeuObject, a seurat object, after running FindNeighbors, FindClusters
## res, which resolution to use for cell clusters, default "seurat_clusters"
## AnnFile, annotation file for each cell when comparing to reference transcriptomes, rownames should be gene IDs, colnames should be cell IDs which must match those in the SeuObject
## AnnCol, the column containing cell type annotation in AnnFile
## cutoff, the ratio of cells of the assigned cell types in a given cell cluster should be larger than this cutoff, default: 0.8
    ifelse(is.null(res), Clusters <- "seurat_clusters", Clusters <- paste0("integrated_snn_res.", res))

    SeuObject_Metadata <- SeuObject@meta.data %>% rownames_to_column(var = "matchid")
    CorAnn <- AnnFile %>% rownames_to_column(var = "matchid")
    SeuObject_CellType <- Reduce(function(x, y, ...) {merge(x, y, by = "matchid", all = FALSE, ...)}, list(SeuObject_Metadata, CorAnn))
    SeuObject_CellTypeSec <- CellIndAgg(SeuObject_CellType, res, AnnCol, topcutoff = cutoff)
    SeuObject@meta.data$cellident <- sapply(SeuObject@meta.data[[Clusters]], function(x) {SeuObject_CellTypeSec[as.numeric(x), c("topidents")]})
#	SeuObject_Metadata$cellident <- sapply(SeuObject_Metadata[[Clusters]], function(x) {SeuObject_CellTypeSec[as.numeric(x), c("topidents")]})
    return(list(SeuObject_CellTypeSec, SeuObject))
}


## Manually convert to CellDataSet for monocle3
ConvertToCDS <- function(SeuratObject, assay = "RNA", slot = "counts") {
## ConvertToCDS convert a Seurat object to a basic CellDataSet object
## assay, which assay in SeuratObject to use, "RNA" or "integrated"
## slot, which slot in the assay to use, "counts" or "data"
    expression_temp <- as(as.matrix(GetAssayData(SeuratObject, assay = assay, slot = slot)), "sparseMatrix")
    celldata_temp <- data.frame(SeuratObject@meta.data) #%>% select(orig.ident, nCount_RNA, nFeature_RNA)
    featuredata_temp <- data.frame(gene_short_name = row.names(expression_temp), row.names = row.names(expression_temp))
    rawCDS_temp <- new_cell_data_set(expression_data = expression_temp, cell_metadata = celldata_temp, gene_metadata = featuredata_temp)
    return(rawCDS_temp)
}


##calculate the coefficient of variation (relative standard deviation, which is the ratio of the standard deviation {σ} to the mean {μ} (or its absolute value, {|μ|})) for each gene in the reference datasets.
CVcal <- function(dataset, idcol = 0, cols, valcutoff = NULL){
# CVcal calculates the coefficient of variation for each gene in the given datasets.
# idcol, column number of unique IDs, 0 will be the rownames of the dataset.
# cols, column containing values for CV calulation, a gene will be removed from the calculation if there are NAs within these colums.
# valcutoff, cutoff for values to be used for CV calulation.
    CVs <- c()
	filtereddataset <- dataset[complete.cases(dataset[, cols]),]
	print(sprintf("For %s, the coefficient of variation for %s features will be calculated among %s samples.", deparse(substitute(dataset)), nrow(filtereddataset), ncol(filtereddataset)))
    for (nfeatures in 1:nrow(filtereddataset)){
        if (is.null(valcutoff)) {
            ifelse(mean(as.numeric(filtereddataset[nfeatures, cols])) == 0, CVs <- append(CVs, NA), CVs <- append(CVs, sd(as.numeric(filtereddataset[nfeatures, cols]))/mean(as.numeric(filtereddataset[nfeatures, cols]))))
        } else {
            ifelse(mean(as.numeric(filtereddataset[nfeatures, cols])) >= valcutoff, CVs <- append(CVs, sd(as.numeric(filtereddataset[nfeatures, cols]))/mean(as.numeric(filtereddataset[nfeatures, cols]))), CVs <- append(CVs, NA))
        }
    }
    ifelse(idcol == 0 , names(CVs) <- rownames(filtereddataset), names(CVs) <- c(sapply(filtereddataset[, idcol], as.character)))
	print(sprintf("The coefficients of variation for a total of %s were calculated for %s.", length(CVs), deparse(substitute(dataset))))
    return(CVs)
}


ClusterProfiler <- function(seuobject, data, res, norm = FALSE, scale.factor = 10000){
## ClusterProfiler generate cell cluster transcriptomes from the input Seurat object.
## seuobject, a seurat object, after running FindNeighbors, FindClusters
## data, "counts" for raw counts or "data" for normalized counts in @assays$RNA
## res, which resolution to use for cell clusters, default "seurat_clusters", should be present in @meta.data in the seuobject
## res, which resolution to use for cell clusters, default "seurat_clusters", should be present in @meta.data in the seuobject
## norm, whether to normalized the expression counts for each cell, default FALSE. If true, the expression counts for each cell will be divided by the total counts for that cell and multiplied by the scale.factor.	
    counts <- slot(seuobject@assays$RNA, data)
	if (norm) {
		print(sprintf("Feature counts in each cell will be normalized by the total counts for each cell and multiplied by %s.", scale.factor))
		counts <- apply(counts, 2, function(featurecounts){featurecounts/sum(featurecounts)}) * scale.factor
	}
	
    clstrres <- paste0("integrated_snn_res.", res)
    clstrprofile <- matrix(nrow = nrow(counts), ncol = 0)
	clstrids <- sort(as.numeric(unique(seuobject@meta.data[[clstrres]])))-1
    
	for (cluster in clstrids) {
        seumetadata <- seuobject@meta.data[seuobject@meta.data[[clstrres]] == cluster, ]
        clstrcells <- as.matrix(counts[ ,colnames(counts) %in% rownames(seumetadata)])
        #cat("Currently processing ", cluster, "\t", nrow(seumetadata), "\t", ncol(clstrcells), "\n", sep = "")
		if (ncol(clstrcells) > 0){
#			print(sprintf("Currently processing C%s in %s, which has %s cells...", cluster, deparse(substitute(seuobject)), ncol(clstrcells)))
			print(sprintf("Currently processing C%s, which has %s cells...", cluster, ncol(clstrcells)))
			clstrprofile <- cbind(clstrprofile, rowSums(clstrcells)/ncol(clstrcells))
		} else {
#			print(sprintf("C%s in %s does not have any cells, please check your input. The expression profile for this cluster will be set to all 0s.", cluster, deparse(substitute(seuobject))))
			print(sprintf("C%s does not have any cells, please check your input. The expression profile for this cluster will be set to all 0s.", cluster))
			clstrprofile <- cbind(clstrprofile, rep(0, times = nrow(counts)))
		}
    }
    colnames(clstrprofile) <- paste0("C", clstrids)
    return(clstrprofile)
}


ClstrProfileIdent <- function(clstrprofile, ref, startcol, idcol, features = NULL){
##ClstrProfileIdent assign cell types to cell clusters based on correlation (spearman) of cell cluster transcriptome to existing reference transcriptomes.
##clstrprofile, a dataframe where each rownames are gene ids, column reprensents the aggregated transcriptome profile of a cell cluster.
##ref, the reference transcriptomes, each column represents the transcriptome profile of a cell type (this input could have extra columns for feature ids before the transcriptomes).
##startcol, the column number where the transcriptome profiles for cell types start in the reference set.
##idcol, id column in ref, e.g. "AtID", column "AtID" in the ref will be used to combine the clstrprofile and ref.
##features, only use the features in this list to calculate the correlation, should be a vector, e.g. c("gene1", "gene2", ...)
    ClusterProfile_wRef <- merge(ref, clstrprofile, by.x = idcol, by.y = 0)
#    ifelse(is.null(features), ClusterProfile_wRef <- ClusterProfile_wRef, ClusterProfile_wRef <- ClusterProfile_wRef[ClusterProfile_wRef[[idcol]] %in% features, ])
    if(!is.null(features)){ClusterProfile_wRef <- ClusterProfile_wRef[ClusterProfile_wRef[[idcol]] %in% features, ]}
    cat(sprintf("%s features were used for correlation calculation.", nrow(ClusterProfile_wRef)), "\n", sep = "")    
  
    mergedtable_cell <- suppressWarnings(sapply((ncol(ref)+1):ncol(ClusterProfile_wRef), function(i) sapply(startcol:ncol(ref), function(j) cor.test(ClusterProfile_wRef[, i], ClusterProfile_wRef[, j], method = "spearman")[c("p.value", "estimate")]))) ##Each column represents a cell; each row represents a reference cell type.
    mergedtable_cell_cor <- mergedtable_cell[seq(2,nrow(mergedtable_cell),2),]
    mergedtable_cell_maxcor <- sapply(1:ncol(mergedtable_cell_cor), function(i) max(as.numeric(mergedtable_cell_cor[,i])))
    mergedtable_cell_pvalue <- mergedtable_cell[seq(1,nrow(mergedtable_cell)-1,2),]
    mergedtable_cell_cellident <- sapply(1:ncol(mergedtable_cell_cor), function(i) colnames(ClusterProfile_wRef)[startcol:ncol(ref)][which(as.numeric(mergedtable_cell_cor[,i]) == max(as.numeric(mergedtable_cell_cor[,i])))])
    mergedtable_cell_maxp <- sapply(1:ncol(clstrprofile), function(i) as.numeric(mergedtable_cell_pvalue[,i])[which(as.numeric(mergedtable_cell_cor[,i]) == max(as.numeric(mergedtable_cell_cor[,i])))])

    mergedtable_cellcortable <- data.frame(clusters = colnames(ClusterProfile_wRef)[(ncol(ref)+1):ncol(ClusterProfile_wRef)], maxcor = as.character(mergedtable_cell_maxcor), maxp = as.character(mergedtable_cell_maxp), cellident = as.character(mergedtable_cell_cellident))

    return(mergedtable_cellcortable)
}


MultiClusterProfiler <- function(seuobject, groupcol, data, res, norm = FALSE, scale.factor = 10000){
## MultiClusterProfiler generates cell cluster transcriptomes for each group in the groupcol in the input Seurat object, and merge them into one table.
## seuobject, a seurat object, after running FindNeighbors, FindClusters
## groupcol, which coloum in @meta.data contains information for cell grouping
## data, "counts" for raw counts or "data" for normalized counts in @assays$RNA
## res, which resolution to use for cell clusters, default "seurat_clusters", should be present in @meta.data in the seuobject
## norm, whether to normalized the expression counts for each cell, default FALSE. If true, the expression counts for each cell will be divided by the total counts for that cell and multiplied by the scale.factor.	
	clstrres <- paste0("integrated_snn_res.", res)
	clstrids <- sort(as.numeric(unique(seuobject@meta.data[[clstrres]])))-1 ##length(unique(seuobject@meta.data[[clstrres]]))
	
	ClusterProfiler <- function(seuobject, data, res, norm = FALSE, scale.factor = 10000){
		counts <- slot(seuobject@assays$RNA, data)
		if (norm) {
			print(sprintf("Feature counts in each cell will be normalized by the total counts for each cell and multiplied by %s.", scale.factor))
			counts <- apply(counts, 2, function(featurecounts){featurecounts/sum(featurecounts)}) * scale.factor
		}
		clstrprofile <- matrix(nrow = nrow(counts), ncol = 0)
		for (cluster in clstrids) {
			seumetadata <- seuobject@meta.data[seuobject@meta.data[[clstrres]] == cluster, ]
			clstrcells <- as.matrix(counts[ ,colnames(counts) %in% rownames(seumetadata)])
			if (ncol(clstrcells) > 0){
				print(sprintf("Currently processing C%s, which has %s cells...", cluster, ncol(clstrcells)))
				clstrprofile <- cbind(clstrprofile, rowSums(clstrcells)/ncol(clstrcells))				
			} else {
				print(sprintf("C%s does not have any cells, please check your input. The expression profile for this cluster will be set to all 0s.", cluster))
				clstrprofile <- cbind(clstrprofile, rep(0, times = nrow(counts)))
			}
		}
		colnames(clstrprofile) <- paste0("C", clstrids)
		return(clstrprofile)
	}
	
	groups <- Seurat::FetchData(object = seuobject, vars = groupcol)
	
	subclstrprofiles <- list()
	for (group in unlist(unique(groups))){
		subseuobject <- seuobject[, which(x = groups == group)]
		print(sprintf("Currently processing %s in '%s' (a total of %s cells)...", group, groupcol, ncol(subseuobject)))
		subclstrprofiles[[group]] <- ClusterProfiler(subseuobject, data = data, res = res, norm = norm, scale.factor = scale.factor) %>% as.data.frame() %>% rename_with(function(x){paste0(group, "_" , x)})
	}
	multiprofiles <- Reduce(function(x, y, ...) {dplyr::full_join(x %>% rownames_to_column(var = "featureid"), y %>% rownames_to_column(var = "featureid"), by = "featureid", all = FALSE, ...) %>% column_to_rownames(var = "featureid")}, subclstrprofiles)
	return(multiprofiles)
}


## "ExpressionProfiler" can replace both "ClusterProfiler" and "MultiClusterProfiler"
ExpressionProfiler <- function(seuobject, assays = NULL, slot = "counts", group.by = "seurat_clusters", split.by = NULL, norm = FALSE, scale.factor = 10000, logspace = FALSE){
## ExpressionProfiler generates averaged expression profiles for each identity class in group.by for each category in split.by in the input Seurat object, and merge them into one table.
## seuobject, a seurat object.
## assays, which assays to use, default is all assays.
## slot, slot(s) to use, "counts" for raw counts or "data" for normalized counts in @assays$RNA. If slot is set to "data", this function assumes that the data has been log normalized (natural-log transformed using log1p) and therefore feature values are exponentiated prior to averaging so that averaging is done in non-log space, the averages were then again natural-log transformed using log1p. Otherwise, if slot is set to either "counts" or "scale.data", no exponentiation is performed prior to averaging.
## group.by, name of the coloum in @meta.data containing information for identity class, default is "seurat_clusters"
## split.by, name of the coloum in @meta.data containing information of categories to split the data
## norm, whether to normalize the data by dividing feature counts for each cell by the total counts for that cell and multiplied by the scale.factor
## logspace, whether to transformed normalize counts using natural-log with log1p
#If the data has been "LogNormalize"d, "counts" will contain raw counts and "data" will contain normalized counts in @assays$RNA; if the data has been "sctransform"ed, "counts" and "data" will be the same, both contain raw counts, the normalized values from "sctransform" will be in @assays$SCT@data. 

	if (is.null(group.by)) {
		identities <- "pseudobulk"
	} else if (is.null(levels(seuobject@meta.data[[group.by]]))){
		identities <- sort(unique(seuobject@meta.data[[group.by]]))
	} else {
		identities <- levels(seuobject@meta.data[[group.by]])
	}
#ifelse returns a value with the same shape as test which is filled with elements selected from either yes or no depending on whether the element of test is TRUE or FALSE.
#ifelse function in R only returns the first element of the vector if there is only one test value. 
#> ifelse(1>0, c(2, 3), c(4, 5)) 
#[1] 2 

	IdentProfiler <- function(seuobject, slot = "counts", group.by = "seurat_clusters", norm = FALSE, scale.factor = 10000, logspace = FALSE){
		print(sprintf("Slot '%s' from @assays$RNA will be used.", slot))
		if(slot == "data") {
			counts <- expm1(slot(seuobject@assays$RNA, slot))
		} else if (slot == "counts") { ##slot %in% c("counts", "scale.data")
			counts <- slot(seuobject@assays$RNA, slot)
			if (norm) {
				print(sprintf("Raw counts in each cell will be normalized by the total counts for each cell and multiplied by %s.", scale.factor))
				counts <- apply(counts, 2, function(featurecounts){featurecounts/sum(featurecounts)}) * scale.factor
			}
		} else {
			stop(sprintf("'%s' is not valid. Please use valid slot, 'data' or 'counts'!", slot))
		}
    
		identprofile <- matrix(nrow = nrow(counts), ncol = 0)
		if (is.null(group.by)) {
			identprofile <- cbind(identprofile, rowSums(counts)/ncol(counts))
		} else {
			for (identity in identities) {
				identcells <- seuobject@meta.data[seuobject@meta.data[[group.by]] == identity, ]
				if (nrow(identcells) > 0){
					identexp <- as.matrix(counts[ ,colnames(counts) %in% rownames(identcells)])
					print(sprintf("Currently processing %s, which has %s cells...", identity, ncol(identexp)))
					# if (slot == "data") {
						# identprofile <- cbind(identprofile, log1p(rowSums(identexp)/ncol(identexp)))
					# } else if (slot == "counts") {
						# identprofile <- cbind(identprofile, rowSums(identexp)/ncol(identexp))
					# }
					identprofile <- cbind(identprofile, rowSums(identexp)/ncol(identexp))
				} else {
					print(sprintf("'%s' does not have any cells, please check your input. The expression profile for this identity will be set to all 0s.", identity)) ##deparse(substitute(seuobject))
					identprofile <- cbind(identprofile, rep(0, times = nrow(counts)))
				}
			}
		}
		colnames(identprofile) <- identities
		
		if (logspace) {
			identprofile <- log1p(identprofile)
		}
		return(identprofile)
	}

	if (is.null(split.by)){
		print(sprintf("Generating expression profiles for %s: %s...", ifelse(is.null(group.by), "entire dataset", group.by), paste0(identities, collapse = ",")))
		expressionprofiles <- IdentProfiler(seuobject, slot = slot, group.by = group.by, norm = norm, scale.factor = scale.factor, logspace = logspace)
	} else {
#		groups <- ifelse(is.null(levels(seuobject@meta.data[[split.by]])), sort(unique(seuobject@meta.data[[split.by]])), levels(seuobject@meta.data[[split.by]]))
		if(is.null(levels(seuobject@meta.data[[split.by]]))) {
			cats <- sort(unique(seuobject@meta.data[[split.by]]))
		} else {
			cats <- levels(seuobject@meta.data[[split.by]])
		}
	
		print(sprintf("Generating expression profiles for %s: %s for %s: %s ...", ifelse(is.null(group.by), "entire dataset", group.by), paste0(identities, collapse = ","), split.by, paste0(cats, collapse = ",")))
	
		subcatprofiles <- list()
		for (subcat in cats){
			subseuobject <- seuobject[, which(x = seuobject@meta.data[[split.by]] == subcat)]
			print(sprintf("Currently processing '%s' in '%s' (a total of %s cells)...", ifelse(is.null(group.by), "entire dataset", group.by), subcat, ncol(subseuobject)))
			subcatprofiles[[subcat]] <- IdentProfiler(subseuobject, slot = slot, group.by = group.by, norm = norm, scale.factor = scale.factor, logspace = logspace) %>% as.data.frame() %>% rename_with(function(x){paste0(subcat, "_" , x)})
		}
		expressionprofiles <- Reduce(function(x, y, ...) {dplyr::full_join(x %>% rownames_to_column(var = "featureid"), y %>% rownames_to_column(var = "featureid"), by = "featureid", all = FALSE, ...) %>% column_to_rownames(var = "featureid")}, subcatprofiles)
	}
	return(expressionprofiles)
}

if (FALSE) {
IntegratedCor <- function(exptable, groupids, cormethod = "spearman"){
## IntegratedCor calculates pairwise correlation for each feature in the input expression table.
## exptable, the transcriptomic profile table generated by MultiClusterProfiler including all groups in the groupcol.
## groupids, group id used for each group in a list, e.g. list("1" = "Ath", "2" = "Esa", "3" = "Sir", "4" = "Spa"), make sure the numbers are unique.
## cormethod, correlation method will be used in cor.test, one of "pearson", "kendall", or "spearman", default: "spearman".
	
	corlist <- list()
	numgroup <- length(groupids)
	for (groupi in 1:(numgroup-1)){
		for (groupj in (groupi+1):numgroup){
			idi <- exptable %>% dplyr::select(contains(groupids[[groupi]])) %>% colnames() %>% sapply(., function(x){strsplit(x, groupids[[groupi]])[[1]][2]})
			idj <- exptable %>% dplyr::select(contains(groupids[[groupj]])) %>% colnames() %>% sapply(., function(x){strsplit(x, groupids[[groupj]])[[1]][2]})
			if (all(idi == idj)) {
				print(sprintf("Calculating correlation for %s and %s...", groupids[[groupi]], groupids[[groupj]]))
				cortable <- suppressWarnings(sapply(1:nrow(exptable), function(rownum){
					geneexpi <- exptable %>% dplyr::select(contains(groupids[[groupi]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupi-1)*numclstr)])
					geneexpj <- exptable %>% dplyr::select(contains(groupids[[groupj]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupj-1)*numclstr)])
					cor.test(geneexpi, geneexpj, method = cormethod)[c("estimate", "p.value")] %>% unlist() #%>% as.numeric()
				})) 
				corlist[[paste0(groupids[[groupi]], "_vs_", groupids[[groupj]])]] <- cortable %>% t() %>% as.data.frame() %>% setNames(., paste0(groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , c("cor", "pvalue"))) #rename_with(function(x){paste0(groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , x)})
			} else {
				warning(sprintf("Calculating correlation for %s and %s, but the profiles from them are not matching!", groupids[[groupi]], groupids[[groupj]]), noBreaks. = TRUE)
				## will try to contiune even when the ids are not matching.
				cortable <- suppressWarnings(sapply(1:nrow(exptable), function(rownum){
					geneexpi <- exptable %>% dplyr::select(contains(groupids[[groupi]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupi-1)*numclstr)])
					geneexpj <- exptable %>% dplyr::select(contains(groupids[[groupj]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupj-1)*numclstr)])
					cor.test(geneexpi, geneexpj, method = cormethod)[c("estimate", "p.value")] %>% unlist() #%>% as.numeric()
				})) 
				corlist[[paste0("mis_", groupids[[groupi]], "_vs_", groupids[[groupj]])]] <- cortable %>% t() %>% as.data.frame() %>% setNames(., paste0("mis_", groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , c("cor", "pvalue"))) #rename_with(function(x){paste0("mis_", groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , x)})
				
			}
		}
	}
	mergedcortable <- Reduce(function(x, y, ...) {dplyr::full_join(x %>% rownames_to_column(var = "matchid"), y %>% rownames_to_column(var = "matchid"), by = "matchid", ...) %>% column_to_rownames(var = "matchid")}, corlist)
	rownames(mergedcortable) <- rownames(exptable)
	return(mergedcortable)
}
}

IntegratedCor <- function(exptable, groupids, cormethod = "spearman"){
## IntegratedCor calculates pairwise correlation for each feature in the input expression table.
## exptable, the transcriptomic profile table generated by MultiClusterProfiler including all groups in the groupcol.
## groupids, group id used for each group in a list, e.g. list("1" = "Ath", "2" = "Esa", "3" = "Sir", "4" = "Spa"), make sure the numbers are unique.
## cormethod, correlation method will be used in cor.test, one of "pearson", "kendall", or "spearman", default: "spearman".
	
	corlist <- list()
	numgroup <- length(groupids)
	for (groupi in 1:(numgroup-1)){
		for (groupj in (groupi+1):numgroup){
			idi <- exptable %>% dplyr::select(contains(groupids[[groupi]])) %>% colnames() %>% sapply(., function(x){strsplit(x, groupids[[groupi]])[[1]][2]})
			idj <- exptable %>% dplyr::select(contains(groupids[[groupj]])) %>% colnames() %>% sapply(., function(x){strsplit(x, groupids[[groupj]])[[1]][2]})
			if (all(idi == idj)) {
				print(sprintf("Calculating correlation for %s and %s...", groupids[[groupi]], groupids[[groupj]]))
				cortable <- suppressWarnings(sapply(1:nrow(exptable), function(rownum){
					geneexpi <- exptable %>% dplyr::select(contains(groupids[[groupi]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupi-1)*numclstr)])
					geneexpj <- exptable %>% dplyr::select(contains(groupids[[groupj]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupj-1)*numclstr)])
					if (!any(is.na(geneexpi)) & !any(is.na(geneexpj))) {
						cor.test(geneexpi, geneexpj, method = cormethod)[c("estimate", "p.value")] %>% unlist() #%>% as.numeric()
					} else {
						return(c(NA, NA)) ## return is not needed here, but is kept to make the code clean.
					}
				})) 
				corlist[[paste0(groupids[[groupi]], "_vs_", groupids[[groupj]])]] <- cortable %>% t() %>% as.data.frame() %>% setNames(., paste0(groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , c("cor", "pvalue"))) #rename_with(function(x){paste0(groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , x)})
			} else {
				warning(sprintf("Calculating correlation for %s and %s, but the profiles from them are not matching!", groupids[[groupi]], groupids[[groupj]]), noBreaks. = TRUE)
				## will try to contiune even when the ids are not matching.
				cortable <- suppressWarnings(sapply(1:nrow(exptable), function(rownum){
					geneexpi <- exptable %>% dplyr::select(contains(groupids[[groupi]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupi-1)*numclstr)])
					geneexpj <- exptable %>% dplyr::select(contains(groupids[[groupj]])) %>% slice(rownum) %>% as.numeric() #as.numeric(exptable[rownum, c(1:numclstr + (groupj-1)*numclstr)])
					if (!any(is.na(geneexpi)) & !any(is.na(geneexpj))) {
						cor.test(geneexpi, geneexpj, method = cormethod)[c("estimate", "p.value")] %>% unlist() #%>% as.numeric()
					} else {
						return(c(NA, NA)) ## return is not needed here, but is kept to make the code clean.
					}
				})) 
				corlist[[paste0("mis_", groupids[[groupi]], "_vs_", groupids[[groupj]])]] <- cortable %>% t() %>% as.data.frame() %>% setNames(., paste0("mis_", groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , c("cor", "pvalue"))) #rename_with(function(x){paste0("mis_", groupids[[groupi]], "_vs_", groupids[[groupj]], "_" , x)})
				
			}
		}
	}
	mergedcortable <- Reduce(function(x, y, ...) {dplyr::full_join(x %>% rownames_to_column(var = "matchid"), y %>% rownames_to_column(var = "matchid"), by = "matchid", ...) %>% column_to_rownames(var = "matchid")}, corlist)
	rownames(mergedcortable) <- rownames(exptable)
	return(mergedcortable)
}


CumSumPerc <- function(num.vector, perc = 0.5, interval = "top", use.abs = FALSE){
## CumSumPerc is a function to extract top (or bottom) values whose sum exceeds a cetain percentage of the total sum of the input numeric vector.
## num.vector, a vector of numbers, which will be used to extract top (or bottom) values whose sum exceeds a cetain percentage of the total sum. 
## perc, percentage cutoff of the cumulative sum, default: 0.5
## interval, whether to return the top or the bottom values whose sum meets the cutoff, default: top
## use.abs, whether or not to convert every element in num.vector to absolute value during the extraction, default: FALSE
	
	num.vector <- num.vector[!is.na(num.vector)]
	if (use.abs) {
		num.vector.used <- abs(num.vector)
	} else {
		num.vector.used <- num.vector
	}

	if (interval == "top") {
		larger.vector <- num.vector[order(num.vector.used)][cumsum(sort(num.vector.used)) >= sum(num.vector.used) * (1 - perc)]
		
		if (use.abs) {
			larger.vector.used <- abs(larger.vector)
		} else {
			larger.vector.used <- larger.vector
		}
		print(sprintf("The output contains %s top contributors that cumulatively contributes to %1.2f%% of the total of the input.", length(larger.vector), sum(larger.vector.used)/sum(num.vector.used) * 100))
		
		return(larger.vector)
	} else if (interval == "bottom") {
		#smaller.vector <- num.vector[order(num.vector.used)][cumsum(sort(num.vector.used)) <= sum(num.vector.used) * perc]
		smaller.vector <- num.vector[order(num.vector.used)][1:(sum(cumsum(sort(num.vector.used)) <= sum(num.vector.used) * perc) + 1)]
		
		if (use.abs) {
			smaller.vector.used <- abs(smaller.vector)
		} else {
			smaller.vector.used <- smaller.vector
		}
		print(sprintf("The output contains %s bottom contributors that cumulatively contributes to %1.2f%% of the total of the input.", length(smaller.vector), sum(smaller.vector.used)/sum(num.vector.used) * 100))
		
		return(smaller.vector)
	} else {
		print(sprintf("'%s' is not recognized. Please use 'top' or 'bottom' for interval.", interval))
	}
}


##Comparisons of various types of normality tests
###Shapiro-Wilk test of normality (stats::shapiro.test): the number of non-missing values must be between 3 and 5000.
###D'Agostino test (moments::agostino.test): sample size must be between 8 and 46340
FindSigOutliers <- function(nums, test = "lilliefors", pcut = 0.05) {
## FindSigOutliers is a function to identify elements in a vector that contribute to non-normality of the vector (normality test is done with nortest::lillie.test)
## nums, a named vector-like object
## test, normality test to be used, "anderson-darling" for Anderson-Darling Test, "lilliefors" for Lilliefors (Kolmogorov-Smirnov) Test
## pcut, p-value cutoff
 
	## check if it is a named vector.
	if (is.null(names(nums)) | !is.vector(nums)) {
		stop(sprintf("%s is not a valid input.", deparse(substitute(nums))))
	}
 
	## if there are NAs.
	if (any(is.na(nums))) {
		print(sprintf("There are %s NAs in the input, all of which will be removed prior to normality test, leaving %s non-NAs.", sum(is.na(nums)), sum(!is.na(nums))))
		nums <- nums[!is.na(nums)]
	}
	
	## select normality test.
	if (test %in% c("lilliefors", "anderson-darling")) {
		print(sprintf("%s test will be used for normality test.", test))
		norm.test <- switch(test, "lilliefors" = nortest::lillie.test, "anderson-darling" = nortest::ad.test)
	} else {
		stop(sprintf("%s is NOT valid, please use a valid normality test.", test))
	}
	
	## initial normality test, nortest::lillie.test, fBasics::dagoTest.
	normtest <- norm.test(nums)

	## There is no way to have a mixture of character and numeric class items in a matrix. (There is also no way of have a column in a data.frame which is a mixture of different atomic classes.)
	SigOutliers <- matrix(nrow = 0, ncol = 4, dimnames = list(NULL, c("id", "weight", "stat", "pvalue")))
	while (normtest$p.value <= pcut) {
		SigOutliers <- rbind(SigOutliers, c(names(nums[which.max(abs(nums))]), as.numeric(nums[which.max(abs(nums))]), normtest$statistic, normtest$p.value))
		
		## remove the elements with the largest absolute loading.
		premax <- nums[which.max(abs(nums))]
		nums <- nums[-which.max(abs(nums))]
		#invisible(ifelse(max(abs(nums)) == max(nums), nums <- nums[nums != max(abs(nums))], nums <- nums[nums != -max(abs(nums))]))
		#print(sprintf("%s elements left...", length(nums)))
		
		## normality test after removing the largest.
		normtest <- norm.test(nums)
		#if (normtest$p.value <= pcut) {SigOutliers <- rbind(SigOutliers, c(names(premax), as.numeric(premax), normtest$statistic, normtest$p.value))}
		
		if (length(nums) <= 4) {
			break
		}
	}
	
#	return(data.frame(SigOutliers, padj = p.adjust(SigOutliers[, "pvalue"], method = "BH")))
	return(data.frame(SigOutliers))
}


SpecIndex <- function(expprofile) {
## SpecIndex calculates expression specificity index ("τ" defined in 10.1093/bioinformatics/bti042) among tissues/cell types for the given expression profile.
## expprofile, the normalized expression profile for a feature as a vector.
	return(sum(1 - (expprofile/max(expprofile)))/(length(expprofile) - 1))
}


AddSpecIndex <- function(exptable) {
## AddSpecIndex calculates expression specificity index ("τ" defined in 10.1093/bioinformatics/bti042) among tissues/cell types for each feature in the input expression table, and add them as a column as "specidx".
## exptable, the normalized expression table where row names are features, column names are tissues/cell types.
	exptable <- as.data.frame(exptable)
	maxexp <- apply(exptable, 1, max)
	exptable_tmp <- 1 - sweep(x = exptable, MARGIN = 1, STATS = maxexp, FUN = "/")
	exptable$specidx <- apply(exptable_tmp, 1, sum)/(ncol(exptable) - 1)
	return(exptable)
}


#####################################################################################################################################################################################
#####################################################################################################################################################################################
## pseudobulkDE
simplepseudo <- function(seuobj, anncol) {
## similar function can be found here: https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R; https://github.com/MarioniLab/scran/blob/master/R/pseudoBulkDGE.R
## anncol: the column in the meta.data of the seuobj for cell type annotation, could be cell type or cluster, better as factor 
	
	## Quickly Creating Pseudobulks from Single Cell Gene Expression Data with Linear Algebra: https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
	meta <- seuobj@meta.data %>% dplyr::mutate(anninfo = as.character(!!(as.name(anncol)))) ## rlang::UQ(rlang::sym(anncol)); rlang::UQ(str2lang(anncol))
	mm <- model.matrix(~ 0 + anninfo, data = meta)
	## multiple row of mat.sparse (gene) onto column of the model matrix (cell-type annotation)
	simplepseudomat <- seuobj@assays$RNA@counts %*% mm ### Seurat::GetAssayData(seuobj, assay = "RNA", slot = "counts") %*% mm
	colnames(simplepseudomat) <- gsub("anninfo", "", colnames(simplepseudomat))
	return(simplepseudomat)
}



simplepseudoR <- function(seuobj, anncol, repcol) {
## similar function can be found here: https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R; https://github.com/MarioniLab/scran/blob/master/R/pseudoBulkDGE.R
## anncol: the column in the meta.data of the seuobj for cell type annotation, could be cell type or cluster, better as factor 
## repcol: the column in the meta.data of the seuobj for replicate information, e.g. "R1", "R2", better as factor
	
	## Quickly Creating Pseudobulks from Single Cell Gene Expression Data with Linear Algebra: https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
	meta <- seuobj@meta.data %>% dplyr::mutate(anninfo = as.character(!!(as.name(anncol))), repinfo = as.character(!!(as.name(repcol)))) ## rlang::UQ(rlang::sym(anncol)); rlang::UQ(str2lang(anncol))
	mm <- model.matrix(~ 0 + anninfo:repinfo, data = meta)
	## multiple row of mat.sparse (gene) onto column of the model matrix (cell-type annotation)
	simplepseudomatR <- seuobj@assays$RNA@counts %*% mm %>% as.data.frame() %>% rename_with(~ gsub("anninfo|repinfo", "", .x)) ### Seurat::GetAssayData(seuobj, assay = "RNA", slot = "counts") %*% mm
#	colnames(simplepseudomatR) <- gsub("anninfo|repinfo", "", colnames(simplepseudomatR)) ## simplepseudomatR %>% rename_with(~ gsub("anninfo|repinfo", "", .x))
	return(simplepseudomatR)
}




pseudobulks <- function(seuobj, anncol, repcol, concol){
## similar function can be found here: https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R; https://github.com/MarioniLab/scran/blob/master/R/pseudoBulkDGE.R
## anncol: the column in the meta.data of the seuobj for cell type annotation, could be cell type or cluster, better as factor 
## repcol: the column in the meta.data of the seuobj for replicate information, e.g. "R1", "R2", better as factor 
## concol: the column in the meta.data of the seuobj for condition information, e.g. "control", "treatment", better as factor 
	coninfo <- seuobj@meta.data %>% pull(!!as.name(concol)) %>% unique()
	
	pseudobulks <- coninfo %>% purrr::map( ~ {
		print(sprintf("Generating pseudobulks for %s", .))
		subseuobj <- subset(seuobj, subset = !!as.name(concol) == .)
		return(simplepseudoR(subseuobj, anncol, repcol))
	}) %>% setNames(coninfo)
	
	return(as.data.frame(do.call(cbind, pseudobulks)))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## pseudobulkDE to identify cell type markers
pseudoDE.findmarkers <- function(seuobj, anncol, repcol, ident, method = "edgeR", pcut = 0.05, logFC.cutoff = 0, params = list()) {
## similar function can be found here: https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R; https://github.com/MarioniLab/scran/blob/master/R/pseudoBulkDGE.R
## anncol: the column in the meta.data of the seuobj for cell type annotation, could be cell type or cluster, better as factor 
## repcol: the column in the meta.data of the seuobj for replicate information, e.g. "R1", "R2", better as factor
## ident, identity class in anncol to define markers for
	
	if (ident %in% seuobj@meta.data[[anncol]]) {
		seuobj@meta.data <- seuobj@meta.data %>% mutate(annR = ifelse((!!as.name(anncol)) == ident, ident, "others"))
	} else {
		stop(sprintf("Please make sure '%s' in '%s' in the metadata.", ident, anncol))
	}
	
	pseudobulks <- simplepseudoR(seuobj, "annR", repcol)

	if (tolower(method) == "edger"){
		library(edgeR)
		DE <- tryCatch({ 
			designtable <- data.frame(row.names = colnames(pseudobulks), groupR = colnames(pseudobulks)) %>% mutate(group = factor(gsub("\\:.*", "", groupR), levels = c("others", ident)), rep = gsub(".*\\:", "", groupR))
			DGEList <- DGEList(counts = pseudobulks, samples = designtable)
			DGEList <- DGEList[filterByExpr(DGEList, group = DGEList$samples$group), , keep.lib.sizes = FALSE]
			design <- model.matrix(~ group, data = designtable)
			
			# set defaults on params
			if (!"test" %in% names(params)) {
				params$test = "LRT"
			}
			if (!"method" %in% names(params)) {
				params$method = "TMM"
			}
			if (!"trend.method" %in% names(params)) {
				params$trend.method = "locfit"
			}
			if (!"robust" %in% names(params)) {
				params$robust = FALSE
			}
			if (!"tagwise" %in% names(params)) {
				params$tagwise = TRUE
			}
			if (!"prior" %in% names(params)) {
				params$prior = NULL
			}
			
			# catch some parameter type errors
			if (!is.null(params$prior)) {
				params$prior = as.numeric(params$prior)
			}
			if (params$tagwise == "FALSE") {
				params$tagwise = FALSE
			}
			
			DGEList = calcNormFactors(DGEList, method = params$method)
			if (params$robust) {
				DGEList <- estimateGLMRobustDisp(DGEList, design, trend.method = params$trend.method)
			} else {
				DGEList <- estimateDisp(DGEList, design, trend.method = params$trend.method, tagwise = params$tagwise, prior.df = params$prior)
			}
			
			if (params$test == "LRT") {
				fit <- glmFit(DGEList, design = design)
				test <- glmLRT(fit)
			} else {
				# QLF: default
				fit <- glmQLFit(DGEList, design)
				test <- glmQLFTest(fit, coef = -1)
			}
			res <- topTags(test, n = Inf) %>% as.data.frame() %>% rownames_to_column("gene") %>% mutate(group = ident, test = "pseudobulk_edgeR") %>% filter(FDR <= pcut, abs(logFC) >= logFC.cutoff)
			## "tryCatch" returns the value of the last expression evaluated, e.g. if the last expression is "print", the "print" message will be returned. 
		}, error = function(e) {
			message(e)
			return(data.frame())
		})
	}
	
	# DE %>% rename(any_of(c(
		# "p_val" = "p.value",  ## DESeq2
		# "p_val" = "p.value",  ## t/wilcox
		# "p_val" = "P.Value",  ## limma
		# "p_val" = "PValue"  , ## edgeR
		# "p_val_adj" = "padj", ## DESeq2/t/wilcox
		# "p_val_adj" = "adj.P.Val",      ## limma
		# "p_val_adj" = "FDR",            ## edgeER
		# "avg_logFC" = "log2FoldChange", ## DESEeq2
		# "avg_logFC" = "logFC" ## limma/edgeR
	# )))
	
	return(DE)
}


pseudoDE.findallmarkers <- function(seuobj, anncol, repcol, method = "edgeR", pcut = 0.05, logFC.cutoff = 0, params = list()) {
## similar function can be found here: https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R; https://github.com/MarioniLab/scran/blob/master/R/pseudoBulkDGE.R
## anncol: the column in the meta.data of the seuobj for cell type annotation, could be cell type or cluster, better as factor 
## repcol: the column in the meta.data of the seuobj for replicate information, e.g. "R1", "R2", better as factor
## ident, identity class in anncol to define markers for
	anncat <- seuobj@meta.data %>% pull(!!as.name(anncol)) %>% unique()
	
	DElist <- anncat %>% purrr::map( ~ {
		print(sprintf("Finding markers for group %s", .))
		pseudoDE.findmarkers(seuobj, anncol, repcol, ., method, pcut, logFC.cutoff, params)
	}) #%>% setNames(anncat)

	return(as.data.frame(do.call(rbind, DElist)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## pseudobulkDE to identify treatment responsive DEGs
simplepseudoDE <- function(seuobj, anncol, ident, concol, repcol, cond.1, cond.2, method = "edgeR", pcut = 0.05, logFC.cutoff = 0, params = list()) {
## pseudo.findDEG is to identify diffferential expression genes between conditions for a given cell type using pseudobulks (similar function can be found here: https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R; https://github.com/MarioniLab/scran/blob/master/R/pseudoBulkDGE.R)
## anncol: the column in the meta.data of the seuobj for cell type annotation, could be cell type or cluster, better as factor
## ident, identity class in anncol to identify DEGs for, could be cell type or cluster, better as factor
## concol: the column in the meta.data of the seuobj for condition annotation, better as factor
## repcol: the column in the meta.data of the seuobj for replicate information, e.g. "R1", "R2", better as factor
## cond.1, cond.2: condition class in concol to identify DEGs, the DE results will be cond.2/cond.1

#(only two levels allowed), e.g. "Ctrl", "Stress", better as factor if specific comparison order is desired (e.g. Stress/Ctrl or Ctrl/Stress) otherwise the default results will be levels(concol)[2]/levels(concol)[1]
	if (! all(c(anncol, concol, repcol) %in% colnames(seuobj@meta.data))) {
		stop(sprintf("Please make sure '%s', '%s', '%s' in the metadata.", anncol, concol, repcol))
	}
	
	
	if (ident %in% seuobj@meta.data[[anncol]]) {
		subseuobj <- subset(seuobj, subset = !!as.name(anncol) == ident)
	} else {
		stop(sprintf("Please make sure '%s' in '%s' in the metadata.", ident, anncol))
	}
	
	if (all(c(cond.1, cond.2) %in% subseuobj@meta.data[[concol]])) {
		subseuobj <- subset(subseuobj, subset = !!as.name(concol) %in% c(cond.1, cond.2))
	} else {
		stop(sprintf("Please make sure '%s' and '%s' in '%s' in the metadata.", cond.1, cond.2))
	}
	
	subseuobj@meta.data[[concol]] <- factor(subseuobj@meta.data[[concol]])
	count.tmp <- table(subseuobj@meta.data[[concol]])
	print(sprintf("Calculating DEGs for %s, the results will be %s (%s cells)/%s (%s cells)...", ident, cond.2, count.tmp[cond.2], cond.1, count.tmp[cond.1]))
	
	pseudobulks <- simplepseudoR(subseuobj, concol, repcol)

	if (tolower(method) == "edger"){
		library(edgeR)
		DE <- tryCatch({ 
			designtable <- data.frame(row.names = colnames(pseudobulks), conR = colnames(pseudobulks)) %>% mutate(condition = factor(gsub("\\:.*", "", conR), levels = c(cond.1, cond.2)), rep = gsub(".*\\:", "", conR))
			DGEList <- DGEList(counts = pseudobulks, samples = designtable)
			DGEList <- DGEList[filterByExpr(DGEList, group = DGEList$samples$condition), , keep.lib.sizes = FALSE]
			design <- model.matrix(~ condition, data = designtable)
			
			# set defaults on params
			if (!"test" %in% names(params)) {
				params$test = "LRT"
			}
			if (!"method" %in% names(params)) {
				params$method = "TMM"
			}
			if (!"trend.method" %in% names(params)) {
				params$trend.method = "locfit"
			}
			if (!"robust" %in% names(params)) {
				params$robust = FALSE
			}
			if (!"tagwise" %in% names(params)) {
				params$tagwise = TRUE
			}
			if (!"prior" %in% names(params)) {
				params$prior = NULL
			}
			
			# catch some parameter type errors
			if (!is.null(params$prior)) {
				params$prior = as.numeric(params$prior)
			}
			if (params$tagwise == "FALSE") {
				params$tagwise = FALSE
			}
			
			DGEList = calcNormFactors(DGEList, method = params$method)
			if (params$robust) {
				DGEList <- estimateGLMRobustDisp(DGEList, design, trend.method = params$trend.method)
			} else {
				DGEList <- estimateDisp(DGEList, design, trend.method = params$trend.method, tagwise = params$tagwise, prior.df = params$prior)
			}
			
			if (params$test == "LRT") {
				fit <- glmFit(DGEList, design = design)
				test <- glmLRT(fit)
			} else {
				# QLF: default
				fit <- glmQLFit(DGEList, design)
				test <- glmQLFTest(fit, coef = -1)
			}
			res <- topTags(test, n = Inf) %>% as.data.frame() %>% rownames_to_column("gene") %>% mutate(group = ident, comparison = paste0(cond.2, "/", cond.1), test = "pseudobulk_edgeR") %>% filter(FDR <= pcut, abs(logFC) >= logFC.cutoff)
			## "tryCatch" returns the value of the last expression evaluated, e.g. if the last expression is "print", the "print" message will be returned. 
		}, error = function(e) {
			message(e)
			return(data.frame())
		})
	}
	
	return(DE)
}


pseudo.findDEG <- function(seuobj, anncol, concol, repcol, cond.1, cond.2 = NULL, ident = NULL, method = "edgeR", pcut = 0.05, logFC.cutoff = 0, params = list()) {
## pseudo.findDEG is to identify diffferential expression genes between conditions for a given cell type using pseudobulks (similar function can be found here: https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R; https://github.com/MarioniLab/scran/blob/master/R/pseudoBulkDGE.R)
## anncol: the column in the meta.data of the seuobj for cell type annotation, could be cell type or cluster, better as factor
## concol: the column in the meta.data of the seuobj for condition annotation, better as factor
## repcol: the column in the meta.data of the seuobj for replicate information, e.g. "R1", "R2", better as factor
## cond.1, cond.2: condition class in concol to identify DEGs, the DE results will be cond.2/cond.1; if there are more than one cond.2, the DE results will be cond.n/cond.1 for each cond.2
## ident, identity class in anncol to identify DEGs for, could be cell type or cluster, better as factor

#(only two levels allowed), e.g. "Ctrl", "Stress", better as factor if specific comparison order is desired (e.g. Stress/Ctrl or Ctrl/Stress) otherwise the default results will be levels(concol)[2]/levels(concol)[1]

	if (is.null(ident)) {
		ident <- unique(seuobj@meta.data[[anncol]])
		print(sprintf("Calculating DEGs for %s ...", paste(ident, collapse = ",")))
	}
	if (is.null(cond.2)) {
		conds <- unique(seuobj@meta.data[[concol]])
		cond.2 <- conds[conds != cond.1]
		print(sprintf("DE results for %s will be %s/%s...", paste(ident, collapse = ","), paste(cond.2, collapse = ","), cond.1))
	}
	
	
	DElist <- purrr::map(cond.2, function(cond.n) {
		DElist.tmp <- purrr::map(ident, function(ident.n) {
			#print(sprintf("Finding markers for group %s", ident.n, cond.n))
			simplepseudoDE(seuobj, anncol, ident.n, concol, repcol, cond.1, cond.n, method, pcut, logFC.cutoff, params)
		})
		as.data.frame(do.call(rbind, DElist.tmp))
	}) %>% setNames(cond.2)
	
	return(DElist)
}





