#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(tidyverse)
library(readxl)
library(Augur)


calcProbfromChildren <- function(p1, p2, l1 = 1,l2 = 1){
    (p1 * l2 + p2 * l1)/(l1 + l2)
}

#####################################################################################################################################################################################
Orth1v1_AthEsaSirSpa <- read.csv("Orth1v1_AthEsaSirSpa.csv", header = TRUE) %>% mutate(Sir = gsub("_", "-", Sir)) %>% column_to_rownames(var = "FinalOGs") %>% select(Ath:Spa) 
branchLengths <- c(SpaSir_Spa = 25.339, SpaSir_Sir = 25.339, SpaSirEsa_SpaSir = 3.157, SpaSirEsa_Esa = 28.496, SpaSirEsaAth_SpaSirEsa = 4.684, SpaSirEsaAth_Ath = 33.180) ## unit in Mya (million years ago)  ### Data extracted from Fig4(c) in Huang et al, Resolution of Brassicaceae Phylogeny Using Nuclear Genes Uncovers Nested Radiations and Supports Convergent Morphological Evolution, Molecular Biology and Evolution (DOI: 10.1093/molbev/msv226) using WebPlotDigitizer.


analysis.dir <- "path/to/analysis/dir"
norexp.dir <- "path/to/pseudobulk/expression/dir"
DE.dir <- "path/to/differential/expression/dir"

#####################################################################################################################################################################################
## Identify DEGs between cell types pairwisely using FindMarkers (only DEGs will be used for clustering)
Ath_wAnn <- readRDS("IntgratedAth.RDS")
Esa_wAnn <- readRDS("IntgratedEsa.RDS")
Sir_wAnn <- readRDS("IntgratedSir.RDS")
Spa_wAnn <- readRDS("IntgratedSpa.RDS")


for (speciesid in c("Ath", "Esa", "Sir", "Spa")){
    seuobj <- get(paste0(speciesid, "_wAnn")) %>% subset(., subset = intspace_treatment == "Ctrl")
    Idents(object = seuobj) <- "intspace_celltype"
    celltypes <- as.character(unique(Idents(object = seuobj)))
	print(sprintf("%s control sample has %s features and %s cells, %s cell types: %s.", speciesid, dim(seuobj)[1], dim(seuobj)[2], length(celltypes), paste(celltypes, collapse = ",")))
    DEres <- matrix(nrow = 0, ncol = 8, dimnames = list(NULL, c("geneid", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "species", "comp")))
    for (compair in apply(combn(celltypes, 2), 2, paste, collapse = "-")) { ## generate all pairwise combinations: https://stackoverflow.com/questions/48961247/how-do-i-create-all-unique-pairwise-combinations-of-my-sample-dataset
        pairDE <- tryCatch({
            ## Returns the value of the last expression evaluated. If the last expression is "print", the "print" message will be returned.
            pairDE_tmp <- FindMarkers(seuobj, assay = "RNA", slot = "data", ident.1 = strsplit(compair, "-")[[1]][1], ident.2 = strsplit(compair, "-")[[1]][2], test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, pseudocount.use = 1, min.diff.pct = -Inf, min.cells.group = 10) %>% rownames_to_column(var = "geneid") %>% mutate(species = speciesid, comp = compair)
        }, error = function(e) {
            message(e)
            pairDE_tmp <- matrix(nrow = 0, ncol = 8, dimnames = list(NULL, c("geneid", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "species", "comp")))
            return(pairDE_tmp)
        })
        print(sprintf("In %s: there are %s DEGs between %s.", speciesid, nrow(pairDE), compair))
        DEres <- rbind(DEres, pairDE)
    }
    print(sprintf("In %s: there are total %s DEGs between pairs of cell types.", speciesid, length(unique(DEres$geneid))))
    DEres %>% saveRDS(., file = paste0(speciesid, "Ctrls_CellType_pairDE.RDS"))
}

#####################################################################################################################################################################################
## Cluster the expression of the DEGs identified above using Mfuzz
## AthEsaSirSpa_pairwiseDE is a dataframe containing the expression of the DEGs identified above
AthEsaSirSpa_pairwiseDE_clstrs <- as.matrix(AthEsaSirSpa_pairwiseDE) %>% Biobase::ExpressionSet(assayData = .) %>% Mfuzz::standardise(.) %>% mfuzz(., c = 36, m = mestimate(.))

#####################################################################################################################################################################################
## Identify DEGs between cell types pairwisely using FindMarkers (only DEGs will be used for clustering)
AthEsaSirSpa_pairwiseDE_clstrmem <- as.data.frame(AthEsaSirSpa_pairwiseDE_clstrs$membership)

probs_df <- lapply(1:nrow(Orth1v1_AthEsaSirSpa_Sun), function(idx) {
    genes <- as.character(Orth1v1_AthEsaSirSpa_Sun[idx, ])
    if (all(genes %in% rownames(AthEsaSirSpa_pairwiseDE_clstrmem))) {
        gene_probs <- AthEsaSirSpa_pairwiseDE_clstrmem[genes, ]
        rownames(gene_probs) <- substr(rownames(gene_probs), 1, 2)
        #df <- reshape2::melt(as.matrix(gene_probs), varnames = c("node", "cluster"), value.name = "prob") %>% spread(key = node, value = prob)
        df <- t(gene_probs) %>% as.data.frame()
        
        ## calculate internal probabilities by weigthed average of its children
        df$SpaSir <- calcProbfromChildren(p1 = df$Sp, p2 = df$Si, l1 = branchLengths["SpaSir_Spa"], l2 = branchLengths["SpaSir_Sir"])
        df$SpaSirEsa <- calcProbfromChildren(p1 = df$SpaSir, p2 = df$Th, l1 = branchLengths["SpaSirEsa_SpaSir"], l2 = branchLengths["SpaSirEsa_Esa"])
        df$SpaSirEsaAth <- calcProbfromChildren(p1 = df$SpaSirEsa, p2 = df$AT, l1 = branchLengths["SpaSirEsaAth_SpaSirEsa"], l2 = branchLengths["SpaSirEsaAth_Ath"])
        
        df$OGid <- rownames(Orth1v1_AthEsaSirSpa_Sun)[idx] # genes[grepl("AT", genes)]
        df
    }
}) %>% bind_rows()

prob_same_df <- probs_df %>% group_by(OGid) %>% dplyr::summarise(Spa_Sir = sum(Sp * Si), SpaSir_Esa = sum(Th * SpaSir), SpaSirEsa_Ath = sum(AT * SpaSirEsa)) %>% gather(key = "comparison", value = "prob_same", -1)
prob_same_df_called <- filter(prob_same_df, prob_same < 0.01)
DfcallsPerGene <- prob_same_df_called %>% group_by(OGid) %>% dplyr::summarise(Spa_Sir = any(comparison %in% "Spa_Sir"), SpaSir_Esa = any(comparison %in% "SpaSir_Esa"), SpaSirEsa_Ath = any(comparison %in% "SpaSirEsa_Ath"), n_sig = length(comparison))
noBranchAssignableGenes <- filter(DfcallsPerGene, (Spa_Sir == SpaSir_Esa & Spa_Sir == TRUE) | (SpaSir_Esa == SpaSirEsa_Ath & SpaSir_Esa == TRUE))



## function to map a called change to one of its children by pairwise comparsion to "closest joint relative"
mapToBranch <- function(comp, gene4c, amb.gene.table, alpha = 0.01){
## comp: which comparison/node to look at
## gene4c: which gene should be used for the comparion
## amb.gene.table: a table for genes that cannot be assigned to a branch
    # genes that have called neighboring nodes cannot be used in a pairwise comparison due to changes in the line to comparsion node
    if(gene4c %in% amb.gene.table$OGid) {return(NA)}
    
    # on the branch to one of the following nodes a branch has happened:
    diff_nodes <- strsplit(comp, "_" )[[1]]
  
    # we use the other child of the parent of the called node as comparison
    # notes: 1) for changes between Arabidopsis and the rest we do not have a comparison node; 2) by the tree topology only end nodes are used in the comparison but could also be internal nodes
    comparisonNode <- ifelse(comp == "Spa_Sir", "Esa", ifelse(comp == "SpaSir_Esa", "Ath", ifelse(comp == "SpaSirEsa_Ath", "none", "wrong comparison")))
    stopifnot(comparisonNode != "wrong comparison")

    if(comparisonNode == "none") {
        return("last node")
        ## Try to include Ath branch
        # get probabilites of nodes for the comparsion
#         pp <- filter(probs_df, OGid == gene4c) %>% dplyr::rename(Ath = AT, Esa = Th, Sir = Si, Spa = Sp) #probs_df %>% dplyr::rename(Ath = AT, Esa = Th, Sir = Si, Spa = Sp)
#         prob1 <- sum(pp["SpaSirEsa"])
#         prob2 <- sum(pp["SpaSirEsaAth"])
#         changeAt <- ifelse(prob1 < alpha & prob2 > alpha, "SpaSirEsa", ifelse(prob1 > alpha & prob2 < alpha, "SpaSirEsaAth", ifelse(prob1 < alpha & prob2 < alpha, "last-both", "last-none")))
        return(changeAt)
    } else {
        # get probabilites of nodes for the comparsion
        pp <- filter(probs_df, OGid == gene4c) %>% dplyr::rename(Ath = AT, Esa = Th, Sir = Si, Spa = Sp) #probs_df %>% dplyr::rename(Ath = AT, Esa = Th, Sir = Si, Spa = Sp)
        prob1 <- sum(pp[diff_nodes[1]] * pp[comparisonNode])
        prob2 <- sum(pp[diff_nodes[2]] * pp[comparisonNode])

        changeAt <- ifelse(prob1 < alpha & prob2 > alpha, diff_nodes[1], ifelse(prob1 > alpha & prob2 < alpha, diff_nodes[2], ifelse(prob1 < alpha & prob2 < alpha, "both", "none")))
        return(changeAt)
    }
}


prob_same_df_called$on_branch <- sapply(1:nrow(prob_same_df_called), function(i) {mapToBranch(prob_same_df_called$comparison[i], prob_same_df_called$OGid[i], noBranchAssignableGenes)})


alpha <- 0.01
prob_same_df_called$pattern <- sapply(1:nrow(prob_same_df_called), function(i) {
    if(prob_same_df_called$OGid[i] %in% noBranchAssignableGenes$OGid) {
        return(NA)
    } else {
        pp <- filter(probs_df, OGid == prob_same_df_called$OGid[i]) %>% dplyr::rename(Ath = AT, Esa = Th, Sir = Si, Spa = Sp)
        AthEsa <- ifelse(sum(pp["Ath"]*pp["Esa"]) < alpha, 1, 0)
        AthSir <- ifelse(sum(pp["Ath"]*pp["Sir"]) < alpha, 1, 0)
        AthSpa <- ifelse(sum(pp["Ath"]*pp["Spa"]) < alpha, 1, 0)
        EsaSir <- ifelse(sum(pp["Esa"]*pp["Sir"]) < alpha, 1, 0)
        EsaSpa <- ifelse(sum(pp["Esa"]*pp["Spa"]) < alpha, 1, 0)
        SirSpa <- ifelse(sum(pp["Sir"]*pp["Spa"]) < alpha, 1, 0)
        paste0(AthEsa, "_", AthSir, "_", AthSpa, "_", EsaSir, "_", EsaSpa, "_", SirSpa) 
    }
})

#####################################################################################################################################################################################
q("no")

