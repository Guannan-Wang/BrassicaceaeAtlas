#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(tidyverse)
library(readxl)
library(Augur)

source("scRNAseq_funs.R")


#####################################################################################################################################################################################
## orthDEGs_AggregatedExpSet.RDS contains expression for the expression of all orthologs under different conditions
## orthDEGs_AggregatedExpSet.mb contains aerie-derived cluster categorization for all orthologs

stdnorexp.table <- readRDS(file.path(analysis.dir, paste0("orthDEGs_AggregatedExpSet.RDS"))) %>% Biobase::exprs() %>% as.data.frame()
cluster.table <- read.table(file.path(analysis.dir, paste0("orthDEGs_AggregatedExpSet.mb")), sep = "\t", header = TRUE, fill = TRUE) %>% select(UID, starts_with("Centroid.")) %>% select(-last_col()) %>% column_to_rownames(var = "UID") %>% select(starts_with("Centroid.")) %>% mutate(max.membership = do.call(pmax, c(select(., everything()), na.rm = TRUE)), cluster = names(.)[max.col(.)]) #%>% filter(max.membership >= 0.5)


count.list <- list()
for (clstrid in 0:(length(unique(cluster.table$cluster)) - 1)) {
    genes.used <- cluster.table %>% filter(cluster == paste0("Centroid.", clstrid), max.membership >= 0.5) %>% rownames()
    stdnorexp.tmp <- stdnorexp.table[genes.used, ]

    ### check correlation of individual member to the means of the cluster
    if (length(genes.used) > 0)  {
        genes.filtered <- lapply(rownames(stdnorexp.tmp), function(OGID) {
            cor.tmp <- cor.test(as.numeric(stdnorexp.tmp[OGID, ]), colMeans(stdnorexp.tmp), method = "pearson")
            data.frame(row.names = OGID, "cor" = as.numeric(cor.tmp$estimate), "pvalue" = cor.tmp$p.value) 
        }) %>% dplyr::bind_rows() %>% filter(cor >= 0.8, pvalue < 0.01) %>% rownames()
        
        ### counts of orthDEGs from different cell types for each cluster
        count.list[[paste0("C", clstrid)]] <- stdnorexp.table[genes.filtered,] %>% rownames_to_column(var = "UID") %>% mutate(treatment = gsub("\\..*", "", UID), celltype = word(UID, 2, sep = "\\."), OGID = gsub(".*\\.", "", UID)) %>% group_by(treatment, celltype) %>% dplyr::summarize(counts = n()) %>% dplyr::rename(!! paste0("C", clstrid, ".counts") := "counts")
        print(sprintf("For aggregated orthDEGs: there are %s clusters; for cluster %s, there are %s members with membership cutoff >= 0.5, %s left after filtering.", length(unique(cluster.table$cluster)) - 1, clstrid, length(genes.used), length(genes.filtered)))
    } else {
        print(sprintf("For aggregated orthDEGs: there are %s clusters; for cluster %s, there are %s members with membership cutoff >= 0.5, none left after filtering.", length(unique(cluster.table$cluster)) - 1, clstrid, length(genes.used)))
    }
}



clstr.countable <- Reduce(function (x, y, ...) {full_join(x, y, by = c("treatment", "celltype"))}, count.list) %>% replace(is.na(.), 0)

total <- clstr.countable %>% ungroup() %>% select(ends_with(".counts")) %>% sum(., na.rm = TRUE)
ABA.total <- clstr.countable %>% filter(treatment == "ABA") %>% ungroup() %>% select(ends_with(".counts")) %>% sum(., na.rm = TRUE)
mMNa.total <- clstr.countable %>% filter(treatment == "mMNa") %>% ungroup() %>% select(ends_with(".counts")) %>% sum(., na.rm = TRUE)

clstr.pertable <- clstr.countable %>% ungroup() %>% dplyr::mutate(celltype.total = rowSums(select(., ends_with(".counts")), na.rm = TRUE)) %>% group_by(treatment) %>% dplyr::mutate(across(ends_with(".counts"), ~ sum(., na.rm = TRUE), .names = "{col}.sum"), total = ifelse(treatment == "ABA", ABA.total, mMNa.total)) %>% dplyr::mutate(across(ends_with(".counts"), ~ . / get(paste0(cur_column(), ".sum")) * 100, .names = "{col}.clstrper"), background.counts.clstrper = celltype.total/sum(celltype.total, na.rm = TRUE) * 100)

hg.res <- clstr.pertable %>% dplyr::mutate(across(ends_with(".counts"), ~ phyper(., celltype.total, total - celltype.total, get(paste0(cur_column(), ".sum")), lower.tail = FALSE), .names = "{col}.pvalue")) %>% group_by(treatment) %>% dplyr::mutate(across(ends_with(".pvalue"), ~ p.adjust(., method = "BH"), .names = "{col}.adj"), background.counts.pvalue.adj = 1) %>% select(treatment, celltype, ends_with(".pvalue.adj"))

#####################################################################################################################################################################################
q("no")

