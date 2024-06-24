#!/usr/bin/env Rscript

library(Matrix)
library(Seurat)
library(tidyverse)
library(readxl)
library(Augur)


source("scRNAseq_funs.R")

RowSD <- function(x, rmNAs = TRUE) {
    if(isTRUE(rmNAs)) {
        sqrt(rowSums((x - rowMeans(x, na.rm = rmNAs))^2, na.rm = rmNAs)/(sum(!is.na(x)) - 1))
    } else (
        sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
    )
}

#####################################################################################################################################################################################
### from Chaudhary et al., 10.1534/g3.119.400957
redefined.subG <- list(Csa.subG1 = c("Csa14", "Csa07", "Csa19", "Csa04", "Csa08", "Csa11"), Csa.subG2 = c("Csa03", "Csa16", "Csa01", "Csa06", "Csa13", "Csa10", "Csa18"), Csa.subG3 = c("Csa17", "Csa05", "Csa15", "Csa09", "Csa20", "Csa02", "Csa12"))

Csa.refined.triplets <- readxl::read_excel("Chaudhary et al._Camelina sativa_triplets.xlsx", sheet = "TableS8", skip = 1) %>% dplyr::rename("AthID" = "ArabidopsisGeneIdentifier", "Csa.subG1" = "Cs-G1", "Csa.subG2"  = "Cs-G2", "Csa.subG3"  = "Cs-G3") %>% mutate(Csa.subG1 = gsub("\\..*", "", Csa.subG1), Csa.subG2 = gsub("\\..*", "", Csa.subG2), Csa.subG3 = gsub("\\..*", "", Csa.subG3)) %>% select(-Cneglecta )%>% dplyr::mutate(num.homeologs = rowSums(across(all_of(c("Csa.subG1", "Csa.subG2", "Csa.subG3")), ~ .x %in% c("x")))) %>% dplyr::mutate(num.homeologs = 3 - num.homeologs) %>% filter(num.homeologs != 0) %>% mutate(subG1.cats = ifelse(Csa.subG1 == "x", NA, ifelse(num.homeologs == 1, "singlet", ifelse(num.homeologs == 2, "doublet", "triplet"))), subG2.cats = ifelse(Csa.subG2 == "x", NA, ifelse(num.homeologs == 1, "singlet", ifelse(num.homeologs == 2, "doublet", "triplet"))), subG3.cats = ifelse(Csa.subG3 == "x", NA, ifelse(num.homeologs == 1, "singlet", ifelse(num.homeologs == 2, "doublet", "triplet"))))


paml.dir <- "path/to/paml/dir"
norexp.dir <- "path/to/pseudobulk/expression/dir"
DE.dir <- "path/to/differential/expression/dir"

#####################################################################################################################################################################################
Csa.pseudoExp <- readRDS(file.path(norexp.dir, "220812_Csa_NorExp.RDS"))
Csa.norexp <- Csa.pseudoExp[["Csa_pseudoCT"]]



complete.triplets.exp <- Chaudhary.complete.triplets %>% pivot_longer(cols = contains("Csa.subG"), names_to = "Csa.subG", values_to = "geneID") %>% filter(geneID != "x", num.homeologs == 3) %>% 
    left_join(., Csa.norexp %>% rownames_to_column(var = "geneID"), by = "geneID") %>%
    left_join(., Csa.pseudoSum.vst %>% dplyr::rename_with(~ gsub(":", "_", .x)) %>% rownames_to_column(var = "geneID"), by = "geneID") %>% 
    select(HomeologID, AthID, num.homeologs, Csa.subG, geneID, starts_with("CsaC_"), ends_with("_Ctrl"))

homeolog.expDiv <- lapply(c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Phloem", "Xylem"), function(celltypeid) {
    complete.triplets.exp %>% select(HomeologID, AthID, num.homeologs, Csa.subG, geneID, starts_with(paste0(celltypeid, "_"))) %>% drop_na(any_of(c(paste0(celltypeid, "_Ctrl")))) %>% dplyr::rename(c("vstExp" = paste0(celltypeid, "_Ctrl"))) %>% 
    transform(rev.num.HemG = ave(geneID, HomeologID, FUN = length)) %>% filter(rev.num.HemG == 3) %>% 
    select(HomeologID, AthID, Csa.subG, geneID, vstExp) %>%
    pivot_wider(names_from = "Csa.subG", values_from = c("geneID", "vstExp"), names_glue = "{Csa.subG}_{.value}") %>% 
    dplyr::mutate(vstExpAvg = rowMeans(select(., ends_with("vstExp"))), vstExpSD = RowSD(select(., ends_with("vstExp")))) %>%
    dplyr::mutate(vstExpCV = vstExpSD/vstExpAvg, celltype = celltypeid)
}) %>% data.table::rbindlist() %>% filter(vstExpCV != 0)






Csa.expSpe <- Csa.norexp %>% select(starts_with("CsaC_")) %>% AddSpecIndex(.) %>% dplyr::rename(c("CsaC.specidx" = "specidx"))

complete.triplets.expSpe <- Chaudhary.complete.triplets %>% pivot_longer(cols = contains("Csa.subG"), names_to = "Csa.subG", values_to = "geneID") %>% filter(geneID != "x", num.homeologs == 3) %>% 
    left_join(., Csa.expSpe %>% rownames_to_column(var = "geneID"), by = "geneID") %>%
    select(HomeologID, AthID, num.homeologs, Csa.subG, geneID, CsaC.specidx)

complete.triplets.allspeDiv <- Chaudhary.complete.triplets %>% pivot_longer(cols = contains("Csa.subG"), names_to = "Csa.subG", values_to = "geneID") %>% filter(geneID != "x", num.homeologs == 3) %>% 
    left_join(., Csa.expSpe %>% filter(if_any(starts_with("CsaC_"), ~.x > 0.01)) %>% rownames_to_column(var = "geneID"), by = "geneID") %>% ## exp > 0.1 (4783); exp > 0.01 (11500)
    select(HomeologID, AthID, num.homeologs, Csa.subG, geneID, CsaC.specidx) %>% drop_na(CsaC.specidx) %>%
    dplyr::rename(c("tau" = "CsaC.specidx")) %>%
    transform(rev.num.HemG = ave(geneID, HomeologID, FUN = length)) %>% filter(rev.num.HemG == 3) %>%
    select(HomeologID, AthID, Csa.subG, geneID, tau) %>%
    pivot_wider(names_from = "Csa.subG", values_from = c("geneID", "tau"), names_glue = "{Csa.subG}_{.value}") %>%
    dplyr::mutate(tauAvg = rowMeans(select(., ends_with("tau"))), tauSD = RowSD(select(., ends_with("tau")))) %>%
    dplyr::mutate(tauCV = tauSD/tauAvg)
	





Csa.DE <- readRDS(file.path(DE.dir, "220816_Csa_pseudobulkDE.RDS"))

Csa.fullDE <- full_join(
    Csa.DE[["5uM_ABA"]] %>% select(all_of(c("gene", "group", "logFC", "FDR"))) %>% pivot_wider(names_from = "group", values_from = c("logFC", "FDR"), names_glue = "{group}.{.value}") %>% column_to_rownames(var = "gene") %>% dplyr::rename_with(~ paste0("ABA.", .x)) %>% rownames_to_column(var = "geneID"), 
Csa.DE[["100mM_NaCl"]] %>% select(all_of(c("gene", "group", "logFC", "FDR"))) %>% pivot_wider(names_from = "group", values_from = c("logFC", "FDR"), names_glue = "{group}.{.value}") %>% column_to_rownames(var = "gene") %>% dplyr::rename_with(~ paste0("NaCl.", .x)) %>% rownames_to_column(var = "geneID") %>% select(-starts_with("NaCl.Endodermis."), -starts_with("NaCl.Phloem."), -starts_with("NaCl.Xylem.")),
    by = "geneID"
)

complete.triplets.resDiv <- lapply(c("ABA", "NaCl"), function(treatmentid) {
    lapply(c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Phloem", "Xylem"), function(celltypeid){
        if (paste0(treatmentid, ".", celltypeid, ".logFC") %in% colnames(Csa.fullDE)) {
            Chaudhary.complete.triplets %>% pivot_longer(cols = contains("Csa.subG"), names_to = "Csa.subG", values_to = "geneID") %>% filter(geneID != "x", num.homeologs == 3) %>% 
            left_join(., Csa.fullDE %>% select(geneID, starts_with(paste0(treatmentid, ".", celltypeid, "."))), by = "geneID") %>%
            select(HomeologID, AthID, num.homeologs, Csa.subG, geneID, starts_with(paste0(treatmentid, ".", celltypeid, "."))) %>%
            ## adding a filter to remove homeologous groups where none of the three members are differentially expressed. 
            group_by(HomeologID) %>% filter(any(.data[[paste0(treatmentid, ".", celltypeid, ".FDR")]] <= 0.01)) %>%
            transform(rev.num.HemG = ave(geneID, HomeologID, FUN = length)) %>% filter(rev.num.HemG == 3) %>%
            ## replacing log2FC of NAs to 0
            replace(is.na(.), 0) %>% 
            group_by(HomeologID) %>%
            dplyr::summarise(resAvg = mean(.data[[paste0(treatmentid, ".", celltypeid, ".logFC")]]), resSD = sd(.data[[paste0(treatmentid, ".", celltypeid, ".logFC")]])) %>% 
            arrange(desc(resSD)) %>% 
            mutate(treatment = treatmentid, celltype = celltypeid)
        }
    }) %>% data.table::rbindlist()
}) %>% data.table::rbindlist()







### HmeG.omegas contains omega values calculated for each homelog group using codeML
ranked.div <- full_join(
    HmeG.omegas %>% filter(omega >= quantile(omega, probs = c(0.99), na.rm = TRUE)) %>% arrange(desc(omega)) %>% dplyr::mutate(top.omega = 1:n()) %>% select(HmeG, top.omega),
    homeolog.expDiv %>% filter(expCV != 0) %>% arrange(desc(expCV)) %>% filter(expCV >= quantile(expCV, probs = c(0.99), na.rm = TRUE)) %>% dplyr::mutate(top.expCV = 1:n()) %>% select(HomeologID, celltype, top.expCV) %>% pivot_wider(names_from = "celltype", values_from = "top.expCV", names_glue = "{celltype}.{.value}") %>% select(HomeologID, all_of(paste0(c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Phloem", "Xylem"), ".top.expCV"))),
    by = c("HmeG" = "HomeologID")
) %>% 
full_join(., complete.triplets.allspeDiv %>% filter(speCV.exp0.01 >= quantile(speCV.exp0.01, probs = c(0.99), na.rm = TRUE)) %>% dplyr::mutate(top.speCV = 1:n()) %>% select(HomeologID, top.speCV),
    by = c("HmeG" = "HomeologID")
) %>% 
full_join(., complete.triplets.resDiv %>% filter(treatment == "ABA") %>% arrange(desc(resSD)) %>% filter(resSD >= quantile(resSD, probs = c(0.99), na.rm = TRUE)) %>% dplyr::mutate(top.resSD = 1:n()) %>% select(HomeologID, celltype, top.resSD) %>% pivot_wider(names_from = "celltype", values_from = "top.resSD", names_glue = "{celltype}.{.value}") %>% dplyr::rename_with(~ paste0("ABA.", .x), ends_with("top.resSD")) %>% select(HomeologID, any_of(paste0("ABA.", c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Phloem", "Xylem"), ".top.resSD"))),
    by = c("HmeG" = "HomeologID")
) %>%
full_join(., complete.triplets.resDiv %>% filter(treatment == "NaCl") %>% arrange(desc(resSD)) %>% filter(resSD >= quantile(resSD, probs = c(0.99), na.rm = TRUE)) %>% dplyr::mutate(top.resSD = 1:n()) %>% select(HomeologID, celltype, top.resSD) %>% pivot_wider(names_from = "celltype", values_from = "top.resSD", names_glue = "{celltype}.{.value}") %>% dplyr::rename_with(~ paste0("NaCl.", .x), ends_with("top.resSD")) %>% select(HomeologID, any_of(paste0("NaCl.", c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Phloem", "Xylem"), ".top.resSD"))),
    by = c("HmeG" = "HomeologID")
) %>% column_to_rownames(var = "HmeG") %>% dplyr::mutate(num.nas = rowSums(is.na(select(., everything()))), num.occ = rowSums(!is.na(select(., everything())))) %>% rownames_to_column(var = "HmeG")

#####################################################################################################################################################################################
q("no")

