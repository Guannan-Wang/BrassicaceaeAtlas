#!/usr/bin/env Rscript

library(tidyverse)

###################################################################################################
analysis.dir <- "path/to/analysis/dir"
OrthDEG.logFC.unipattern <- readRDS(file.path(analysis.dir, "OrthDEG.logFC.unipattern.RDS"))


###################################################################################################
for (treatment in c("5uM_ABA", "100mM_NaCl")) {
	bootstrapping.res <- list() 
	for (iter in 1:10000) {
		print(sprintf("Running iteration %s for %s-treated samples...", iter, gsub("_", " ", treatment)))
		bootstrapping.res[[paste0("r", iter)]] <- Reduce(function(x, y, ...) {dplyr::full_join(x, y, by = "pattern", ...)}, 
			lapply(c("LRC", "Columella", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Pericycle", "Procambium", "Phloem", "Xylem"), function(celltypeid) {
				randomized.OrthDEGs <- OrthDEG.logFC.unipattern[[treatment]] %>% filter(celltype == celltypeid) %>% select(ends_with(".change")) 
				
				### shuffle the DEGs for each column while not perserving the order and position of NAs
				#randomized.OrthDEGs <- randomized.OrthDEGs %>% mutate(across(everything(), ~sample(.)))
				
				### shuffle the DEGs for each column while perserving the order and position of NAs
				for (species in c("Ath", "Esa", "Sir", "Spa")) {
					column <- paste0(species, ".", "change") 
					DEorder <- randomized.OrthDEGs[[column]]
					DEorder[!is.na(DEorder)] <- sample(DEorder[!is.na(DEorder)], replace = FALSE, prob = NULL)
					randomized.OrthDEGs[[column]] <- DEorder
				}
				
				### add a quantity checker here to see if the number of DEGs match every 100 iterations.
				if (iter %% 100 == 0) {
					lapply(c("Ath", "Esa", "Sir", "Spa"), function(speciesid) {
						table(randomized.OrthDEGs[[paste0(speciesid, ".change")]]) %>% as.matrix() %>% t() %>% as.data.frame() %>% mutate(species = speciesid, celltype = celltypeid)
					}) %>% data.table::rbindlist(fill = TRUE) %>% column_to_rownames(var = "species") %>% print()
				}
				
				### count the frequencies of different patterns.
				full_join(randomized.OrthDEGs %>% replace(is.na(.), 0) %>% unite("pattern", Ath.change:Spa.change, sep = "_", remove = TRUE, na.rm = FALSE) %>% group_by(pattern) %>% dplyr::summarize(!!paste0(celltypeid, ".wNA.counts") := n()), randomized.OrthDEGs %>% drop_na() %>% unite("pattern", Ath.change:Spa.change, sep = "_", remove = TRUE, na.rm = FALSE) %>% group_by(pattern) %>% dplyr::summarize(!!paste0(celltypeid, ".woNA.counts") := n()), by = "pattern")
			})
		)
	}
	saveRDS(bootstrapping.res, file = paste0("OrthDEG_", treatment, "_bootstrapping.RDS"))
	print(strrep("-", 100))
}

print("Completed bootstrapping!")

q("no")



