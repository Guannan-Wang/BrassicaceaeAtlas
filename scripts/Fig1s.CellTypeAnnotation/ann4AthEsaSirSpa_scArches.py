#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse, sys, re, os, io, operator 
import warnings
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import matplotlib.pyplot as plt


def scarches_scanvi_annotator(adata_path, ref_key, ref_value, labels_key, model_affix, batch_key = "orig.ident", hvg_num = 2000):
	'''
	If your reference data is labeled (cell-type labels) and you have unlabeled or labeled query then use scArches scANVI
	ref_key, ref_value: "ref_value" in column "ref_key" from .obs will be selected as the reference dataset, ref_value should be list, e.g. ["study1"], ["study1", "study2"]
	labels_key: column in .obs for cell type annotations
	'''
	sc.set_figure_params(dpi_save = 300, frameon = False, figsize = (4, 4))

	adata = sc.read(adata_path)
	
	### 1. identify highly-variable genes
	adata_4HVGs = adata.copy()
	sc.pp.normalize_total(adata_4HVGs, target_sum = 1e4)
	sc.pp.log1p(adata_4HVGs)
	sc.pp.highly_variable_genes(adata_4HVGs, n_top_genes = hvg_num, batch_key = batch_key, subset = True)
	print("When n_top_genes = {}, {} highly-variable genes were selected.".format(hvg_num, len(adata_4HVGs.var_names)))
	
	### 1.1 save highly-variable genes
	with open(os.path.splitext(adata_path)[0] + "_" + str(hvg_num) + "HVGs.txt", "w") as hvg_fout:
		_ = hvg_fout.write(", ".join(adata_4HVGs.var_names))
	
	### 2. extract data for highly-variable genes
	adata_HVGs = adata[:, adata_4HVGs.var_names]
	ref_HVGs = adata_HVGs[adata_HVGs.obs[ref_key].isin(ref_value)].copy()
	query_HVGs = adata_HVGs[~adata_HVGs.obs[ref_key].isin(ref_value)].copy()
	
	import scarches as sca
	### 3. create SCANVI model and train it on fully labelled reference dataset
	sca.models.SCVI.setup_anndata(ref_HVGs, batch_key = None, labels_key = labels_key)
	ref_HVGs_vae = sca.models.SCVI(ref_HVGs, n_layers = 2, encode_covariates = True, deeply_inject_covariates = False, use_layer_norm = "both", use_batch_norm = "none")
	ref_HVGs_vae.train()
	ref_HVGs_scanvae = sca.models.SCANVI.from_scvi_model(ref_HVGs_vae, unlabeled_category = "Unknown") ### unlabeled_category: Value used for unlabeled cells in labels_key used to setup AnnData with scvi.
	
	print("For reference dataset: Labelled Indices: {}; Unlabelled Indices: {}".format(len(ref_HVGs_scanvae._labeled_indices), len(ref_HVGs_scanvae._unlabeled_indices)))
	
	ref_HVGs_scanvae.train(max_epochs = 20, use_gpu = False)
	ref_HVGs_scanvae.save(model_affix + "_ref_model/", overwrite = True) ## save the trained model, can be loaded with scarches.models.SCANVI.load("ref_model/")
	
	### 3.1. create anndata file of latent representation and compute UMAP
	ref_HVGs_latent = sc.AnnData(ref_HVGs_scanvae.get_latent_representation())
	ref_HVGs_latent.obs[labels_key] = ref_HVGs.obs[labels_key].tolist()
	ref_HVGs_latent.obs[batch_key] = ref_HVGs.obs[batch_key].tolist()
	ref_HVGs_latent.obs["predictions"] = ref_HVGs_scanvae.predict()
	print("Accuracy at n_top_genes = {}: {}".format(hvg_num, np.mean(ref_HVGs_latent.obs["predictions"] == ref_HVGs_latent.obs[labels_key])))
	
	sc.pp.neighbors(ref_HVGs_latent, n_neighbors = 10)
	sc.tl.leiden(ref_HVGs_latent)
	sc.tl.umap(ref_HVGs_latent)
	sc.pl.umap(ref_HVGs_latent, color = [batch_key, labels_key], frameon = False, wspace = 0.6, show = False, save = os.path.split(adata_path)[1].split(".")[0] + "_" + str(hvg_num) + "HVGs_ref.png")
	
	### 4. perform surgery on reference model and train on query dataset without cell type labels
	query_HVGs.obs["cCellType"] = ref_HVGs_scanvae.unlabeled_category_
	query_HVGs_model = sca.models.SCANVI.load_query_data(query_HVGs, model_affix + "_ref_model/", freeze_dropout = True)
	query_HVGs_model._unlabeled_indices = np.arange(query_HVGs.n_obs)
	query_HVGs_model._labeled_indices = []
	print("For query dataset: Labelled Indices: {}; Unlabelled Indices: {}".format(len(query_HVGs_model._labeled_indices), len(query_HVGs_model._unlabeled_indices)))
	
	query_HVGs_model.train(max_epochs = 100, plan_kwargs = dict(weight_decay = 0.0), check_val_every_n_epoch = 10, use_gpu = False)
	query_HVGs_model.save(model_affix + "_surgery_model/", overwrite = True)
	
	### 5. Get latent representation of reference + query dataset and compute UMAP
	full_HVGs = ref_HVGs.concatenate(query_HVGs)
	full_HVGs_latent = sc.AnnData(query_HVGs_model.get_latent_representation(adata = full_HVGs))
	full_HVGs_latent.obs[labels_key] = full_HVGs.obs[labels_key].tolist()
	full_HVGs_latent.obs[batch_key] = full_HVGs.obs[batch_key].tolist()
	full_HVGs_latent.obs["predictions"] = query_HVGs_model.predict(adata = full_HVGs)
	
	#### save predictions
	predictions = full_HVGs_latent.obs
	predictions.index = full_HVGs.obs.index
	predictions.to_csv(os.path.split(adata_path)[1].split(".")[0] + "_" + str(hvg_num) + "HVGs_scArchesAnnotaion.txt", sep = "\t", index = True)
	
	sc.pp.neighbors(full_HVGs_latent, n_neighbors = 10)
	sc.tl.leiden(full_HVGs_latent)
	sc.tl.umap(full_HVGs_latent)
	sc.pl.umap(full_HVGs_latent, color = [batch_key, labels_key, "predictions"], frameon = False, wspace = 0.6, save = os.path.split(adata_path)[1].split(".")[0] + "_" + str(hvg_num) + "HVGs_full.png")

###################################################################################################################################################################################


for num in list(range(2000, 16000, 2000)):
	try:
		print("### Processing n_top_genes = {}...".format(num))
		scarches_scanvi_annotator("AthEsaSirSpa.h5ad", ref_key = "orig.ident", ref_value = ["AthEsaSirSpa_AthCR1", "AthEsaSirSpa_AthCR2"], labels_key = "cCellType", model_affix = "scANVI" + str(num) + "HVGs", batch_key = "orig.ident", hvg_num = num)
		print("-" * 100)
	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		print("When n_top_genes = {}, error type {} occurs at line {}: {} ...".format(num, exc_type, exc_tb.tb_lineno, e))	
		continue







