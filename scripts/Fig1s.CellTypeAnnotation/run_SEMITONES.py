#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse, sys, re, os, io, operator, os.path
from decimal import Decimal ## This is for float comparison when accuracy is important!
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import anndata
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import pairwise_kernels
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.cell_selection import from_gui
from SEMITONES.cell_selection import from_2D_embedding
from SEMITONES.cell_selection import get_cells_from_gui
from SEMITONES.enrichment_scoring import calculate_escores
from SEMITONES.enrichment_scoring import permute, sig_interval
from SEMITONES.support_funcs import load_sparse_h5, save_sparse_h5, pairwise_similarities
from SEMITONES.support_funcs import sig_dictionary


synopsis = "\n\
#############################################################################################################################################\n\
#run_SEMITONES.py runs SEMITONES pipeline to identify highly enriched features in single-cell omics data without prior clustering.          #\n\
#                                                                                                                                           #\n\
#                                     Copyleft by wangguannan2014@gmail.com 20211006.                                                       #\n\
#############################################################################################################################################"


#####################################################################################################################################################################
def run_escores_calculation(scaledata, transpace, refcells = None, present = False, sig = True, ncpu = 1):
	"""
	run_escores_calculation
	scaledata (required), an array object where the columns are samples (i.e. cells) and the rows are features (i.e. genes). Accepts pandas dataframes, numpy arrays, and scipy compressed sparse row matrix, data should be normalized.
	transpace (required), a matrix-like object which contains the embedding of the original data in higher dimensions (e.g. 50D PC/umap), where rows are the cells and the columns are the dimensions.
	refcells (optional), a list of reference cells to be used in SEMITONES, given as the column indices in scaledata, default: None (all cells will be used).
	ncpu (optional), number of CPUs to be used, default: 1. 
	present (optional), whether the escores file is present, default: False
	sig (optional), whether perform significance testing, default: True
	"""
	sim = pairwise_similarities(transpace.to_numpy(), query = refcells, metric = "rbf", metric_params = {"gamma": 1.0/transpace.shape[1]})
	
	if not present:
		refcellids = scaledata.columns[refcells]
		escores = calculate_escores(scaledata.transpose(), query = refcellids, S = sim, scale_exp = False, ncpu = ncpu)
	else:
		escores = None
	
	if sig:
		n = 100 #if n is None else n
		seed = 42 #if seed is None else seed
		axis = 0 #if axis is None else axis
		pscaledata = permute(scaledata.transpose().values, n = n, seed = seed, axis = axis)
		pescores = calculate_escores(pscaledata, query = refcells, S = sim, scale_exp = False, ncpu = ncpu)
	else:
		pescores = None
		
	return escores, pescores
	
	

#####################################################################################################################################################################
#####################################################################################################################################################################
## Add arguments, "argcomplete" can be used to further manipulate the arguments.
## The difference between positional and optional arguments: "positional" and "--optional".
## "metavar" only changes the displayed name - the name of the attribute on the parse_args() object is still determined by the dest value.
## "dest" allows a custom attribute name to be provided for optional arguments, NOT for positional arguments.
## When "nargs" is specified, arguments from the command line will be gathered together into a list, e.g. nargs = 1 produces a list of one item. 
## Variable number of arguments can be set with the * character (* nargs expects 0 or more arguments, which will be gathered into a list). argparse: nargs='*' positional argument doesn't accept any items if preceded by an option and another positional. nargs = "+" like nargs = "*", but requires at least one value.
## "args" has all arguments as its attributes. All attributes can be verified with print(args.dest), e.g. print(args.Input), print(args.ClusteringMethod), print(args.Cc).
parser = argparse.ArgumentParser(description = synopsis, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("fscaledata", metavar = None, help = "an array file where the columns are samples (i.e. cells) and the rows are features (i.e. genes). Please make sure the data is in pandas dataframes. Data should be normalized.", action = "store", nargs = None, const = None, default = None, type = str, choices = None) ## Make this input as a list separated by comma
parser.add_argument("ftranspace", metavar = None, help = "a matrix-like file which contains the embedding of the original data in higher dimensions (e.g. 50D PC/umap), where rows are the cells and the columns are the dimensions.", action = "store", nargs = None, const = None, default = None, type = str, choices = None)
#parser.add_argument("fescores", metavar = None, help = "name for the output file.", action = "store", nargs = None, const = None, default = None, type = str, choices = None)
parser.add_argument("-frefcells", dest = None, help = "a file contain the list of reference cells to be used in SEMITONES, given as the column indices in fscaledata, default: None (all cells will be used).", action = "store", nargs = None, const = None, default = None, type = str, choices = None)
parser.add_argument("-ncpu", dest = None, help = "number of CPUs to be used, (default: %(default)s)", action = "store", nargs = None, const = None, default = "1", type = int, choices = None)

args = parser.parse_args() 

## for interactive prompt testing: passing a list of strings to parse_args(), e.g. parser.parse_args(["x", "y", "-out", "a,b,c,d,e"]) 	

if __name__ == "__main__":
	scaledata = pd.read_table(args.fscaledata, delimiter = "\t", header = 0, index_col = 0)
	print("The data in {} contains {} features and {} cells.".format(args.fscaledata, scaledata.shape[0], scaledata.shape[1]))
	transpace = pd.read_table(args.ftranspace, delimiter = "\t", header = 0, index_col = 0)
	print("The higher embedding in {} contains {} cells and {} dimensions.".format(args.ftranspace, transpace.shape[0], transpace.shape[1]))
	
	if args.frefcells:
		frefcells = open(args.frefcells, "r")
		refcells = [int(indice) for indice in frefcells.readlines()]
	else:
		refcells = np.arange(scaledata.shape[1])

	escores, pescores = run_escores_calculation(scaledata, transpace, refcells = refcells, present = os.path.isfile(os.path.splitext(args.fscaledata)[0] + "_escores.txt"), sig = True, ncpu = args.ncpu)	
	escores.to_csv(os.path.splitext(args.fscaledata)[0] + "_escores.txt", sep = "\t", columns = None, header = True, index = True) if escores is not None else None
	pescores.to_csv(os.path.splitext(args.fscaledata)[0] + "_pescores.txt", sep = "\t", columns = None, header = True, index = True)

		

