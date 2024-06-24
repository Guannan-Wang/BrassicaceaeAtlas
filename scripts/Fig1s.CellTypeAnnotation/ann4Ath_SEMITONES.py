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



#######################################################################################################################################
AthMarkers = pd.read_excel("ArabidopsisCellTypeMarkers.xlsx", header = 0, index_col = None)


## Ath_Ctrls_RNAscaledata_escores.txt and Ath_Ctrls_RNAscaledata_pescores.txt were generated using run_SEMITONES.py
Ath_Ctrls_escores = pd.read_table("Ath_Ctrls_RNAscaledata_escores.txt", sep = "\t", header = 0, index_col = 0)
Ath_Ctrls_pescores = pd.read_table("Ath_Ctrls_RNAscaledata_pescores.txt", sep = "\t", header = 0, index_col = 0)

AthMarkers_noStele = AthMarkers[AthMarkers["CellTypeGroup"] != "Stele"] ## 454 rows × 5 columns
Ath_Ctrls_escores_wOnlyMarkers_noStele = Ath_Ctrls_escores.loc[[AthMarker_noStele for AthMarker_noStele in AthMarkersList_noStele if AthMarker_noStele in Ath_Ctrls_escores.index], :]


## Use different cutoffs for n_sds, and extract only the marker genes.
celltypepercell_noStele = []  ## create a list for cell types for each cells with different n_sds (based on the marker with the highest enrichment score).
MarkerCellType_noStele = dict(zip(AthMarkers_noStele["Gene"], AthMarkers_noStele["CellTypeGroup"]))

for nsds in [5, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50]:    
    ## sig_interval returns a dictionary {query cell: (lower, upper} (e.g. {'0': (-0.0338, 0.0342), '1': (-0.0352, 0.0356), ...}) of enrichment score significance cut-off below (lower) and above (upper) which the scores are significant at a certain standard deviation (n_sds) away from the mean of the permutation enrichment scores.
    escore_cutoffs = sig_interval(Ath_Ctrls_pescores, n_sds = nsds)
    
    ## identify the marker with the highest enrichment score for each cell.
    sigmarkerspercell = dict()
    for col in range(len(Ath_Ctrls_escores_wOnlyMarkers_noStele.columns)):
        sigmarkerspercell[Ath_Ctrls_escores_wOnlyMarkers_noStele.columns[col]] = [escore if escore_cutoffs[str(col)][0] < escore < escore_cutoffs[str(col)][1] else None for escore in list(Ath_Ctrls_escores_wOnlyMarkers_noStele.iloc[:, col])]
    maxmarkerpercell = pd.DataFrame(data = sigmarkerspercell, index = Ath_Ctrls_escores_wOnlyMarkers_noStele.index).idxmax(axis = 0, skipna = True)
    
    ## map signifcant markers with highest enrichment scores to their corresponding cell types. A cell is annotated as the label which its marker has the highest escore.
    maxmarkercelltype = dict()
    for cell in maxmarkerpercell.to_dict():
        maxmarkercelltype[cell] = MarkerCellType_noStele[maxmarkerpercell.to_dict()[cell]]    
    
    ## add the dictionaries, which containing cell type information for each cell at different n_sds cutoffs, to a list.
    print("Cell type (without markers for 'Stele') for each cell at n_sds = {} has been identified.".format(nsds))
    celltypepercell_noStele.append(maxmarkercelltype)


celltypes = pd.DataFrame(celltypepercell_noStele, index = ["marker_nsds" + str(nsds) for nsds in [5, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50]])
celltypes.transpose().to_csv("Ath_Ctrls_RNAscaledata_CellTypes_noStele.txt", sep = "\t", columns = None, header = True, index = True)



