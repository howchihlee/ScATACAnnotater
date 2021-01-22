#!/usr/bin/env python
# coding: utf-8
import os
import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc

import SortedNoDupeBedOverlap as bdO
import stats_utils
import io_utils 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('count_input', type=str, help='path to 10x input folder')
    parser.add_argument('ref_folder', type=str, help='path to a folder of cell type specific peaks')
    parser.add_argument('-o', '--output', default = 'score.csv', help='name of the output file')
    parser.add_argument('-t', '--threshold', type=float, default = 0.5, action='store', 
                        help='threshold to compute bed overlap. a float number ranges from 0 - 1.')
    parser.add_argument('--cpu', type=int, default = 1, action='store', help='number of cpus')
    parser.add_argument('--verbose', action='store_true', help='flag to print messages')
    ## for plotting
    parser.add_argument('--umap_fig', action='store_true', help='flag to skip ploting umap')
    parser.add_argument('-u', '--umap_input', default = '', type=str, 
                        help='path to a csv file of umap coordinates (cell x 2). The file  should have an index column of cell barcodes')
    parser.add_argument('-f', '--output_fig', default = 'score.png', help='name of the output png figure')
    
    configs = parser.parse_args()
        
    return configs



if __name__ == "__main__":      
    #print(args.accumulate(args.integers))
    #check_input(configs)
    configs = parse_args()
    
    count_input = configs.count_input
    verbose = configs.verbose
    ref_folder = configs.ref_folder
    threshold = configs.threshold
    n_cpu = configs.cpu
    output_file = configs.output
    
    mat, id2bc, id2peak = io_utils.read_atac_count_10x(count_input, verbose = verbose)
    ref_subtype_peaks, ref_bk_peak = io_utils.read_reference_data(ref_folder, verbose = verbose)
    
    if verbose: print("Calculating Bed Overlaps")
    bk_intersect = set(bdO.BedOverlap(sorted(id2peak), sorted(ref_bk_peak), threshold))
    
    ref_subtype_peaks_intersect = []
    for c, peaks in ref_subtype_peaks:
        intersect = bdO.BedOverlap(sorted(id2peak), sorted(peaks), threshold)
        ref_subtype_peaks_intersect.append( (c, set( intersect )) )
    id2ct = [c for c, _ in ref_subtype_peaks_intersect ]
    
    if verbose: 
        print("Calculating Fischer Exact Scores")
    scores = stats_utils.compute_enrichment_score(mat, bk_intersect, id2peak, ref_subtype_peaks_intersect, num_cores = n_cpu)

    if verbose: print("Outputting scores")
    score_df = io_utils.make_dataframe(scores, index = id2bc, columns = id2ct) 
    score_df.to_csv(output_file)
    
    if not configs.umap_fig:
        umap_input = configs.umap_input
        output_fig = configs.output_fig
        
        if len(umap_input) > 0:
            df_umap = pd.read_csv(umap_input, index_col = 0)
            data_obj = sc.AnnData(df_umap, 
                  obs = pd.DataFrame([], index = df_umap.index),
                  var = pd.DataFrame([], index = df_umap.columns),
                  )  
            data_obj.obsm['X_umap'] = df_umap.values   

        else:
            
            ## pca implementation using sparse matrix to save memory
            pcas = stats_utils.normalized_pca_from_sparse_mat(mat, n_pc = 20) 
            
            data_obj = sc.AnnData(pcas, 
                  obs = pd.DataFrame([], index = id2bc),
                  )        

            sc.pp.neighbors(data_obj, use_rep='X')
            sc.tl.umap(data_obj) 

        data_obj.obs = score_df

        sc.pl.umap(data_obj, color = score_df.columns, cmap = 'RdBu_r', vmax = 5, show = False)
        plt.savefig(output_fig)
    
    sys.exit(0)
