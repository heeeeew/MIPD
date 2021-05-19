import argparse
import time
import numpy as np
import pandas as pd
import networkx as nx
import os
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ttest_1samp, ttest_ind, ranksums, ttest_rel, wilcoxon
from scipy.special import expit, logit
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from sklearn import metrics
from sklearn.utils import shuffle
import pickle

def sigmoid(x):
    return 1 / (1 +np.exp(-x))

def EM(eset_T, eset_N, edges, genes, DEGscore):

    print(SM.shape)
    print(f"normal : {eset_N.shape}")
    print(f"tumor : {eset_T.shape}")
    print(f"genes : {genes.shape}")
    print(f"edges : {edges.shape}") 

    PCC_N = 1. - squareform(pdist(eset_N.T, 'correlation'))
    PCC_T = 1. - squareform(pdist(eset_T[genes].T, 'correlation'))
    avg = eset_T[genes].mean(axis=0).values.reshape(-1)
    std = eset_T[genes].std(axis=0).values.reshape(-1)
    G = nx.DiGraph()
    G.add_edges_from(edges)    
    A = nx.to_numpy_array(G, nodelist=genes)
    all_W = A * np.abs(PCC_T) * (0.5*np.abs(PCC_T - PCC_N))
    samples = np.array(eset_T.index)
    print('start')
    for test_index in range(0,eset_T.shape[0]):
        start = time.time()
        u = eset_T[genes].iloc[test_index].values
        save_name = samples[test_index]
        P_barcode = save_name[:12]
        
        ## Weighted adjacency matrix for PageRank
        start2 = time.time()

        nor_vector = ((u - avg)/std).reshape(-1,1)
        ERR = expit(np.dot(nor_vector,nor_vector.T) * np.sign(PCC_T))

        W = all_W * ERR
        W_denominator = W.sum(axis=1).reshape(-1)
        W_denominator[W_denominator==0] = 1.
        
        # dest, src
        W = W.T / W_denominator
        W_sum = W.sum(axis=0)   
        D = DEGscore.loc[save_name].values
        D = D.reshape(-1)
        D[W_sum == 0] = 1   
        W = (W*(1-D)) + np.diag(D.reshape(-1))
        
        ## PageRank
        Init_score = 1e4
        Epoch = 5
      
        S2 = np.ones(len(genes))
        I = np.diag(Init_score*S2)
        
        for _ in range(Epoch): 
            I = np.dot(W, I)
        I = pd.DataFrame(I, index = genes, columns = genes)
        
        with open(OUTPUT_PATH + save_name+".pkl",'wb') as f:
            pickle.dump(I,f) 
           
        end = time.time()
        print(f"[{test_index}] : {end-start}sec")
