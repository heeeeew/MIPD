import time
import os
import numpy as np
import pandas as pd
import pickle
import networkx as nx
from sklearn import metrics
from scipy.spatial.distance import pdist, squareform, expit, logit
from scipy.stats import ttest_1samp
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.model_selection import KFold
from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestClassifier
from . import prepare.py
from . import onesampttest.py
from . import DEGscore.py
from . import fold_genes.py
from . import Effect_Matrix.py
from . import train_model_test.py

DATA_PATH = ""
EM_PATH = ""
c_type = ""

eset_T = pd.read_csv(f"{DATA_PATH}{c_type}_eset_T.csv",index_col=0) 
eset_N = pd.read_csv(f"{DATA_PATH}{c_type}_eset_N.csv",index_col=0) 
mset = pd.read_csv(f"{DATA_PATH}{c_type}_mset_T.csv",index_col=0) 
CNV = pd.read_csv(f"{DATA_PATH}{c_type}_CNV.csv",index_col='Patient_barcode') 
methyl = pd.read_csv(f"{DATA_PATH}{c_type}_methylation.csv",index_col=0)  
edges = pd.read_csv(f"{DATA_PATH}network.csv")
KO = pd.read_csv(f"{DATA_PATH}/{c_type}_driver.csv")['Gene_symbol'].values

eset_T,eset_N,SM,edges,genes = prepare.prepare(eset_T, eset_N, mset, CNV, methyl, edges)
samples = eset_T.index.values
T_result = onesampttest.one_sample_t_test(eset_T,eset_N,genes)
DEGscore = DEGscore.DEGscore(genes,T_result)
train_genes, test_genes = fold_genes(genes, KO)
Effect_Matrix.EM(EM_PATH, eset_T, eset_N, edges, genes, DEGscore)
train_P, train_N = train_model_test.make_train_data(EM_PATH, train_genes, SM, genes, samples)
clfs = train_model_test.learning (train_P,train_N)
drivers = train_model_test.get_driver(genes,test_genes,clfs,samples)
