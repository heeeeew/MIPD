import numpy as np
import pandas as pd
import pickle
from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestClassifier
import time

def make_train_data (EM_PATH, train_genes, SM, genes, samples):
    N_column = {}
    P_column = {}

    for fold in range(5):
        N_column[fold] = pd.DataFrame()
        P_column[fold] = pd.DataFrame()
     
    i = 0    
     
    for samp in samples:
        start = time.time()
        P_barcode = samp[:12]
        
        with open(EM_PATH +samp+".pkl",'rb') as f:
            I = pickle.load(f).loc[genes,genes] 
        somatic = SM.loc[P_barcode].values

        for fold in range(5):
            P_gene = set(train_genes[fold]['KO']) & set(genes[somatic==1])
            N_gene = shuffle(train_genes[fold]['notKO'])[:len(P_gene)]

            N_column[fold] = pd.concat([N_column[fold],I[N_gene]],axis=1)       
            P_column[fold] = pd.concat([P_column[fold],I[P_gene]],axis=1)
           
        end = time.time()
        print(f"[{i}] {P_barcode} : {end - start} sec")
        i += 1
    return train_P,train_N

    
def learning (train_P,train_N):
    model = {}
    for fold in range(5):
        start = time.time()
        train_P = P[fold].values.T
        train_N = N[fold].values.T  
        print(f"P {train_P.shape}, N {train_N.shape}")
        
        train_data = np.concatenate([train_P,train_N])
        train_label = np.concatenate([np.zeros(train_P.shape[0]),np.zeros(train_N.shape[0])+1])
        train_data,train_label = shuffle(train_data,train_label)
        print(f"data {train_data.shape} label {train_label.shape}")
        
        
        clf = RandomForestClassifier(random_state=119)
        clf.fit(train_data,train_label)
        model[fold] = clf
        end = time.time() 
        print(f"{fold} : {end - start}sec")
    return model

    
def get_driver(genes,test_genes,clf,samples):
    tmp_gene = {}
    for fold in range(5):
        tmp_gene[fold] = np.array(sorted(list(set(test_genes[fold]['KO']) | set(test_genes[fold]['notKO']))))
    test_genes = tmp_gene
        
    P_genes = {}
    i = 0

    for samp in samples:
        start = time.time()
        P_barcode = samp[:12]
        tmp_gene = set()
        with open(EM_PATH +samp+".pkl",'rb') as f:
            I = pickle.load(f).loc[genes,genes] 
        
        for fold in range(5):
            test_gene = test_genes[fold]
            test_data = I[test_gene].T.values
            pred = clf[fold].predict(test_data)
            tmp_gene = tmp_gene | set(test_gene[pred==0])
        P_genes[samp] = sorted(list(tmp_gene))
        end = time.time()
        print(f"[{i}] {end - start}sec")
        i+=1
    return P_genes
    
