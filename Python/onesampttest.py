import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp

def one_sample_t_test(eset_T,eset_N,genes):
    eset_T = eset_T[genes]  
    eset_N = eset_N[genes]

    print(f"normal : {eset_N.shape}")
    print(f"tumor : {eset_T.shape}")
    print(f"genes : {genes.shape}")
    T_result = []
    for gene in genes:
        result = []
        for j in range(eset_T.shape[0]):
            tmp = ttest_1samp(eset_N[gene].values,eset_T[gene][j])
            result.append(tmp[0])
        T_result.append(result)
    T_result = pd.DataFrame(T_result, index = genes, columns = eset_T.index)
    return T_result