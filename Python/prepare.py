import numpy as np
import pandas as pd

def prepare(eset_T, eset_N, mset, CNV, methyl, edges):
    T_samples = eset_T.columns.values
    M_samples = mset.columns.values
    C_samples = CNV.index.values
    L_samples = methyl.index.values

    print(f"L : {L_samples.shape}")
    print(f"T : {T_samples.shape}")
    print(f"C : {C_samples.shape}")
    print(f"M : {M_samples.shape}")
    print(f"N : {eset_N.shape}")

    tmp = []
    for samp in T_samples:
        tmp.append(samp[:12])
    T_samples = tmp
    eset_T = eset_T.T
    eset_N = eset_N.T
    Aliquot = eset_T.index.values
    eset_T['Aliquot_barcode'] = Aliquot
    eset_T['Patient_barcode'] = T_samples
    eset_T = eset_T.set_index('Patient_barcode')

    tmp = []
    for samp in L_samples:
        tmp.append(samp[:12])
    L_samples = tmp
    methyl.index = L_samples

    samples = sorted(list(set(T_samples) & set(M_samples) & set(C_samples) & set(L_samples)))

    print("get samples")

    eset_T = eset_T.loc[samples]
    eset_T = eset_T.set_index('Aliquot_barcode')

    T_genes = eset_T.columns.values

    tmp_T = eset_T.values
    tmp_T = tmp_T.astype(np.bool)
    tmp_T = tmp_T.sum(axis=0)
    T_cut = int(eset_T.shape[0]*0.2)
    T_genes = set(T_genes[tmp_T>T_cut])
    edges_genes = set(edges['source']) | set(edges['destination'])
    SM_genes = set(mset.index) | set(CNV.columns) | set(methyl.columns)
    genes = sorted(list(T_genes & edges_genes & SM_genes))

    print("get genes")

    tmp = []
    for src,dest in edges.values:
        if src in genes and dest in genes:
            tmp.append([src,dest])
            
    edges = pd.DataFrame(tmp,columns=['src','dest'])
    print("get edges")

    eset_T = eset_T[genes]
    eset_N = eset_N[genes]

    mset = mset.T.loc[samples]
    CNV = CNV.loc[samples]
    methyl = methyl.loc[samples]
    tmp = methyl.values.astype(int)
    methyl = pd.DataFrame(tmp,index=methyl.index,columns=methyl.columns)
    CNV = CNV.drop('Aliquot_barcode',axis=1)
    c_tmp = []
    cut = CNV.shape[1]
    for samp in samples:
        c_row = np.array(CNV.loc[samp])
        if c_row.shape[0] != cut:
            tmp = c_row[0]
            for row in c_row:
                tmp = tmp|row
            c_row = tmp
        c_tmp.append(c_row)
    CNV = pd.DataFrame(c_tmp,index=samples,columns=CNV.columns)

    mset = pd.DataFrame(mset.values.astype(int),index=mset.index,columns=mset.columns)

    tmp_L = pd.DataFrame(np.zeros([len(samples),len(genes)]).astype(int),index=samples,columns=genes)
    tmp_M = pd.DataFrame(np.zeros([len(samples),len(genes)]).astype(int),index=samples,columns=genes)
    tmp_C = pd.DataFrame(np.zeros([len(samples),len(genes)]).astype(int),index=samples,columns=genes)

    L_genes = sorted(list(set(methyl.columns) & set(genes)))
    M_genes = sorted(list(set(mset.columns) & set(genes)))
    C_genes = sorted(list(set(CNV.columns) & set(genes)))

    tmp_L[L_genes] = methyl[L_genes]
    tmp_M[M_genes] = mset[M_genes]
    tmp_C[C_genes] = CNV[C_genes]

    SM = tmp_L | tmp_M | tmp_C


    print("make SM")
    genes = pd.DataFrame(genes,columns=['GeneSymbol'])

    print(f"normal : {eset_N.shape}")
    print(f"tumor : {eset_T.shape}")
    print(f"genes : {len(genes)}")
    print(f"edges : {edges.shape}")
    
    return(eset_T,eset_N,SM,edges,genes)