# MPD

MPD is a program that identifies individual driver genes.

The details of the method are described in (url)

### Experiment environment
python = 3.6.2

### Library
+ numpy
+ pandas
+ pickle
+ sklearn
+ scipy
+ networkx
+ time

### Experiment example
The code below is the 'main.py' in the Python folder and you can proceed as described below when experimenting with your dataset.

```python
# Path where the dataset exists
DATA_PATH = ""
# Path to save Effect Matrix
EM_PATH = ""
# Cancer type
c_type = ""

# Expression files with only tumor samples
eset_T = pd.read_csv(f"{DATA_PATH}{c_type}_eset_T.csv",index_col=0) 
# Expression files with only normal samples
eset_N = pd.read_csv(f"{DATA_PATH}{c_type}_eset_N.csv",index_col=0) 
# File expressing binary with or without mutataion
mset = pd.read_csv(f"{DATA_PATH}{c_type}_mset_T.csv",index_col=0) 
# File expressing binary with or without CNV
CNV = pd.read_csv(f"{DATA_PATH}{c_type}_CNV.csv",index_col='Patient_barcode') 
# File expressing metyhlation information
methyl = pd.read_csv(f"{DATA_PATH}{c_type}_methylation.csv",index_col=0)  
# Network file
edges = pd.read_csv(f"{DATA_PATH}network.csv")
# File expressing known driver gene
KO = pd.read_csv(f"{DATA_PATH}/{c_type}_driver.csv")['Gene_symbol'].values

# Data preprocessing process
eset_T,eset_N,SM,edges,genes = prepare.prepare(eset_T, eset_N, mset, CNV, methyl, edges)
samples = eset_T.index.values
# 1-sample t-test
T_result = onesampttest.one_sample_t_test(eset_T,eset_N,genes)
# The process of determining the weight of self-loop by the value of t, which is the result of the t-test
DEGscore = DEGscore.DEGscore(genes,T_result)
# The process of classifying genes to be used for train and test
train_genes, test_genes = fold_genes(genes, KO)
# The process of calculating and storing the effect matrix
Effect_Matrix.EM(EM_PATH, eset_T, eset_N, edges, genes, DEGscore)
# Process of identifying driver gene by sample by train and test model
train_P, train_N = train_model_test.make_train_data(EM_PATH, train_genes, SM, genes, samples)
clfs = train_model_test.learning (train_P,train_N)
drivers = train_model_test.get_driver(genes,test_genes,clfs,samples)
```
