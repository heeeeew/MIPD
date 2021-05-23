# MPD

MPD는 개인별 드라이버 유전자를 식별하기 위한 프로그램입니다.

이 방법에 대한 자세한 내용은 링크를 참조하세요 (url)

### 실험 환경
python = 3.6.2

### 사용 라이브러리
+ numpy
+ pandas
+ pickle
+ sklearn
+ scipy
+ networkx
+ time

### 실험 예제
아래의 코드는 Python 폴더 내의 main.py이며 당신의 데이터셋으로 실험을 진행할때 아래의 설명대로 진행하면 된다.

```
# 데이터셋이 존재하는 경로
DATA_PATH = ""
# Effect Matrix를 저장할 경로
EM_PATH = ""
# 암 종류
c_type = ""

# tumor sample만 존재하는 expression 파일
eset_T = pd.read_csv(f"{DATA_PATH}{c_type}_eset_T.csv",index_col=0) 
# normal sample만 존재하는 expression 파일
eset_N = pd.read_csv(f"{DATA_PATH}{c_type}_eset_N.csv",index_col=0) 
# mutataion 유무가 binary를 표현된 파일
mset = pd.read_csv(f"{DATA_PATH}{c_type}_mset_T.csv",index_col=0) 
# CNV 유무가 binary로 표현된 파일 
CNV = pd.read_csv(f"{DATA_PATH}{c_type}_CNV.csv",index_col='Patient_barcode') 
# metyhlation 정보
methyl = pd.read_csv(f"{DATA_PATH}{c_type}_methylation.csv",index_col=0)  
# network file
edges = pd.read_csv(f"{DATA_PATH}network.csv")
# 알려진 driver gene이 암종별로 표현된 파일
KO = pd.read_csv(f"{DATA_PATH}/{c_type}_driver.csv")['Gene_symbol'].values

# 데이터 전처리과정
eset_T,eset_N,SM,edges,genes = prepare.prepare(eset_T, eset_N, mset, CNV, methyl, edges)
samples = eset_T.index.values
# 1-sample t-test
T_result = onesampttest.one_sample_t_test(eset_T,eset_N,genes)
# t-test 결과로 나온 t값으로 self-loop의 가중치를 결정하는 과정
DEGscore = DEGscore.DEGscore(genes,T_result)
# train 및 test에 사용할 유전자를 구분하는 과정
train_genes, test_genes = fold_genes(genes, KO)
# effect matrix를 계산하여 저장하는 과정
Effect_Matrix.EM(EM_PATH, eset_T, eset_N, edges, genes, DEGscore)
# model을 train 및 test하여 sample별 driver gene을 식별하는 과정
train_P, train_N = train_model_test.make_train_data(EM_PATH, train_genes, SM, genes, samples)
clfs = train_model_test.learning (train_P,train_N)
drivers = train_model_test.get_driver(genes,test_genes,clfs,samples)
```
