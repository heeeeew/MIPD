import pandas as pd
from sklearn.preprocessing import MinMaxScaler

def DEGscore(gene,T_result):
    T_result = T_result.T
    scaler = MinMaxScaler(feature_range=(0.1,0.9))
    score = []
    for i in range(T_result.shape[0]):
        D = np.abs(T_result.iloc[i].values)
        D = D.reshape(-1,1)
        D = scaler.fit_transform(D)
        score.append(D.reshape(-1))

    DEGscore = pd.DataFrame(score, index = T_result.index,columns = genes)
    return DEGscore