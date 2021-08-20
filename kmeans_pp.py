import numpy as np
import pandas as pd
import sys

import spkmeans
np.random.seed(0)

def init_data(src: str):
    datapoints = pd.read_csv(src, index_col=0, header=None).to_numpy()
    return spkmeans.prepareData(datapoints, datapoints.shape[0], datapoints.shape[1]) # datapoints, count, dim


def algo(data, k):
    centeroids = [np.random.choice(data.shape[0])]
    for _ in range(k-1):
        probs = np.zeros(shape=data.shape[0])
        for i, point in enumerate(data):
            probs[i] = ((np.array([data[x] for x in centeroids])-point) ** 2).sum(axis=1).min()
        centeroids.append(np.random.choice(data.shape[0], p=probs/np.sum(probs)))
    return centeroids

if __name__ == "__main__":
    assert len(sys.argv) >= 4

    K = int(sys.argv[1])
    src = sys.argv[-2]

    data = init_data(src)
    assert K <= data.shape[0]
    centers = algo(data,K)

    res = mykmeanspp.fit([data[i] for i in centers], data, MAX_ITER, data.shape[0], data.shape[1] ,K)
    res = np.round(res,decimals=4)
    print(str(centers)[1:-1].replace(" ",""))
    print("\n".join(",".join(f"{coord}" for coord in obs) for obs in res), end="")