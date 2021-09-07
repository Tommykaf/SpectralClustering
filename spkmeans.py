import numpy as np
import pandas as pd
import sys

import spkmeans
np.random.seed(0)


def read_data(src: str):
    return pd.read_csv(src, header=None).to_numpy()


def pretty_print(matrix: np.ndarray):
    # res = np.round(matrix, decimals=4)
    format_zero = lambda x : 0 if 0 > x > -0.00005 else x
    print("\n".join(",".join(f"{format_zero(coord):.4f}" for coord in row) for row in matrix), end="")


def algo(data, k):
    centeroids = [np.random.choice(data.shape[0])]
    for _ in range(k-1):
        probs = np.zeros(shape=data.shape[0])
        for i, point in enumerate(data):
            probs[i] = ((np.array([data[x] for x in centeroids])-point) ** 2).sum(axis=1).min()
        centeroids.append(np.random.choice(data.shape[0], p=probs/np.sum(probs)))
    return centeroids


if __name__ == "__main__":
    #spkmeans K GOAL FILENAME
    assert len(sys.argv) == 4

    K = int(sys.argv[1])
    assert K >= 0
    goal = sys.argv[2].lower()
    assert goal in {"spk", "ddg", "wam", "lnorm", "jacobi"}
    src = sys.argv[3]

    data_matrix = read_data(src)
    count = data_matrix.shape[0]
    dim = data_matrix.shape[1]
    assert len(data_matrix) > K

    # switch(goal)
    if goal == 'jacobi':
        eigen_array, V = spkmeans.jacobi(data_matrix, count)
        print(eigen_array)
        pretty_print(V)
    elif goal == 'wam':
        pretty_print(spkmeans.WAM(data_matrix, count, dim))
    elif goal == 'ddg':
        pretty_print(spkmeans.DDG(data_matrix, count, dim))
    elif goal == 'lnorm':
        pretty_print(spkmeans.lNorm(data_matrix, count, dim))
    elif goal == 'spk':
        t_array = spkmeans.prepareData(data_matrix, count, dim, K)
        K = t_array.shape[1]
        centers_indices = algo(t_array, K)
        print(str(centers_indices)[1:-1].replace(" ","")) # [1, 2, 3, 4] => 1,2,3,4
        pretty_print(spkmeans.fit([t_array[i] for i in centers_indices], t_array, count, dim, K))  