import numpy as np
import tools as tl
from copulas.multivariate import GaussianMultivariate


def get_cor_copulas(data, gene_list):
    cor_o = np.zeros((13, 13))
    copula = GaussianMultivariate()
    for i in range(13):
        for j in range(i+1, 13):
            df_0 = data[[gene_list[i], gene_list[j]]]
            df = tl.filter_outliers(df_0)
            copula.fit(df)
            synthetic = copula.sample(len(df))
            cor_o[i, j] = synthetic.corr().iloc[0, 1]
    return (cor_o)
