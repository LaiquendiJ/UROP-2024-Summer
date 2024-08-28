import numpy as np
import pandas as pd
gene_length = pd.read_csv("mt-datasets/gene_length.csv")


def RPM(df):
    df = df[gene_length.columns]

    # Divide the read counts by the length of each gene
    # df_norm_1 = df.div(gene_length.iloc[0])

    # df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    # Drop non-MT values
    df_norm_1 = df.drop(columns=['non-MT'])

    row_sums = df_norm_1.sum(axis=1)

    # Drop zero lines to aviod dividing something by zero
    zero_indices = row_sums.index[row_sums == 0].tolist()
    df_norm_1 = df_norm_1.drop(index=zero_indices)

    row_sums = row_sums.drop(index=zero_indices)
    df_normalized = df_norm_1.div(row_sums, axis=0)

    # Get the interested mtRNA columns
    interested_genes = ["MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
                        "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4",
                        "MT-ND4L", "MT-ND5", "MT-ND6", "MT-ATP6",
                        "MT-ATP8"]
    # df_normalized_fin = round(df_normalized[interested_genes] *1000)
    df_normalized_fin = round(df_normalized[interested_genes] * 1000000)
    return (df_normalized_fin)


def filter_outliers(df, l_range=3, u_range=3):
    df_i = df.values
    data1 = df_i[:, 0]  # Access the first column
    data2 = df_i[:, 1]  # Access the second column

    # Calculate the IQR for data1
    Q1_1 = np.percentile(data1, 25)
    Q3_1 = np.percentile(data1, 75)
    IQR_1 = Q3_1 - Q1_1
    bol_1 = (data1 < (Q3_1 + u_range * IQR_1)
             ) & (data1 > (Q3_1 - l_range * IQR_1))

    # Calculate the IQR for data2
    Q1_2 = np.percentile(data2, 25)
    Q3_2 = np.percentile(data2, 75)
    IQR_2 = Q3_2 - Q1_2
    bol_2 = (data2 < (Q3_2 + u_range * IQR_2)
             ) & (data1 > (Q3_2 - l_range * IQR_2))

    # Combine the boolean conditions
    bol = bol_1 & bol_2

    # Return the filtered DataFrame
    result = pd.DataFrame({'x': data1[bol], 'y': data2[bol]})
    result.columns = df.columns

    return result
