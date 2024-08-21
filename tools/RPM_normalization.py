import pandas as pd
from bioinfokit.analys import norm
gene_length = pd.read_csv("mt-datasets/gene_length.csv")


def RPM(df):
    df = df[gene_length.columns]

    # Divide the read counts by the length of each gene
    # df_norm_1 = df.div(gene_length.iloc[0])
    
    #df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
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
