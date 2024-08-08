import numpy as np
import pandas as pd


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
