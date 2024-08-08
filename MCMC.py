import numpy as np
import pandas as pd
import pymc as pm
import aesara.tensor as at

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
# Read the CSV files
data1_cd4 = pd.read_csv("small_dataset/Donor1_CD4_Genes.csv")
data1_cd8 = pd.read_csv("small_dataset/Donor1_CD8_Genes.csv")
data2_cd4 = pd.read_csv("small_dataset/Donor2_CD4_Genes.csv")
data2_cd8 = pd.read_csv("small_dataset/Donor2_CD8_Genes.csv")

# List of protein-coding genes
protein_coding_genes = ["MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", 
                        "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", 
                        "MT-ND4L", "MT-ND5", "MT-ND6", "MT-ATP6", 
                        "MT-ATP8"]

selected_data1_cd4 = data1_cd4[protein_coding_genes + ["non-MT"]]
selected_data1_cd8 = data1_cd8[protein_coding_genes + ["non-MT"]]
selected_data2_cd4 = data2_cd4[protein_coding_genes + ["non-MT"]]
selected_data2_cd8 = data2_cd8[protein_coding_genes + ["non-MT"]]

print(selected_data1_cd4.head())

# Load the data
df = selected_data1_cd4
x = df.iloc[:, -1].values  # x is the last column
J = df.shape[0]
V = 13
K = 2
# Define the data for the model
data_modelling = {
    'df': df,
    'x': x,
    'V': V,
    'J': J,
    'e_0': 0.01,
    'f_0': 0.01,
    'h': 0.01,
    'K': K
}

# Define the model
with pm.Model() as model:
    # Priors for the model parameters
    r = pm.Gamma('r', alpha=data_modelling['e_0'], beta=data_modelling['h'], shape=J)

    alpha = pm.Gamma('alpha', alpha=data_modelling['e_0'], beta=data_modelling['f_0'], shape=V)
    gamma = pm.Gamma('gamma', alpha=data_modelling['e_0'], beta=data_modelling['f_0'], shape=K)
    h = pm.Gamma('h', alpha=data_modelling['e_0'], beta=data_modelling['f_0'])

    # Initialize beta
    beta = pm.Normal('beta', mu=0, sigma=at.sqrt(1/alpha), shape=V)
    
    # Initialize Phi and theta
    Phi = pm.MvNormal('Phi', mu=at.zeros(V), cov=at.eye(V), shape=(V, V))
    theta = pm.Normal('theta', mu=0, sigma=at.sqrt(1/gamma), shape=(J, V))

    # Logistic function for p
    phi = pm.Deterministic('phi', at.dot(data_modelling['x'], beta.T) + at.dot(theta, Phi.T))
    p = pm.Deterministic('p', pm.math.sigmoid(phi))
    
    # Likelihood
    n = pm.NegativeBinomial('n', n=r, p=p, observed=df.values)

# Sampling from the model
with model:
    trace = pm.sample(draws=5000, tune=1000, chains=1, cores=1)

# Analyzing the trace
pm.summary(trace, var_names=['r', 'p', 'phi'])
