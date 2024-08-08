import numpy as np
import pandas as pd
import pymc as pm

import statsmodels.api as sm
import plotly.express as px
import plotly.figure_factory as ff
from arch import arch_model

from ipywidgets import HBox, VBox
from scipy.optimize import fmin, minimize
from scipy.stats import t
from scipy.stats import norm
from math import inf

import bs4 as bs
import requests
import yfinance as yf
import datetime

from scipy.stats import gamma
from scipy.stats import norm
from scipy.stats import t
from scipy.stats import beta
import matplotlib.pyplot as plt

from statsmodels.distributions.empirical_distribution import ECDF

from copulas.visualization import scatter_2d
from copulas.visualization import dist_1d
from copulas.multivariate import GaussianMultivariate

import tools as tl

# Read the CSV files
data1_cd4 = pd.read_csv("small_dataset/Donor1_CD4_Genes.csv")
df_0 = data1_cd4[["MT-CO1", "MT-CO2"]]

normalized_gibbs = pd.read_csv("gibbs_results.csv")
observed = df_0.values
observed_Y = np.array(df_0[["MT-CO2"]])
name_1 = data1_cd4.columns.values[1]

protein_coding_genes = ["MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", 
                        "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", 
                        "MT-ND4L", "MT-ND5", "MT-ND6", "MT-ATP6", 
                        "MT-ATP8"]
normalized_gibbs.columns=protein_coding_genes
df_gib = normalized_gibbs[["MT-CO1", "MT-CO2"]]
df_gib = tl.filter_outliers(df_gib)
name_1 = protein_coding_genes[0]
name_2 = protein_coding_genes[1]
df = tl.filter_outliers(df_0)
df_gib = tl.filter_outliers(df_gib)
scatter_2d(df)