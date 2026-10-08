import numpy as np
import csv
import random
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import seaborn as sns
import pandas as pd
df = pd.read_csv("../Data/GSE150910_gene-level_count_file.csv", index_col=0)
print(df.head())