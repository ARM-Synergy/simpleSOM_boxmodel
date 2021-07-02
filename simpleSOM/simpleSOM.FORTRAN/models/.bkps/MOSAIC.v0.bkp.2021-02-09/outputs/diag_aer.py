
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt

df = pd.read_csv('diag.aer',header=None,delim_whitespace=True)

df1 = df.iloc[:,:13]

df1 = df.iloc[:,13:]



print(df1)
