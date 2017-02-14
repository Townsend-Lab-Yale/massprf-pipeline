import pandas as pd
import numpy as np

f = open('/Users/jaystanley/Documents/School/Grad/townsend/ricemassprf/junruiricemap.csv')
df = pd.read_csv(f)

for x in df:
    print(list(df[x][(pd.notnull(df[x]))]))