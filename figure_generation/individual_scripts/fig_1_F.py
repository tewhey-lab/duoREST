import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster import hierarchy
from scipy import stats

table = 'processed_data/benchmark_GM12878_results.txt'
df = pd.read_table(table)

#### filtering out negative controls with significance
df = df[~df['S'].str.contains('random') | (df['S'].str.contains('random') & (df['padj_S_run1'] > 0.001))]
df = df.dropna(subset=['log2FC_ES', 'log2FC_SE'])

#annoA = {'1.CTCF_nonTAD','2.CTCF_TAD','3.GFI1','4.REST', '5.YY1', '6.active'}

### P's rank is based on the P's activity on run2
dfp = df.sort_values(by='E',ascending=True)
dfp = pd.Series(dfp['E'])
dfp = dfp.drop_duplicates()
lp = dfp.values.tolist()

### make an empty dataframe
output = pd.DataFrame(columns=['S', 'E', 'log2FC_E_run1', 'log2FC_E_run2','zES', 'zSE'])

### calc z scores

for p in lp:
    df1 = df[df['E'] == p].reset_index(drop=True)
    rndm = df1[df1['annotation_S'] == '7.random_genomic']
    rndm = rndm.dropna(subset=['log2FC_ES'])

    if len(rndm) >0:
        df1['zES'] = (df1['log2FC_ES'] - rndm['log2FC_ES'].mean()) / rndm['log2FC_ES'].std(ddof=0)
        df1['zSE'] = (df1['log2FC_SE'] - rndm['log2FC_SE'].mean()) / rndm['log2FC_SE'].std(ddof=0)

        df1 = df1.loc[:,['S','E','log2FC_E_run1', 'log2FC_E_run2', 'zES', 'zSE']]

        output = pd.concat([output, df1], axis=0)

    else:
        pass


key = 'REST'

output = output[output['S'].str.contains(key)]
output = output.dropna(subset=['zES'])

xmin, xmax = -20,20
fig = plt.figure(figsize=(16,16))

plt.subplot(2,1,1)
sns.boxplot(data = output, x='E', y='zES', color="darkgrey",fliersize=0)
plt.xticks(rotation=30)
plt.ylim(-4.5,4.5)
plt.hlines([0],xmin, xmax, "0.7", linestyles='dashed', zorder=0)
plt.ylabel('zscore ES')

plt.subplot(2,1,2)
sns.boxplot(data = output, x='E', y='zSE', color="darkgrey",fliersize=0)
plt.xticks(rotation=30)
plt.ylim(-5,5)
plt.hlines([0],xmin, xmax, "0.7", linestyles='dashed', zorder=0)
plt.ylabel('zscore SE')

plt.savefig('Fig1F_box.REST.png')
