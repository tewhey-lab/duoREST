### plotting Fig 1C
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from PIL import Image
import math
import numpy as np
import scipy

table = 'processed_data/benchmark_GM12878_results.txt'
df = pd.read_table(table)

#### filtering out negative controls with significance
df = df[~df['S'].str.contains('random') | (df['S'].str.contains('random') & (df['padj_S_run1'] > 0.001))]
df = df[~df['S'].str.contains('random') | (df['S'].str.contains('random') & (df['padj_S_run2'] > 0.001))]

## calc Pearson's r
df = df.drop_duplicates(subset=['log2FC_ES'])
df1 = df.dropna(subset=['log2FC_ES', 'log2FC_SE'])
cor = df1.iloc[:,[8,14]]
cor_p = cor.corr(method='pearson')
cor_p = cor_p.iloc[0,1]

r1 = 'r = '+ '{:.2f}'.format(cor_p)

## plot
fig = plt.figure(figsize=(4,4), dpi=300)
ax = fig.add_subplot(111, xticks=[-2,0,2,4,6],yticks=[-2,0,2,4,6])
plt.xlim((-2.9,7.5))
plt.ylim((-2.9,7.5))
plt.axhline(y=0, color='lightgrey', linestyle='dotted', zorder=0)
plt.axvline(x=0, color='lightgrey', linestyle='dotted', zorder=0)

## scatter plot
plt.scatter(y=df['log2FC_ES'], x=df['log2FC_SE'], s=1, alpha=0.1,color='0.3', zorder=1)

## 45 degree line
x2 = np.arange(-10,10,0.1)
y2 = x2
plt.plot(x2,y2, color='coral', alpha=0.8, zorder=4, linestyle='dotted')
plt.text(-2.3,6.5,r1, color='0.2')

plt.ylabel('Expression log2 (Fold Change)', color='0.2')
plt.xlabel('Expression log2 (Fold Change)', color='0.2')

k= '0.4'
ax.spines['right'].set_color(k)
ax.spines['left'].set_color(k)
ax.spines['top'].set_color(k)
ax.spines['bottom'].set_color(k)
ax.tick_params(axis='x', colors=k)
ax.tick_params(axis='y', colors=k)

plt.savefig('Fig1C.png')
