import pandas as pd
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
import numpy as np
from scipy import stats
import math


df = pd.read_table('norm_count.target1.txt')
targets = ['EPB41L1','AAR2', 'DLGAP4', 'MYL9', 'TGIF2', 'RAB5IF', 'SLA2']
df = df[df['Gene_name'].isin(targets)].reset_index(drop=True)

df1 = pd.DataFrame(index=[], columns=['Gene', 'norm_UMI_count', 'sample'])

for i in ['NT_r1', 'NT_r2', 'NT_r3', 'NT_r4', 'NT_r5']:
    k = df.loc[:,['Gene_name',i]]
    k['sample'] = 'ctrl'
    k.columns = df1.columns
    df1 = pd.concat([df1,k], axis=0)

for i in ['target1_r1', 'target1_r2', 'target1_r3', 'target1_r4', 'target1_r5']:
    k = df.loc[:,['Gene_name',i]]
    k['sample'] = 'KO'
    k.columns = df1.columns
    df1 = pd.concat([df1,k], axis=0)

df1 = df1.reset_index(drop=True)

for i in range(len(df1)):
    g = df1.loc[i,'Gene']
    if g == 'EPB41L1':
        h = 1
    elif g == 'AAR2':
        h = 2
    elif g == 'DLGAP4':
        h = 3
    elif g == 'MYL9':
        h = 4
    elif g == 'TGIF2':
        h = 5
    elif g == 'RAB5IF':
        h = 6
    else:
        h = 7
    df1.loc[i,'order'] = h

df1 = df1.sort_values(by='order')

fig = plt.figure(figsize=(2 ,1.75), dpi=300)
x = sns.boxplot(x="Gene", y="norm_UMI_count", hue='sample', data=df1, whis=[50, 50], showfliers=False, showbox=False, showcaps=False)
x = sns.swarmplot(x='Gene', y='norm_UMI_count', hue='sample', data=df1, dodge=True, size=3)
x.set(ylim=(0, 1100))
plt.xticks(rotation=15, fontsize='xx-small')
plt.yticks(fontsize='xx-small')
plt.xlabel('')
plt.ylabel('norm. UMI count', fontsize='x-small')
plt.legend('')

plt.title('target1')
x.figure.savefig('plots/fig_5_A_BRBseq_target1.pdf')
