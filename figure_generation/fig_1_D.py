### plotting duo and log additive model

import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from PIL import Image
import math
import numpy as np
import scipy
from sklearn import linear_model as lm

table = 'processed_data/benchmark_GM12878_results.txt'
df = pd.read_table(table)

#### filtering out negative controls with significance
df = df[~df['S'].str.contains('random') | (df['S'].str.contains('random') & (df['padj_S_run1'] > 0.001))]
df = df[~df['S'].str.contains('random') | (df['S'].str.contains('random') & (df['padj_S_run2'] > 0.001))]

df = df.dropna(subset=['log2FC_ES', 'log2FC_SE'])

# making single x single term
df['ES_ExS'] = df['log2FC_S_run2']*df['log2FC_E_run2']
df['SE_ExS'] = df['log2FC_S_run1']*df['log2FC_E_run1']

# linear regression
run2 = df.loc[:,['log2FC_S_run1', 'log2FC_E_run1','SE_ExS']].as_matrix()
SE = df['log2FC_SE'].as_matrix()
run3 = df.loc[:,['log2FC_S_run2', 'log2FC_E_run2', 'ES_ExS']].as_matrix()
ES = df['log2FC_ES'].as_matrix()

linSE = lm.LinearRegression()
linES = lm.LinearRegression()
linSE.fit(run2,SE)
linES.fit(run3,ES)

SEscore = linSE.score(run2,SE)
ESscore = linES.score(run3,ES)
print linSE.coef_
print linES.coef_

## add a column of additive model to plot
df['AP_SE_fit'] = linSE.intercept_ + (df['log2FC_S_run1'] * linSE.coef_[0]) + (df['log2FC_E_run1']*linSE.coef_[1]) + (df['SE_ExS']*linSE.coef_[2])
df['AP_ES_fit'] = linES.intercept_ + (df['log2FC_S_run2'] * linES.coef_[0]) + (df['log2FC_E_run2']*linES.coef_[1]) + (df['ES_ExS']*linSE.coef_[2])


SEscore = 'r^2 = '+ '{:.3f}'.format(linSE.score(run2,SE))
ESscore = 'r^2 = '+ '{:.3f}'.format(linES.score(run3,ES))
SEfunc = 'SE = '+'{:.3f}'.format(linSE.intercept_)+' + '+'{:.3f}'.format(linSE.coef_[0])+'S + '+'{:.3f}'.format(linSE.coef_[1])+'E + '+'{:.3f}'.format(linSE.coef_[2])+'ExS'
ESfunc = 'ES = '+'{:.3f}'.format(linES.intercept_)+' + '+'{:.3f}'.format(linES.coef_[0])+'S + '+'{:.3f}'.format(linES.coef_[1])+'E + '+'{:.3f}'.format(linES.coef_[2])+'ExS'

for key in ['ES', 'SE']:
    x = 'AP_%s_fit' % key
    y = 'log2FC_%s' % key

    ## plot
    fig = plt.figure(figsize=(4,4), dpi=300)
    ax = fig.add_subplot(111, xticks=[-2,0,2,4,6],yticks=[-2,0,2,4,6])

    plt.xlim((-2.9,7.5))
    plt.ylim((-2.9,7.5))

    plt.scatter(y=df[y], x=df[x], s=1, alpha=0.1,color='0.3', zorder=1)
    plt.axhline(y=0, color='lightgrey', linestyle='dotted', zorder=0)
    plt.axvline(x=0, color='lightgrey', linestyle='dotted', zorder=0)

    ## 45 degree line
    x2 = np.arange(-10,10,0.1)
    y2 = x2
    plt.plot(x2,y2, color='coral', alpha=0.8, zorder=4, linestyle='dotted')

    plt.ylabel('Expression log2(fold change)', color='0.2')
    plt.xlabel('log additive model', color='0.2')

    k= '0.4'
    ax.spines['right'].set_color(k)
    ax.spines['left'].set_color(k)
    ax.spines['top'].set_color(k)
    ax.spines['bottom'].set_color(k)
    ax.tick_params(axis='x', colors=k)
    ax.tick_params(axis='y', colors=k)

    plt.savefig('Fig1D.linmodel_%s.png' % key)


with open("Fig1D.linmodel_score.txt", "w") as results:
    results.write('SEscore: %s\n' % SEscore)
    results.write('%s\n' % SEfunc)
    results.write('ESscore: %s\n' % ESscore)
    results.write('%s\n' % ESfunc)
