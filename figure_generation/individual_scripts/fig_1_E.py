import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
import numpy as np
from scipy import stats
import math

table = 'processed_data/benchmark_GM12878_results.txt'
df = pd.read_table(table)

#### filtering out negative controls with significance
df = df[~df['S'].str.contains('random') | (df['S'].str.contains('random') & (df['padj_S_run1'] > 0.001))]
df = df.dropna(subset=['log2FC_ES', 'log2FC_SE'])

annoA = {'2.REST','3.CTCF_TAD','4.CTCF_27ac', '5.GFI1', '6.YY1','7.active'}

## add annotation of Silencers
def dscri(x):
    if 'CTCF_TAD' in x:
        return '3.CTCF_TAD'
    elif 'CTCF_chr' in  x:
        return '4.CTCF_27ac'
    elif 'pos' in x:
        return '7.active'
    elif 'random' in x:
        return '1.random'
    elif 'YY1' in x:
        return '6.YY1'
    elif 'REST' in x:
        return '2.REST'
    elif 'GFI1' in x:
        return '5.GFI1'
    else:
        return '8.Blocker'

disc = lambda x: dscri(x)
df['annotation_S']= df['S'].apply(disc)
df = df[~(df['annotation_S']=='8.Blocker')]

df = df.sort_values(['annotation_S'], ascending=True)
low = df[df['E'] == 'En02']
high = df[df['E'] == 'En19']


for library in ['ES', 'SE']:
    key = 'log2FC_%s' % library

    ### calc median of background controls
    rndmL = low[(low['annotation_S']=='1.random')]
    rndmL_med = pd.Series(rndmL[key]).median()

    rndmH = high[(high['annotation_S']=='1.random')]
    rndmH_med = pd.Series(rndmH[key]).median()

    ## plot
    fig = plt.figure(figsize=(3,2.5), dpi=300)

    ymin, ymax = -4.5, 8.5
    pl = sns.light_palette("seagreen", n_colors=8)

    ax = fig.add_subplot(121, yticks=[-2,0,2,4,6])
    ax.set_aspect(1.4)
    fig.subplots_adjust(left=0.2, bottom=0.2)
    plt.ylim([-2,7])
    plt.hlines([rndmL_med], ymin, ymax, "coral", linewidth=1, zorder=0,alpha=0.5)
    sns.boxplot(x=low['annotation_S'],y=low[key], zorder=2, fliersize=0,palette=pl)
    plt.ylabel('expression log2 (fold change)', color='0.2',fontsize='xx-small')
    plt.xlabel('', color='0.2',fontsize='xx-small')
    plt.xticks(rotation=40,ha='right',fontsize="xx-small")
    plt.yticks(fontsize='xx-small')
    plt.text(0.5,6,'En02', color='0.2',fontsize='xx-small')

    k='0.4'
    ax.spines['right'].set_color(k)
    ax.spines['left'].set_color(k)
    ax.spines['top'].set_color(k)
    ax.spines['bottom'].set_color(k)
    ax.tick_params(axis='x', colors=k)
    ax.tick_params(axis='y', colors=k)


    ax2 = fig.add_subplot(122, yticks=[-2,0,2,4,6])
    ax2.set_aspect(1.4)
    fig.subplots_adjust(left=0.2, bottom=0.2)
    plt.ylim([-2,7])
    plt.hlines([rndmH_med], ymin, ymax, "coral",linewidth=1, zorder=0,alpha=0.5)
    sns.boxplot(x=high['annotation_S'],y=high[key], zorder=2, fliersize=0,palette=pl)
    plt.ylabel('', color='0.2',fontsize='xx-small')
    plt.xlabel('', color='0.2',fontsize='xx-small')
    plt.xticks(rotation=40,ha='right',fontsize="xx-small")
    plt.yticks(fontsize='xx-small')
    plt.text(0.5,6,'En19', color='0.2',fontsize='xx-small')


    k='0.4'
    ax2.spines['right'].set_color(k)
    ax2.spines['left'].set_color(k)
    ax2.spines['top'].set_color(k)
    ax2.spines['bottom'].set_color(k)
    ax2.tick_params(axis='x', colors=k)
    ax2.tick_params(axis='y', colors=k)

    fig.savefig('Fig1E.%s.box.pdf' % library, transparent=True)
