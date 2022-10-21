import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from PIL import Image
import math
import numpy as np
from scipy import stats


seg = 20.86
enhancer = 'En19'

df = pd.read_table('processed_data/RESTscreen_AF_K562_emVAR.out', sep='\t')
df = df.dropna(subset=['LogSkew'])
df = df[df['ID'].str.contains(enhancer)]
df = df[df['ID'].str.contains(':wP')]


### making a column for SNP position
spl = df['ID'].str.split(':wP',expand=True)
spl.columns = ['sil','SNP_pos']
df = pd.concat([df,spl['SNP_pos']],axis=1)
df['SNP_pos'] = df['SNP_pos'].astype(int)


df = df.reset_index(drop=True)

### calc delta binding score
ref = pd.read_table('processed_data/variants_ref_fimo.txt', sep='\t')
ref = ref[(ref['start']== 91) | (ref['start']== 92)]
ref = ref.sort_values(by='#pattern name', ascending=False) # prioritize MA0138.2
ref = ref.drop_duplicates(subset=['sequence name'])
spl_r = ref['sequence name'].str.split(':R:wP', expand=True)
spl_r.columns = ['SNP', 'position']
ref = pd.concat([ref, spl_r['SNP']],axis=1).set_index('SNP')
ref_score = ref['score'].to_dict()

alt = pd.read_table('processed_data/variants_alt_fimo.txt', sep='\t')
alt = alt[(alt['start']== 91) | (alt['start']== 92)]
alt = alt.sort_values(by='#pattern name', ascending=False)
alt = alt.drop_duplicates(subset=['sequence name'])
spl_a = alt['sequence name'].str.split(':A:wP', expand=True)
spl_a.columns = ['SNP', 'position']
alt = pd.concat([alt, spl_a['SNP']],axis=1).set_index('SNP')
alt_score = alt['score'].to_dict()


df['score_maj'] = 0
df['score_min'] = 0
df = df.reset_index(drop=True)
for i in range(len(df)):
    id = df.loc[i,'ID']
    s = df.loc[i,'SNP']
    if 'R:wP' in id:
        df.loc[i,'score_maj'] = ref_score.get(s,0)
        df.loc[i,'score_min'] = alt_score.get(s,0)
    elif 'A:wP' in id:
        df.loc[i,'score_maj'] = alt_score.get(s,0)
        df.loc[i,'score_min'] = ref_score.get(s,0)


df = df[(df['score_maj'] > 0) & (df['score_min'] > 0)]
df['dBS'] = df['score_min'] - df['score_maj']



### select one result for one SNP (select the closest one to the center)
df['distance'] = abs(101 - df['SNP_pos'])
df = df.sort_values(by = 'distance')
df = df.drop_duplicates(subset=['SNP'])


### select "inside" SNP
df = df[(df['SNP_pos']>=91) & (df['SNP_pos']<=111)]

df =df.dropna(subset=['LogSkew'])


### remove outliers of dBS derived from mis score calc
df = df[(df['dBS'] > -20) & (df['dBS'] < 20)]

### plots
for segment in ['above', 'below']:

    if segment == 'above':
        df1 = df[(df['score_maj'] > seg)| (df['score_min'] > seg)]
    elif segment == 'below':
        df1 = df[(df['score_maj'] < seg) & (df['score_min'] < seg)]

    skew = np.array(df1['LogSkew'].tolist()).reshape(-1,1)
    dbs = np.array(df1['dBS'].tolist()).reshape(-1,1)
    corr, pval = stats.pearsonr(dbs, skew)

    fig = plt.figure(figsize=(2,2), dpi=300)
    ax = fig.add_subplot(111)
    ax.set(ylim=(-4, 4))
    ax.set(xlim=(-11, 11))
    fig = plt.ylabel('LogSkew (Alt-Ref)', color='0.2', fontsize='xx-small')
    fig = plt.xlabel('delta BS (Alt-Ref)', color='0.2',fontsize='xx-small')
    fig =  plt.scatter(y=df1['LogSkew'], x=df1['dBS'], s=2, alpha=0.6,color="0.3", zorder=1)
    fig = plt.text(-10,-3,corr, color='0.2',fontsize='xx-small')

    fig = plt.axhline(y=0, color='0.8', zorder=0)
    fig = plt.axvline(x=0, color='0.8', zorder=0)

    k='0.4'
    ax.spines['right'].set_color(k)
    ax.spines['left'].set_color(k)
    ax.spines['top'].set_color(k)
    ax.spines['bottom'].set_color(k)
    ax.tick_params(axis='x', colors=k)
    ax.tick_params(axis='y', colors=k)

    fig =  plt.savefig('Fig5C_score_skew.%s.pdf' % segment)
