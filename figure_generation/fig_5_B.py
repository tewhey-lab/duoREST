import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from PIL import Image
import math
import numpy as np
import seaborn as sns

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
ref_strand = ref['strand'].to_dict()

alt = pd.read_table('processed_data/variants_alt_fimo.txt', sep='\t')
alt = alt[(alt['start']== 91) | (alt['start']== 92)]
alt = alt.sort_values(by='#pattern name', ascending=False)
alt = alt.drop_duplicates(subset=['sequence name'])
spl_a = alt['sequence name'].str.split(':A:wP', expand=True)
spl_a.columns = ['SNP', 'position']
alt = pd.concat([alt, spl_a['SNP']],axis=1).set_index('SNP')
alt_score = alt['score'].to_dict()
alt_strand = alt['strand'].to_dict()


df['score_maj'] = 0
df['score_min'] = 0
df['strand'] = 'NA'

df = df.reset_index(drop=True)
for i in range(len(df)):
    id = df.loc[i,'ID']
    s = df.loc[i,'SNP']
    if 'R:wP' in id:
        df.loc[i,'score_maj'] = ref_score.get(s,0)
        df.loc[i,'score_min'] = alt_score.get(s,0)
        df.loc[i,'strand'] = ref_strand.get(s,0)


    elif 'A:wP' in id:
        df.loc[i,'score_maj'] = alt_score.get(s,0)
        df.loc[i,'score_min'] = ref_score.get(s,0)
        df.loc[i,'strand'] = alt_strand.get(s,0)


df = df[(df['score_maj'] > 0) & (df['score_min'] > 0)]
df['dBS'] = df['score_min'] - df['score_maj']

### orient motifs in a same direction
df['rel_pos'] = 0
df = df.reset_index(drop=True)
for i in range(len(df)):
    pos = df.loc[i,'SNP_pos']
    str = df.loc[i,'strand']
    if str == '+':
        df.loc[i,'rel_pos'] = pos
    elif str == '-':
        df.loc[i,'rel_pos'] = 202 - pos


### taking absolute value of logSkew
df['abs_skew'] = df['LogSkew'].abs()

df['dist'] = abs(101 - df['rel_pos'])
df = df.sort_values(by='dist')
df = df.drop_duplicates(subset=['SNP'])

df = df[(df['rel_pos'] <= 126) & (df['rel_pos'] >= 76)]


### separate above and below for line plots
above = df[(df['score_min'] > seg)| (df['score_maj'] > seg)]
below = df[(df['score_min'] < seg) & (df['score_maj'] < seg)]



### scatter: position-vs-skew
fig2 = plt.figure(figsize=(3,3), dpi=300)
boxprops = dict(linewidth=0.5, color='0.6')

ax2 = plt.subplot(211)
ax2.set(ylim=(-0.1, 3.5))
ax2.set(xlim=(76, 126))
ax2 = sns.boxplot(y=above['LogSkew'].abs(), x=above['rel_pos'],fliersize=0,color='coral',boxprops=dict(linewidth=1))
ax2 = plt.ylabel('abs_LogSkew_above', color='0.2', fontsize=6)
ax2 = plt.xlabel('')
plt.xticks(rotation=90,fontsize=4)
plt.yticks(fontsize=5)


ax3 = plt.subplot(212)
ax3.set(ylim=(-0.1, 3.5))
ax3.set(xlim=(76, 126))
ax3 = sns.boxplot(y=below['LogSkew'].abs(), x=below['rel_pos'],fliersize=0,color="darkgrey")
ax3 = plt.ylabel('abs_LogSkew_below', color='0.2',fontsize=6)
ax3 = plt.xlabel('SNP position (relative)', color='0.2', fontsize=4)
plt.xticks(rotation=90, fontsize=4)
plt.yticks(fontsize=5)

fig2 =  plt.savefig('Fig5b_abs_skew_box.pdf')
