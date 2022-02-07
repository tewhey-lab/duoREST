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

df = pd.read_table('processed_data/RESTscreen_derived_K562_emVAR.out', sep='\t')
df = df.dropna(subset=['LogSkew'])
df = df[df['ID'].str.contains(enhancer)]

### making a column for SNP position
spl = df['silencer'].str.split(':wP',expand=True)
spl.columns = ['sil','SNP_pos']
df = pd.concat([df,spl['SNP_pos']],axis=1)
df['SNP_pos'] = df['SNP_pos'].astype(int)



#df = df.dropna(subset=['Ref_BS'])


df = df.reset_index(drop=True)

### calc delta binding score
ref = pd.read_table('processed_data/variants_ref_fimo.txt', sep='\t')
spl_r = df['sequence name'].str.split(':R:wP', expand=True)
spl_r.columns = ['SNP', 'position']
ref = pd.concat([ref, spl_r['SNP']],axis=1).set_index('SNP')
ref_score = ref['score'].to_dict()
def dicr(x):
    return ref_score.get(x,0)

alt = pd.read_table('processed_data/variants_alt_fimo.txt', sep='\t')
spl_a = df['sequence name'].str.split(':A:wP', expand=True)
spl_a.columns = ['SNP', 'position']
alt = pd.concat([alt, spl_a['SNP']],axis=1).set_index('SNP')
alt_score = alt['score'].to_dict()
def dica(x):
    return alt_score.get(x,0)

df['score_ref'] = df['SNP'].apply(dicr)
df['score_alt'] = df['SNP'].apply(dica)
df = df.dropna(subset=['score_ref', 'score_alt'])
df['dBS'] = df['score_alt'] - df['scpre_ref']


### checking strand of motif
score = pd.read_table('processed_data/ref_all_fimo.txt', sep='\t')
MA1 = score[score['#pattern name'] == 'MA0138.1']
MA1l = MA1['sequence name'].tolist()
MA2 = score[score['#pattern name'] == 'MA0138.2']
MA2l = MA2['sequence name'].tolist()



## flip sequences with MA0138.1 motif
## 95th "C" of MA0138.1 = 107th "G" of MA0138.2

df['rel_pos'] = 0
df['color'] = '0.3'
df['segment'] = 'NA'

for i in range(len(df)):
    name = df.iloc[i,3]
    pos = df.iloc[i,18]
    if name in MA1l:
        df.loc[i,'rel_pos'] = 202 - pos
    elif name in MA2l:
        df.loc[i,'rel_pos'] = pos
    else:
        df.loc[i,'rel_pos'] = pos

    ## change color of significant emVars
    fdr = df.loc[i,'Skew_logFDR']
    if float(fdr) > 1.301:
        df.loc[i,'color'] = 'red'


### taking absolute value of logSkew
df['abs_skew'] = df['LogSkew'].abs()

df['dist'] = abs(101 - df['rel_pos'])
df = df.sort_values(by='dist')
df = df.drop_duplicates(subset=['SNP'])

### option to select only emVars
#df = df[df['Skew_logFDR'] > 1.301]

df = df[(df['rel_pos'] <= 126) & (df['rel_pos'] >= 76)]

### separate above and below for line plots
above = df[(df['Ref_BS'] > seg)| (df['Alt_BS'] > seg)]
below = df[(df['Ref_BS'] < seg) & (df['Alt_BS'] < seg)]



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
