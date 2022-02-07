import pandas as pd
import sys
import math
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from PIL import Image

### import files
left = pd.read_table('processed_data/halfmotif_left_fimo.txt')
right = pd.read_table('processed_data/halfmotif_right_fimo.txt')


### make dataframe for ones with both L and R halfs
both = list(set(left['sequence_name'].values.tolist())  & set(right['sequence_name'].values.tolist()))

cmn_l = left[left['sequence_name'].isin(both)]
cmn_l = cmn_l.drop_duplicates(subset=['sequence_name']).sort_values(by='sequence_name').set_index('sequence_name')
cmn_l = cmn_l.drop(['motif_id','motif_alt_id','q-value'],axis=1)
cmn_l.columns = ['startL','stopL','strandL','scoreL','p-valueL','mat_seq_L']


cmn_r = right[right['sequence_name'].isin(both)]
cmn_r = cmn_r.drop_duplicates(subset=['sequence_name']).sort_values(by='sequence_name').set_index('sequence_name')
cmn_r = cmn_r.drop(['motif_id','motif_alt_id','q-value'],axis=1)
cmn_r.columns = ['startR','stopR','strandR','scoreR','p-valueR','mat_seq_R']

df = pd.concat([cmn_l,cmn_r],axis=1)


### cut off weak motifs
df['sum_score'] = df['scoreL'] + df['scoreR']
df = df[df['sum_score']>=20.86]

df = df.reset_index()
df['gap'] = 'NA'
df['alignment'] = 'NA'
df['strands'] = df['strandL'] + df['strandR']

## categorize non-canonical motifs
for i in range(len(df)):
    strands = df.loc[i,'strands']

    if df.loc[i,'stopL'] < df.loc[i,'startR']:
        gap = df.loc[i,'startR'] - df.loc[i,'stopL'] -1
        if  strands == '++':
            alignment = 'atypically_spaced'
        elif strands == '--' :
            alignment = 'flipped'
        elif strands == '+-':
            alignment = 'convergent'
        else:
            alignment = 'divergent'

    elif df.loc[i,'stopR'] < df.loc[i,'startL']:
        gap = df.loc[i,'startL'] - df.loc[i,'stopR'] -1
        if strands == '++':
            alignment = 'flipped'
        elif strands == '--' :
            alignment = 'atypically_spaced'
        elif strands == '+-':
            alignment = 'convergent'
        else:
            alignment = 'divergent'

    else:
        gap = 'NA'
        alignment = 'overlap'

    df.loc[i,'gap'] = gap
    df.loc[i,'alignment'] = alignment

df = df[(df['gap']>=0)&(df['gap']<=25)]
df = df.astype({'gap': int})


## plot
bins = np.arange(0, 25, 1)
t = 60
df = df[df['alignment'] != 'overlap']
LR = df[df['alignment'] == 'atypically_spaced']
LR = LR['gap'].values
RL = df[df['alignment'] == 'flipped']
RL = RL['gap'].values
inw = df[df['alignment'] == 'convergent']
inw = inw['gap'].values
outw = df[df['alignment'] == 'divergent']
outw = outw['gap'].values

fig1 = plt.figure(figsize=(4,4), dpi=300)
fig1 = plt.title('gap (bp) between half motifs', color='0.2')

ax1 = plt.subplot(221)
ax1 = plt.title('atypically_spaced',fontsize='xx-small')
ax1 = plt.ylim((0,t))
ax1 = plt.xlim((0,25))
ax1 = plt.hist(LR,bins)
ax1 = plt.xticks([0,25,50],fontsize='xx-small')
ax1 = plt.yticks(fontsize='xx-small')

ax2 = plt.subplot(222)
ax2 = plt.title('flipped',fontsize='xx-small')
ax2 = plt.ylim((0,t))
ax2 = plt.xlim((0,25))
ax2 = plt.hist(RL,bins)
ax2 = plt.xticks([0,25,50],fontsize='xx-small')
ax2 = plt.yticks(fontsize='xx-small')

ax3 = plt.subplot(223)
ax3 = plt.title('convergent',fontsize='xx-small')
ax3 = plt.ylim((0,t))
ax3 = plt.xlim((0,25))
ax3 = plt.hist(inw,bins)
ax3 = plt.xticks([0,25,50],fontsize='xx-small')
ax3 = plt.yticks(fontsize='xx-small')

ax4 = plt.subplot(224)
ax4 = plt.title('divergent',fontsize='xx-small')
ax4 = plt.ylim((0,t))
ax4 = plt.xlim((0,25))
ax4 = plt.hist(outw,bins)
ax4 = plt.xticks([0,25,50],fontsize='xx-small')
ax4 = plt.yticks(fontsize='xx-small')

fig1 = plt.savefig('Fig3C_hist.pdf' )
