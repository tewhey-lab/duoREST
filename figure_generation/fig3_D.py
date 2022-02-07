import pandas as pd
import sys
import matplotlib.pyplot as plt
from PIL import Image
import math
import numpy as np
import seaborn as sns
from scipy import stats
import statsmodels.stats.multitest as multi


l2fc = pd.read_table('processed_data/RESTscreen_derived_K562_results.run1.txt')

### select negative controls and no-canonicals
l2fc = l2fc[l2fc['ID'].str.contains('En19')]
Neg = l2fc[l2fc['project'] == 'NegCtrl']
l2fc = l2fc[l2fc['project'] == 'noMotif']


ids = l2fc['ID'].str.split('^',expand=True)
ids.columns = ['enhancer','silencer']
l2fc = pd.concat([l2fc,ids],axis=1)
l2fc = l2fc.loc[:,['enhancer','silencer','log2FoldChange']].reset_index(drop=True)

### inport half-scores
left = pd.read_table('processed_data/halfmotif_left_fimo.txt')
right = pd.read_table('processed_data/halfmotif_right_fimo.txt')

#### make dictionaries to call half-scores
lsc = left.loc[:,['sequence_name','score']].sort_values('score',ascending=False)
lsc = lsc.drop_duplicates(subset=['sequence_name'])
lsc = lsc.set_index('sequence_name')
dicL = lsc['score'].to_dict()
def Lscore(x):
    return dicL.get(x,0)

rsc = right.loc[:,['sequence_name','score']].sort_values('score',ascending=False)
rsc = rsc.drop_duplicates(subset=['sequence_name'])
rsc = rsc.set_index('sequence_name')
dicR = rsc['score'].to_dict()
def Rscore(x):
    return dicR.get(x,0)

l2fc['score_L'] = l2fc['silencer'].apply(Lscore)
l2fc['score_R'] = l2fc['silencer'].apply(Rscore)

### make a dataframe for ones with both L and R halfs
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
            alignment = 'atypically_spaced'
        elif strands == '--' :
            alignment = 'flipped'
        elif strands == '+-':
            alignment = 'convergent'
        else:
            alignment = 'divergent'

    else:
        gap = 'NA'
        alignment = 'overlap'

    df.loc[i,'gap'] = gap
    df.loc[i,'alignment'] = alignment

df = df[(df['gap']>=0)&(df['gap']<200)]
df = df.astype({'gap': int})

df = df.set_index('sequence_name')


dic_gap = df['gap'].to_dict()
def gap(x):
    return dic_gap.get(x,'NA')

dic_algn = df['alignment'].to_dict()
def algn(x):
    return dic_algn.get(x,'others')

l2fc['gap'] = l2fc['silencer'].apply(gap)
l2fc['alignment'] = l2fc['silencer'].apply(algn)

## select atypically spaced strong non-canonical motifs with <25 bp gap
l2fc = l2fc[l2fc['alignment'] == 'atypically_spaced']
l2fc = l2fc[l2fc['gap'] <= 25]
l2fc = l2fc.sort_values(by='gap').reset_index(drop=True)
l2fc = l2fc[(l2fc['score_L'] + l2fc['score_R']) > 20.86]

## add dummy data to plot gap=5 and gap=6
l2fc = l2fc.loc[:,['ID','alignment','gap','score_L','score_R','log2FoldChange']]
gap5 = pd.Series(['dummy','LR',5,0,0,7], index=l2fc.columns)
gap6 = pd.Series(['dummy','LR',6,0,0,7], index=l2fc.columns)
l2fc =l2fc.append(gap5, ignore_index=True)
l2fc =l2fc.append(gap6, ignore_index=True)

l2fc = l2fc.sort_values(by='gap')

## concatenate long gap motif
for i in range(len(l2fc)):
    if l2fc.loc[i,'gap'] >= 12:
        l2fc.loc[i,'gap'] = '>=12'

## add negative control
Neg['score_L'] = 0
Neg['score_R'] = 0
Neg['gap'] = 'Neg'
Neg['alignment'] = 'NegCtrl'
negmed = Neg['log2FoldChange'].median()
Neg = Neg.loc[:,['ID','alignment','gap','score_L','score_R','log2FoldChange']]

l2fc = pd.concat([Neg,l2fc],axis=0)

l = ['Neg',0,1,2,3,4,5,6,7,8,9,10,11,'>=12']
l2fc = l2fc[l2fc['gap'].isin(l)]

## plot
fig1 = plt.figure(figsize=(2.25,2), dpi=300)
fig1 = plt.xlabel('')
fig1 = plt.xticks(fontsize= 'xx-small',rotation = 30)
fig1 = plt.yticks(fontsize= 'xx-small')
fig1 = plt.ylim((-3,8.5))
fig1 =  sns.boxplot(y=l2fc['log2FoldChange'], x=l2fc['gap'],linewidth=1,zorder=1,showfliers = False)
fig1 = plt.axhline(y=negmed, color='grey', alpha=0.5,linewidth=1,zorder=0)
fig1 = plt.title('nomotifs expression (K562,En19)')
fig1 = plt.savefig('Fig3D_gapbox.pdf')


## calc p and adj values
Neg = Neg['log2FoldChange'].dropna().values
p = []
gap = []
for i in range(26):
    test = l2fc[l2fc['gap'] == i]
    if len(test) > 3:
        test = test['log2FoldChange'].dropna().values
        pval = stats.mannwhitneyu(Neg,test , alternative='two-sided',use_continuity=True).pvalue
        p.append(pval)
        gap.append(i)

adjp = multi.multipletests(np.array(p),alpha=0.05,method='fdr_bh')
adjp = adjp[1].tolist()
res = pd.Series(adjp,index=gap)

res.to_csv('Fig3D_stat.txt',index=True,header=False,sep='\t')
