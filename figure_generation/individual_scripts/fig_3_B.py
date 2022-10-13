import pandas as pd
import sys
import matplotlib.pyplot as plt
from PIL import Image
import math
import numpy as np
import seaborn as sns
from scipy import stats
import statsmodels.stats.multitest as multi

file = 'processed_data/RESTscreen_derived_K562_results.run1.txt'
l2fc = pd.read_table(file)

## select En19
l2fc = l2fc[l2fc['ID'].str.contains('En19')]

## separate results to canonical, negative control, and no-canonical
ref = l2fc[l2fc['ID'].str.contains('_Ref')]
Neg = l2fc[l2fc['project'] == 'NegCtrl']
l2fc = l2fc[l2fc['project'] == 'noMotif']

### process canonical motifs
#### make a dictionary to call score
refsc = pd.read_table('processed_data/ref_all_fimo.txt', sep='\t')
refsc = refsc.set_index('sequence name')
refsc_dic = refsc['score'].to_dict()
def Refscore(x):
    return refsc_dic.get(x,0)

#### select RE1 with strong canonical motif
idref = ref['ID'].str.split('^',expand=True)
idref.columns = ['enhancer','silencer']
ref = pd.concat([ref,idref],axis=1)
ref = ref.loc[:,['enhancer','silencer','log2FoldChange']].reset_index(drop=True)
ref['score'] = ref['silencer'].apply(Refscore)
ref = ref[ref['score'] > 20.86].dropna(subset= ['log2FoldChange'])

#### process nomotifs
ids = l2fc['ID'].str.split('^',expand=True)
ids.columns = ['enhancer','silencer']
l2fc = pd.concat([l2fc,ids],axis=1)
l2fc = l2fc.loc[:,['enhancer','silencer','log2FoldChange']].reset_index(drop=True)


## FIMO result for left and right half motif
left = pd.read_table('processed_data/halfmotif_left_fimo.txt')
right = pd.read_table('processed_data/halfmotif_right_fimo.txt')

### make dictionaries to call half-score
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

### make a table for ones with both L and R halfs
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

### filter weak motifs
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


for i in range(len(l2fc)):
    algn = l2fc.loc[i,'alignment']
    if algn == 'others':
        if l2fc.loc[i,'score_L'] > 10:
            g = 'left_only'
        elif l2fc.loc[i,'score_R'] > 10:
            g = 'right_only'
        else:
            g = 'No_half'
        l2fc.loc[i,'alignment'] = g
        l2fc.loc[i,'gap'] = 0

sortl = ['canonical','No_half','atypically_spaced','flipped','convergent','divergent','left_only','right_only']
l2fc['sort'] = l2fc['alignment'].apply(lambda x: sortl.index(x) if x in sortl else -1)
l2fc = l2fc.sort_values(by='sort')

l2fc = l2fc[l2fc['gap'] <= 100]


## add negative control
Neg['score_L'] = 0
Neg['score_R'] = 0
Neg['gap'] = 'NA'
Neg['alignment'] = 'NegCtrl'
negmed = Neg['log2FoldChange'].median()

## add canonical
ref['score_L'] = 0
ref['score_R'] = 0
ref['gap'] = 'NA'
ref['alignment'] = 'canonical'

l2fc = pd.concat([Neg, ref, l2fc],axis=0)


## plot
fig1 = plt.figure(figsize=(2,2), dpi=300)
fig1 = plt.xticks(fontsize= 'xx-small',rotation = 20)
fig1 = plt.yticks(fontsize= 'xx-small')
fig1 = plt.ylim((-3,8.5))
fig1 = plt.ylabel('log2FoldChange',fontsize='xx-small')
fig1 =  sns.boxplot(y=l2fc['log2FoldChange'], x=l2fc['alignment'],linewidth=1,zorder=1,showfliers = False)
fig1 = plt.axhline(y=negmed, color='grey', alpha=0.5,linewidth=1,zorder=0)
fig1 = plt.title('nomotifs expression (K562,En19)',fontsize='xx-small')
fig1 = plt.savefig('Fig3B_box.pdf')


## calc p and adjp values
p = []

for i in sortl:
    l2fc_nc = l2fc[l2fc['alignment'] == 'NegCtrl']
    l2fc_pd = l2fc[l2fc['alignment'] == i]

    pval = stats.mannwhitneyu(l2fc_nc['log2FoldChange'].dropna(),l2fc_pd['log2FoldChange'].dropna() , alternative='two-sided',use_continuity=True).pvalue
    p.append(pval)
    #print i,'\t',pval


adjp = multi.multipletests(np.array(p),alpha=0.05,method='fdr_bh')
adjp = adjp[1].tolist()
res = pd.Series(adjp,index=sortl)

res.to_csv('Fig3B_stat.txt',index=True,header=False,sep='\t')
