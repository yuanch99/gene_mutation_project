import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from matplotlib.cm import get_cmap
DATAPATH = '../data/'
OUTPUTPATH = '../output/'

hyper_data = pd.read_csv('/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/raw_data/old_hyperdata.csv')

#load cluster labels
biased_pairs_hd = pd.read_csv(DATAPATH + 'biased_pairs_hd.csv',index_col=0)
unbiased_pairs_hd = pd.read_csv(DATAPATH + 'unbiased_pairs_hd.csv',index_col=0)
cross_cluster_hd = pd.read_csv(DATAPATH + 'cross_cluster_hd.csv',index_col=0)

orientation_data = pd.read_csv(DATAPATH + 'mean_std_log2_all_samples_complete.csv',index_col=0)
tmb_discrepancy_all = pd.read_csv(DATAPATH + 'tmb_discrepancy_all.csv',index_col=0)



def pltClusterCompare(c,name,show_legend=False,palette=None):
    figs,axes = plt.subplots()
    ax1 = plt.subplot(1,1,1)
    if palette is None:
        palette = 'viridis'
    else:   
        palette = palette
    sns.scatterplot(
        x='mean_log2_ratio_F1R2_over_F2R1_50plus',
        y='std_log2_ratio_F1R2_over_F2R1_50plus',
        hue=c,
        data=orientation_data,
        ax=ax1,
        legend=show_legend,
        size=10,
        alpha=1,
        palette=palette
    )
    ax1.set_title(name)
    ax1.set_xlabel('Mean')
    ax1.set_ylabel('Std')
    figs.savefig(OUTPUTPATH + 'mean_std_log2_all_samples_' + name + '.png',dpi=300)
    return figs,axes

# sample sequencing batch
cluster_seq_data = pd.read_csv(DATAPATH + 'cluster_seq_data.csv')
cluster_seq_data.replace(np.nan,'Unknown',inplace=True)
seq_batch = []
for i in orientation_data['sample']:
    seq_batch.append(
        cluster_seq_data[cluster_seq_data['CorSampleName'] == i]['BatchID'].values[0])
cmap_seq = get_cmap('hsv',len(set(seq_batch)))
cmap_seq_dict = dict(zip(set(seq_batch),cmap_seq(np.arange(0,cmap_seq.N))))
figs1, axes1 = pltClusterCompare(seq_batch,'seq_batch',False,palette=cmap_seq_dict)

#pie chart of batch for biased samples
seq_batch_biased = []
for i in np.unique(np.append(biased_pairs_hd.values.flatten(),cross_cluster_hd.values.flatten())):
    seq_batch_biased.append(
        cluster_seq_data[cluster_seq_data['CorSampleName'] == i]['BatchID'].values[0])
seq_batch_biased = Counter(seq_batch_biased)
color_seq_batch_biased = []
for i in seq_batch_biased.keys(): color_seq_batch_biased.append(cmap_seq_dict[i])

seq_batch_unbiased = []
for i in np.unique(unbiased_pairs_hd.values.flatten()):
    seq_batch_unbiased.append(
        cluster_seq_data[cluster_seq_data['CorSampleName'] == i]['BatchID'].values[0])
seq_batch_unbiased = Counter(seq_batch_unbiased)
color_seq_batch_unbiased = []
for i in seq_batch_unbiased.keys(): color_seq_batch_unbiased.append(cmap_seq_dict[i])

figs2,axes2 = plt.subplots(figsize=(15,8))
def my_fmt(x):
    return '{:.0f}'.format(x, total*x/100)
ax1 = plt.subplot(1,2,1)
total = sum(seq_batch_biased.values())
ax1.pie(
    seq_batch_biased.values(),
    labels=seq_batch_biased.keys(),
    colors=color_seq_batch_biased,
    autopct=my_fmt)
ax1.set_title(f'Biased and half-biased,n={total}')
ax2 = plt.subplot(1,2,2)
total = sum(seq_batch_unbiased.values())
ax2.pie(
    seq_batch_unbiased.values(),
    labels=seq_batch_unbiased.keys(),
    colors=color_seq_batch_unbiased,
    autopct=my_fmt)
ax2.set_title(f'Unbiased,n={total}')
figs2.savefig(OUTPUTPATH + 'seq_batch_pie.png',dpi=300)


figs3,axes3 = plt.subplots(figsize=(21,8))
tmb_discrepancy_all = pd.read_csv(DATAPATH+'tmb_discrepancy_all.csv')
for i,j in zip(range(1,4),['tmb_discrepancy_withOB','tmb_discrepancy_noOB','tmb_discrepancy_OB_threshold']):
    ax2 = plt.subplot(1,3,i)
    for k,l in enumerate(tmb_discrepancy_all['tumor'].tolist()):
        if tmb_discrepancy_all['bias'].values[k] == 'biased': continue #color = '#e4455e'
        elif tmb_discrepancy_all['bias'].values[k] == 'unbiased': color = '#86a1a9'
        elif tmb_discrepancy_all['bias'].values[k] == 'cross_cluster': continue #color = '#b18f6a'
        diff = (tmb_discrepancy_all[j].values[k]-tmb_discrepancy_all['old_tmb'].values[k])
        if cluster_seq_data[cluster_seq_data['CorSampleName']==l]['BatchID'].values[0] in seq_batch_biased.keys():
            style = 'dotted'
        else:
            style = '-'mod
        plt.plot([0,diff],[k,k],color=color,linestyle=style)
    ax2.set_title(f'{j},midian=%.2f,mean=%.2f'%(np.median(abs(tmb_discrepancy_all[j].values-tmb_discrepancy_all['old_tmb'].values)),np.mean(abs(tmb_discrepancy_all[j].values-tmb_discrepancy_all['old_tmb'].values))))
    ax2.set_xlabel(f'tmb_discrepancy')
    ax2.set_ylabel('tumor')
    ax2.set_ylim(0,len(tmb_discrepancy_all))
    ax2.set_xlim(-400,400)
figs3.savefig(OUTPUTPATH + 'tmb_bad_batches_on_unbiased.png',dpi=300,bbox_inches="tight")

