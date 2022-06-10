import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.cluster import KMeans
import seaborn as sns

SERVER = 'local'

if SERVER == 'HPC':
    DATAPATH1 = '/hpf/largeprojects/tabori/users/yuan/wes/santiago_pipeline_results/from_bam/orientation/'
    DATAPATH2 = '/hpf/largeprojects/tabori/users/yuan/wes/santiago_pipeline_results/from_fastq/orientation/'
    OUTPUTPATH = '/hpf/largeprojects/tabori/users/yuan/wes/santiago_pipeline_results/vcf_validation/output/'
    print('loaded')

    # Load data
    # Load all tsv files ending in '.read_orientation.summary.tsv'
    files = [f for f in os.listdir(DATAPATH1) if f.endswith('.read_orientation.summary.tsv')]
    files = [DATAPATH1+f for f in files]
    files2 = [f for f in os.listdir(DATAPATH2) if f.endswith('.read_orientation.summary.tsv')]
    files2 = [DATAPATH2+f for f in files2]
    files=files+files2

    #remove emplty files
    files = [f for f in files if os.stat(f).st_size != 0]

    # and concatenate them into a single dataframe
    data = pd.concat([pd.read_csv(f, sep='\t') for f in files], ignore_index=True)

if SERVER == 'local':
    DATAPATH = '../data/'
    OUTPUTPATH = '../output/'
    data = pd.read_csv(DATAPATH + 'mean_std_log2_all_samples.csv')
    tumors_normals = pd.read_csv('/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/raw_data/tumors_and_normals.csv',header=None)

#remove rows contain element list
exclued_list = ['D1410','MMR111B1','D1423','MMR120B1','D1577','M139B1','MD0084T8','MD0084B1','D1244','D1805_1','MD101T10','MMR101B1','MD57T1']
data = data[~data.iloc[:,0].isin(exclued_list)]
data.to_csv(DATAPATH + 'mean_std_log2_all_samples_complete.csv')

#plot mean_log2_ratio_F1R2_over_F2R1_50plus vs std_log2_ratio_F1R2_over_F2R1_50plus
figs,axes = plt.subplots()
ax1 = plt.subplot(1,1,1)
ax1.scatter(
    data['mean_log2_ratio_F1R2_over_F2R1_50plus'], 
    data['std_log2_ratio_F1R2_over_F2R1_50plus'], 
    color='#08cad1',
    alpha=0.5)
ax1.set_title('Mean vs Std of log2(F1R2/F2R1 50+)')
ax1.set_xlabel('Mean')
ax1.set_ylabel('Std')
#save figure
figs.savefig(OUTPUTPATH + 'mean_std_log2_all_samples.png',dpi=300)

data.to_csv(OUTPUTPATH + 'mean_std_log2_all_samples.csv', index=False)

#KMeans
data = data.dropna()
km = KMeans(n_clusters=2, random_state=0).fit(
    data[['mean_log2_ratio_F1R2_over_F2R1_50plus', 'std_log2_ratio_F1R2_over_F2R1_50plus']],
    )
labels = km.labels_
#1 is bias, 0 is unbias
biased_data = data[labels.astype(bool).tolist()]
unbiased_data = data[(~labels.astype(bool)).tolist()]

#connect tumor and normal
#check if the tumor and normal are in the same cluster

#remove rows if not in data
not_finished = []
finished = []
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in data.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in data.iloc[:,0].tolist():
            finished.append(i)
            continue
    not_finished.append(i)
tumors_normals.drop(not_finished, inplace=True)




cross_cluster = []
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in biased_data.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in biased_data.iloc[:,0].tolist():
            continue
    elif tumors_normals.iloc[i,0] in unbiased_data.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in unbiased_data.iloc[:,0].tolist():
            continue
    cross_cluster.append(tumors_normals.iloc[i,:].tolist())
    
#plot cross_cluster samples 
fig2,axes2 = plt.subplots(figsize=(10,10))
ax1 = plt.subplot(1,1,1)
ax1.scatter(
    data['mean_log2_ratio_F1R2_over_F2R1_50plus'],
    data['std_log2_ratio_F1R2_over_F2R1_50plus'],
    c=labels,
    cmap='coolwarm',
    alpha=0.5)
ax1.set_title('KMeans Clustering, color = clutser')
ax1.set_xlabel('Mean')
ax1.set_ylabel('Std')
fig2.savefig(OUTPUTPATH + 'mean_std_log2_all_samples_kmeans.png',dpi=300)

for i in cross_cluster:
    tum_loc = data.loc[data['sample']==i[0],['mean_log2_ratio_F1R2_over_F2R1_50plus','std_log2_ratio_F1R2_over_F2R1_50plus']].values.flatten()
    norm_loc = data.loc[data['sample']==i[1],['mean_log2_ratio_F1R2_over_F2R1_50plus','std_log2_ratio_F1R2_over_F2R1_50plus']].values.flatten()
    print(tum_loc,norm_loc)
    ax1.plot([tum_loc[0],norm_loc[0]],[tum_loc[1],norm_loc[1]],'--',color='#9d94ff')
    ax1.annotate(i[0], tum_loc, fontsize=8)
    ax1.annotate(i[1], norm_loc, fontsize=8)
    
#save figure
fig2.savefig(OUTPUTPATH + 'mean_std_log2_all_samples_kmeans_cross_cluster.png',dpi=300)

unbiased_pairs = []
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in unbiased_data.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in unbiased_data.iloc[:,0].tolist():
            unbiased_pairs.append(tumors_normals.iloc[i,:].tolist())
    
biased_pairs = []
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in biased_data.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in biased_data.iloc[:,0].tolist():
            biased_pairs.append(tumors_normals.iloc[i,:].tolist())

half_biased_pairs = []
#this is the same as cross_cluster
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in biased_data.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in unbiased_data.iloc[:,0].tolist():
            half_biased_pairs.append(tumors_normals.iloc[i,:].tolist())

    elif tumors_normals.iloc[i,1] in biased_data.iloc[:,0].tolist():
        if tumors_normals.iloc[i,0] in unbiased_data.iloc[:,0].tolist():
            half_biased_pairs.append(tumors_normals.iloc[i,:].tolist())



##### Hard threshold
unbiased_data_hd = data[(data['mean_log2_ratio_F1R2_over_F2R1_50plus']>-1.5) & (data['std_log2_ratio_F1R2_over_F2R1_50plus']<0.6)]
biased_data_hd = data[(data['mean_log2_ratio_F1R2_over_F2R1_50plus']<=-1.5) | (data['std_log2_ratio_F1R2_over_F2R1_50plus']>=0.6)]
labels_hd = []
for i in range(data.shape[0]):
    if data.iloc[i,4]>-1.5 and data.iloc[i,5]<0.6:
        labels_hd.append(0)
    else:
        labels_hd.append(1)

cross_cluster_hd = []
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in biased_data_hd.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in biased_data_hd.iloc[:,0].tolist():
            continue
    elif tumors_normals.iloc[i,0] in unbiased_data_hd.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in unbiased_data_hd.iloc[:,0].tolist():
            continue
    cross_cluster_hd.append(tumors_normals.iloc[i,:].tolist())

unbiased_pairs_hd = []
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in unbiased_data_hd.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in unbiased_data_hd.iloc[:,0].tolist():
            unbiased_pairs_hd.append(tumors_normals.iloc[i,:].tolist())

biased_pairs_hd = []
for i in range(tumors_normals.shape[0]):
    if tumors_normals.iloc[i,0] in biased_data_hd.iloc[:,0].tolist():
        if tumors_normals.iloc[i,1] in biased_data_hd.iloc[:,0].tolist():
            biased_pairs_hd.append(tumors_normals.iloc[i,:].tolist())

#plot hard threshold
fig3,axes3 = plt.subplots(figsize=(10,10))
ax1 = plt.subplot(1,1,1)
ax1.scatter(
    data['mean_log2_ratio_F1R2_over_F2R1_50plus'],
    data['std_log2_ratio_F1R2_over_F2R1_50plus'],
    c=labels_hd,
    cmap='coolwarm',
    alpha=0.5)
ax1.vlines(x=-1.5, ymin=ax1.get_ylim()[0], ymax=0.6, color='#ff6961', linestyle='--')
ax1.hlines(y=0.6, xmin=-1.5, xmax=ax1.get_xlim()[1], color='#ff6961', linestyle='--')
ax1.set_title('Hard Threshold, mean > -1.5 and std < 0.6')
ax1.set_xlabel('Mean')
ax1.set_ylabel('Std')
fig3.savefig(OUTPUTPATH + 'mean_std_log2_all_samples_hard_threshold.png',dpi=300)

for i in cross_cluster_hd:
    tum_loc = data.loc[data['sample']==i[0],['mean_log2_ratio_F1R2_over_F2R1_50plus','std_log2_ratio_F1R2_over_F2R1_50plus']].values.flatten()
    norm_loc = data.loc[data['sample']==i[1],['mean_log2_ratio_F1R2_over_F2R1_50plus','std_log2_ratio_F1R2_over_F2R1_50plus']].values.flatten()
    ax1.plot([tum_loc[0],norm_loc[0]],[tum_loc[1],norm_loc[1]],'--',color='#9d94ff')
    ax1.annotate(i[0], tum_loc, fontsize=8)
    ax1.annotate(i[1], norm_loc, fontsize=8)
    
#save figure
fig3.savefig(OUTPUTPATH + 'mean_std_log2_all_samples_hardthreshold_cross_cluster.png',dpi=300)

#save data
pd.DataFrame(cross_cluster,columns=['tumor','normal']).to_csv(DATAPATH + 'cross_cluster.csv')
pd.DataFrame(unbiased_pairs,columns=['tumor','normal']).to_csv(DATAPATH + 'unbiased_pairs.csv')
pd.DataFrame(biased_pairs,columns=['tumor','normal']).to_csv(DATAPATH + 'biased_pairs.csv')
pd.DataFrame(cross_cluster_hd,columns=['tumor','normal']).to_csv(DATAPATH + 'cross_cluster_hd.csv')
pd.DataFrame(unbiased_pairs_hd,columns=['tumor','normal']).to_csv(DATAPATH + 'unbiased_pairs_hd.csv')
pd.DataFrame(biased_pairs_hd,columns=['tumor','normal']).to_csv(DATAPATH + 'biased_pairs_hd.csv')
tumors_normals.columns=['tumor','normal']

tumors_normals.to_csv(DATAPATH + 'tumors_normals_finished_completely.csv')

#FFPE vs Frozen
figs4,axes4 = plt.subplots()
ax1 = plt.subplot(1,1,1)
hyper_data = pd.read_csv('/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/raw_data/old_hyperdata.csv')
format_list = []
for i in data['sample']:
    if i not in hyper_data['TumorName'].tolist():
        format = 'blood'
    else:
        if 'FFPE' in hyper_data.loc[hyper_data['TumorName']==i,'TissueFormat'].values[0]:
            format = 'FFPE'
        elif 'Frozen' in hyper_data.loc[hyper_data['TumorName']==i,'TissueFormat'].values[0]:
            format = 'Frozen'
        else:
            format = hyper_data.loc[hyper_data['TumorName']==i,'TissueFormat'].values[0]
    format_list.append(format)

sns.scatterplot(
    x='mean_log2_ratio_F1R2_over_F2R1_50plus',
    y='std_log2_ratio_F1R2_over_F2R1_50plus',
    hue=format_list,
    data=data,
    ax=ax1)
ax1.set_title('FFPE vs Frozen')
ax1.set_xlabel('Mean')
ax1.set_ylabel('Std')
figs4.savefig(OUTPUTPATH + 'mean_std_log2_all_samples_FFPE_vs_Frozen.png',dpi=300)
