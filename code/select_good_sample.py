import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import shapiro
from numpy.random import normal
from scipy.stats import ttest_ind
import seaborn as sns
DATAPATH = '/Users/yuanchang/tabori/codes/mutation_project/gene_mutation_project/vcf_validation/data/orientation_counts'
OUTPUTPATH = '/Users/yuanchang/tabori/codes/mutation_project/gene_mutation_project/vcf_validation/output/orientation_counts'
# get all files name ending with read_orientation.target_count.tsv in the DATAPATH
files = [f for f in os.listdir(DATAPATH) if f.endswith('read_orientation.target_count.tsv')]

sample_good_bad = pd.read_csv('/Users/yuanchang/tabori/codes/mutation_project/gene_mutation_project/vcf_validation/data/sample_good_bad.csv')

# create a dataframe to with sample names
sample_list = pd.DataFrame(columns=[
    'sample_name',
    'mean','std','stat_shapiro','p_shapiro','shapiro','stat_t','p_t','t',
    'mean_log2','std_log2','stat_shapiro_log2','p_shapiro_log2','shapiro_log2','stat_t_log2','p_t_log2','t_log2',
    'tmb_guess'
    ])
for i in range(len(files)):
    # read the file
    df = pd.read_csv(os.path.join(DATAPATH,files[i]),sep='\t')

    # add column with F1R2/F2R1
    df['F1R2/F2R1'] = (df['F1R2']/df['F2R1']-1)

    #remove nan
    df = df.dropna()

    # replace inf with max value in the column
    df['F1R2/F2R1'].replace(
        [np.inf,-np.inf],
        [df.loc[df['F1R2/F2R1'] != np.inf, 'F1R2/F2R1'].max(),
        df.loc[df['F1R2/F2R1'] != np.inf, 'F1R2/F2R1'].min()],
        inplace=True)

    #shapiro test
    # check if the distribution is normal
    stat_shapiro, p_shapiro = shapiro(df['F1R2/F2R1'])
    print('Shapiro Statistics=%.3f, p=%.3f' % (stat_shapiro, p_shapiro))
    alpha = 0.05
    if p_shapiro > alpha:
        status_shapiro = True
    else:
        status_shapiro = False

    #t test
    norm = normal(0,df['F1R2/F2R1'].std(),len(df['F1R2/F2R1']))
    stat_t,p_t = ttest_ind(df['F1R2/F2R1'],norm)
    print('t Statistics=%.3f, p=%.3f' % (stat_t, p_t))
    # interpret

    if p_t > alpha:
        # store pass result to dataframe       
        status_t = True
    else:
        # store fail result to dataframe
        status_t = False


    #same tests for log2_F1R2_over_F2R1
    df = df.dropna()
    stat_shapiro_log2, p_shapiro_log2 = shapiro(df['log2_F1R2_over_F2R1'])
    print('Shapiro Statistics=%.3f, p=%.3f' % (stat_shapiro_log2, p_shapiro_log2))
    if p_shapiro_log2 > alpha:
        status_shapiro_log2 = True
    else:
        status_shapiro_log2 = False
    norm_log2 = normal(0,df['log2_F1R2_over_F2R1'].std(),len(df['log2_F1R2_over_F2R1']))
    stat_t_log2,p_t_log2 = ttest_ind(df['log2_F1R2_over_F2R1'],norm_log2)
    print('t Statistics=%.3f, p=%.3f' % (stat_t_log2, p_t_log2))
    if p_t_log2 > alpha:
        status_t_log2 = True
    else:
        status_t_log2 = False

    #write to dataframe
    #if sample name in sample_good_bad in tumor column or normal column is good, then good = True
    if sample_good_bad[sample_good_bad[['tumor','normal']].eq(files[i].split('.')[0]).any(axis=1)]['guess_from_tmb'].values[0] == 'good':
        good = 'good'
    else:
        good = 'bad'

    sample_list.loc[len(sample_list)] = [
        files[i].split('.')[0],
        df['F1R2/F2R1'].mean(),
        df['F1R2/F2R1'].std(),
        stat_shapiro,p_shapiro,
        status_shapiro,
        stat_t,p_t,
        status_t,
        df['log2_F1R2_over_F2R1'].mean(),
        df['log2_F1R2_over_F2R1'].std(),
        stat_shapiro_log2,p_shapiro_log2,
        status_shapiro_log2,
        stat_t_log2,p_t_log2,
        status_t_log2,
        good
    ]

    # plot the F1R2/F2R1
    _ = plt.figure()
    _ = plt.hist(df['F1R2/F2R1'],bins=500,color='#ff6961',label='F1R2/F2R1',density=False)
    _ = plt.hist(df['log2_F1R2_over_F2R1'],bins=500,color='#ffb480',label='log2_F1R2/F2R1',density=False)
    _ = plt.hist(norm,bins=500,color='#42d6a4',label='Normal',density=False)
    _ = plt.hist(norm_log2,bins=500,color='#59adf6',label='Normal_log2',density=False)
    _ = plt.legend(title=f'mean={df["F1R2/F2R1"].mean():.3f}, std={df["F1R2/F2R1"].std():.3f},\nt_stat={stat_t:.3f}')
    _ = plt.title(files[i].split('.')[0])
    _ = plt.xlabel('F1R2/F2R1')
    _ = plt.ylabel('Density')
    _ = plt.xlim(-3,4)
    #_ = plt.savefig(os.path.join(OUTPUTPATH,'figs',files[i].split('.')[0]+'.png'),dpi=200)
    _ = plt.close()

#save the dataframe to csv
#sample_list.to_csv(os.path.join(OUTPUTPATH,'orientation_counts_normality_test.csv'),index=False)


#plot features
figs,axes=plt.subplots(figsize=(30,20),nrows=3,ncols=6)
axes=axes.flatten()
for i in range(1,17):
    ax=axes[i]
    feature=sample_list.columns[i]
    
    sns.barplot(x='tmb_guess', y=feature,data=sample_list, capsize=.5,ci='sd',palette=['#ff6961','#08cad1'],ax=ax)
    sns.swarmplot(x='tmb_guess', y=feature,data=sample_list, color='b', alpha=0.8,ax=ax,size=4)
    ax.set_xlabel('guess_from_tmb')
    ax.set_ylabel(feature)
    # _ = plt.tight_layout()
    t,p = ttest_ind(sample_list[feature][sample_list['tmb_guess']=='good'],sample_list[feature][sample_list['tmb_guess']=='bad'])
    if p < 0.05:
        ax.set_title(f'{feature}* , p={p:.3f}')
    else:
        ax.set_title(f'{feature} , p={p:.3f}')
#plt.savefig(os.path.join(OUTPUTPATH,'figs','orientation_counts_normality_test_barplot.png'),dpi=400)

