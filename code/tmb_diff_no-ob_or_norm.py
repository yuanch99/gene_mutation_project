from cProfile import label
from turtle import color
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.patches import  mlines
from matplotlib_venn import venn3, venn2

DATAPATH = '../data/'
OUTPUTPATH = '../output/'
hyper_data = pd.read_csv('/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/raw_data/old_hyperdata.csv')
biased_pairs_hd = pd.read_csv(DATAPATH + 'biased_pairs_hd.csv',index_col=0)
unbiased_pairs_hd = pd.read_csv(DATAPATH + 'unbiased_pairs_hd.csv',index_col=0)
cross_cluster_hd = pd.read_csv(DATAPATH + 'cross_cluster_hd.csv',index_col=0)
orientation_data = pd.read_csv(DATAPATH + 'mean_std_log2_all_samples_complete.csv',index_col=0)
tmb_discrepancy_all = pd.read_csv(DATAPATH + 'tmb_discrepancy_all.csv',index_col=0)

###
overlap_data_norm = pd.read_csv(DATAPATH + 'overlap_counts_norm.csv')
overlap_data_noob_hd5001 = pd.read_csv(DATAPATH + 'overlap_counts_no-ob_hd5001.csv')
overlap_data_noob_hd2005 = pd.read_csv(DATAPATH + 'overlap_counts_no-ob_hd2005.csv')
overlap_data_noob_hd2005both = pd.read_csv(DATAPATH + 'overlap_counts_no-ob_hd2005both.csv')
overlap_data_noob_hd5010both = pd.read_csv(DATAPATH + 'overlap_counts_no-ob_hd5010both.csv')
overlap_data_noob_hd2010both = pd.read_csv(DATAPATH + 'overlap_counts_no-ob_hd2010both.csv')
overlap_data_varscan = pd.read_csv(DATAPATH + 'overlap_counts_varscan.csv')
###

tmb_recal = pd.read_csv(DATAPATH + 'tmb_re-cal_norm.csv')
tmb_recal_noob = pd.read_csv(DATAPATH + 'tmb_re-cal_no-ob_hd5001.csv')


bias_status_overlap = []
for i in overlap_data_noob_hd5001.tumor:
    if i in biased_pairs_hd.tumor.values:
        bias_status_overlap.append('biased')
    elif i in unbiased_pairs_hd.tumor.values:
        bias_status_overlap.append('unbiased')
    elif i in cross_cluster_hd.tumor.values:
        bias_status_overlap.append('cross-cluster')

###
overlap_data_norm['bias_status'] = bias_status_overlap
overlap_data_noob_hd5001['bias_status'] = bias_status_overlap
overlap_data_noob_hd2005['bias_status'] = bias_status_overlap
overlap_data_noob_hd2005both['bias_status'] = bias_status_overlap
overlap_data_noob_hd5010both['bias_status'] = bias_status_overlap
overlap_data_noob_hd2010both['bias_status'] = bias_status_overlap
overlap_data_varscan['bias_status'] = bias_status_overlap
###

###
#Remove MD1665T01
overlap_data_norm = overlap_data_norm[overlap_data_norm['tumor']!='MD1665T01']
overlap_data_noob_hd5001 = overlap_data_noob_hd5001[overlap_data_noob_hd5001['tumor']!='MD1665T01']
overlap_data_noob_hd2005 = overlap_data_noob_hd2005[overlap_data_noob_hd2005['tumor']!='MD1665T01']
overlap_data_noob_hd2005both = overlap_data_noob_hd2005both[overlap_data_noob_hd2005both['tumor']!='MD1665T01']
overlap_data_noob_hd5010both = overlap_data_noob_hd5010both[overlap_data_noob_hd5010both['tumor']!='MD1665T01']
overlap_data_noob_hd2010both = overlap_data_noob_hd2010both[overlap_data_noob_hd2010both['tumor']!='MD1665T01']
overlap_data_varscan = overlap_data_varscan[overlap_data_varscan['tumor']!='MD1665T01']
###
#plot tmb discrepancy based on no-ob

def plotTmbDiscrepancy(data,title,tmb_recal_noob=tmb_recal_noob):
    fig, ax = plt.subplots(figsize=(6,11))
    ax1 = plt.subplot(111)
    k=0
    unbias_data = data[data['bias_status']=='unbiased']
    color = '#86a1a9'
    for i,j in enumerate(unbias_data.tumor.values):
        new_tmb = (data[data['tumor']==j]['uni_in_new'].values[0]+data[data['tumor']==j]['overlap'].values[0]-data[data['tumor']==j]['ori_in_overlap'].values[0]-data[data['tumor']==j]['ori_in_uni_new'].values[0])/50
        old_tmb = (data[data['tumor']==j]['uni_in_old'].values[0]+data[data['tumor']==j]['overlap'].values[0])/50
        record_tmb = tmb_discrepancy_all[tmb_discrepancy_all['tumor']==j]['old_tmb'].values[0]
        diff = new_tmb - old_tmb
        # if abs(diff) > 100:
        if abs(diff/new_tmb) > 0.2:
            if max(new_tmb+10,old_tmb+10) < 800:
                ax1.annotate(j,xy=(max(new_tmb+10,old_tmb+10),k))
            elif max(new_tmb+10,old_tmb+10) >= 800:
                ax1.annotate(j,xy=(700,k))
        ax1.plot([0,new_tmb],[k,k],color=color,linewidth=0.75,ls='-')
        ax1.plot([0,old_tmb],[k+0.5,k+0.5],color=color,linewidth=0.75,ls='--')
        ax1.scatter(record_tmb,k+0.5,color='black',s=3,alpha=0.5)
        k += 2

    cross_data = data[data['bias_status']=='cross-cluster']
    color = '#b18f6a'
    for i,j in enumerate(cross_data.tumor.values):
        new_tmb = (data[data['tumor']==j]['uni_in_new'].values[0]+data[data['tumor']==j]['overlap'].values[0]-data[data['tumor']==j]['ori_in_overlap'].values[0]-data[data['tumor']==j]['ori_in_uni_new'].values[0])/50
        old_tmb = (data[data['tumor']==j]['uni_in_old'].values[0]+data[data['tumor']==j]['overlap'].values[0])/50
        record_tmb = tmb_discrepancy_all[tmb_discrepancy_all['tumor']==j]['old_tmb'].values[0]
        diff = new_tmb - old_tmb
        # if abs(diff) > 100:
        if abs(diff/new_tmb) > 0.2:
            if max(new_tmb+10,old_tmb+10) < 800:
                ax1.annotate(j,xy=(max(new_tmb+10,old_tmb+10),k))
            elif max(new_tmb+10,old_tmb+10) >= 800:
                ax1.annotate(j,xy=(700,k))
        ax1.plot([0,new_tmb],[k,k],color=color,linewidth=0.75,ls='-')
        ax1.plot([0,old_tmb],[k+0.5,k+0.5],color=color,linewidth=0.75,ls='--')
        ax1.scatter(record_tmb,k+0.5,color='black',s=3,alpha=0.5)
        k += 2
    ax1.set_xlim([1,800])
    ax1.set_ylim([-1,k])
    ax1.set_xlabel('TMB')
    ax1.set_ylabel('Tumor#')
    ax1.set_yticks([])
    ax1.set_title(title)
    plt.legend(['New TMB','Old TMB','Record TMB'])
    fig.savefig(OUTPUTPATH + f'{title}_unbia_cross.png',dpi=300)

    k=0
    fig, ax = plt.subplots(figsize=(6,6))
    ax2 = plt.subplot(111)
    bias_data = data[data['bias_status']=='biased']
    color = '#e4455e'
    for i,j in enumerate(bias_data.tumor.values):
        new_tmb = (data[data['tumor']==j]['uni_in_new'].values[0]+data[data['tumor']==j]['overlap'].values[0]-data[data['tumor']==j]['ori_in_overlap'].values[0]-data[data['tumor']==j]['ori_in_uni_new'].values[0])/50
        old_tmb = (data[data['tumor']==j]['uni_in_old'].values[0]+data[data['tumor']==j]['overlap'].values[0])/50
        record_tmb = tmb_discrepancy_all[tmb_discrepancy_all['tumor']==j]['old_tmb'].values[0]
        diff = new_tmb - old_tmb
        # if abs(diff) > 100:
        if abs(diff/new_tmb) > 0.2:
            if max(new_tmb+10,old_tmb+10) < 800:
                ax2.annotate(j,xy=(max(new_tmb+10,old_tmb+10),k))
            elif max(new_tmb+10,old_tmb+10) >= 800:
                ax2.annotate(j,xy=(700,k))
        ax2.plot([0,new_tmb],[k,k],color=color,linewidth=0.75,ls='-',label='new')
        ax2.plot([0,old_tmb],[k+0.5,k+0.5],color=color,linewidth=0.75,ls='--',label='old')
        ax2.scatter(record_tmb,k+0.5,color='black',s=3,alpha=0.5)
        k += 1.5
    ax2.set_xlim([1,800])
    ax2.set_ylim([-1,k])
    ax2.set_xlabel('TMB')
    ax2.set_ylabel('Tumor#')
    ax2.set_yticks([])
    ax2.set_title(title)
    plt.legend(['New TMB','Old TMB','Record TMB'])
    fig.savefig(OUTPUTPATH + f'{title}_bia.png',dpi=300)

###
plotTmbDiscrepancy(overlap_data_norm,'tmb_discrepancy_based_on_norm')
plotTmbDiscrepancy(overlap_data_noob_hd5001,'tmb_discrepancy_based_on_noob_hd5001')
plotTmbDiscrepancy(overlap_data_noob_hd2005,'tmb_discrepancy_based_on_noob_hd2005')
plotTmbDiscrepancy(overlap_data_noob_hd2005both,'tmb_discrepancy_based_on_noob_hd2005both')
plotTmbDiscrepancy(overlap_data_noob_hd5010both,'tmb_discrepancy_based_on_noob_hd5010both')
plotTmbDiscrepancy(overlap_data_noob_hd2010both,'tmb_discrepancy_based_on_noob_hd2010both')
plotTmbDiscrepancy(overlap_data_varscan,'tmb_discrepancy_based_on_varscan')
###

def coveragePlot(data,title):
    def getCoverage(data):
        new_coverage_list = []
        old_coverage_list = []
        for i in range(len(data)):
            new_coverage = (data.iloc[i]['overlap']-data.iloc[i]['ori_in_overlap'])/(data.iloc[i]['overlap']-data.iloc[i]['ori_in_overlap']+data.iloc[i]['uni_in_new']-data.iloc[i]['ori_in_uni_new'])
            old_coverage = (data.iloc[i]['overlap']-data.iloc[i]['ori_in_overlap'])/(data.iloc[i]['overlap']-data.iloc[i]['ori_in_overlap']+data.iloc[i]['uni_in_old'])
            new_coverage_list.append(new_coverage)
            old_coverage_list.append(old_coverage)
        return new_coverage_list,old_coverage_list
    new_coverage_list,old_coverage_list = getCoverage(data)
    fig, ax = plt.subplots(figsize=(6,6))
    sns.scatterplot(new_coverage_list,old_coverage_list,hue=data['bias_status'],ax=ax,color='black',s=15)
    ax.vlines(0.75,0.75,1,linestyles='dashed')
    ax.hlines(0.75,0.75,1,linestyles='dashed')
    # k=0
    # for i,j in zip(new_coverage_list,old_coverage_list):
    #     if i < 0.75 or j < 0.75:
    #         ax.annotate(data.iloc[k]['tumor'],xy=(i,j))
    #     k += 1
    ax.set_xlabel('New Coverage')
    ax.set_ylabel('Old Coverage')
    ax.set_title(title)
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    fig.savefig(OUTPUTPATH + f'{title}.png',dpi=300)

coveragePlot(overlap_data_norm,'coverage_based_on_norm')
coveragePlot(overlap_data_noob_hd5001,'coverage_based_on_noob_hd5001')
coveragePlot(overlap_data_noob_hd2005,'coverage_based_on_noob_hd2005')
coveragePlot(overlap_data_noob_hd2005both,'coverage_based_on_noob_hd2005both')
coveragePlot(overlap_data_noob_hd5010both,'coverage_based_on_noob_hd5010both')
coveragePlot(overlap_data_noob_hd2010both,'coverage_based_on_noob_hd2010both')
coveragePlot(overlap_data_varscan,'coverage_based_on_varscan')

