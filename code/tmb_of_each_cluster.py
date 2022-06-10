from tkinter.ttk import Style
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.patches import  mlines

DATAPATH = '../data/'
OUTPUTPATH = '../output/'

#load cluster labels
biased_pairs_hd = pd.read_csv(DATAPATH + 'biased_pairs_hd.csv',index_col=0)
unbiased_pairs_hd = pd.read_csv(DATAPATH + 'unbiased_pairs_hd.csv',index_col=0)
cross_cluster_hd = pd.read_csv(DATAPATH + 'cross_cluster_hd.csv',index_col=0)
biased_pairs = pd.read_csv(DATAPATH + 'biased_pairs.csv',index_col=0)
unbiased_pairs = pd.read_csv(DATAPATH + 'unbiased_pairs.csv',index_col=0)
cross_cluster = pd.read_csv(DATAPATH + 'cross_cluster.csv',index_col=0)

#load old hyper data
hyper_data = pd.read_csv('/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/raw_data/old_hyperdata.csv')

#load new tmb
new_tmb = pd.read_csv(DATAPATH + 'coverage_and_tmb.csv')

#check if all the samples are in the hyper_data
for i,j in new_tmb.iloc[:,[0,1]].values:
    if i not in hyper_data['TumorName'].tolist():
        print(i,j)

orientation_data = pd.read_csv(DATAPATH + 'mean_std_log2_all_samples_complete.csv',index_col=0)

# tmb discrepencies between old and new (norm)
tmb_discrepancy = []
for i in orientation_data['sample'].tolist():
    if i in hyper_data['TumorName'].tolist():
        tmb_old = hyper_data[hyper_data['TumorName'] == i]['Mut/MB'].values[0]
        tmb_new = new_tmb[new_tmb['tumor'] == i]['tmb_snvs'].values[0]
    else:
        paired_tumor = new_tmb[new_tmb['normal'] == i]['tumor'].values[0]
        tmb_old = hyper_data[hyper_data['TumorName'] == paired_tumor]['Mut/MB'].values[0]
        tmb_new = new_tmb[new_tmb['tumor'] == paired_tumor]['tmb_snvs'].values[0]
    tmb_discrepancy.append(abs(tmb_new - tmb_old))

figs1,axes1 = plt.subplots(figsize=(10,9))
ax1 = plt.subplot(111)
plt.scatter(
    orientation_data['mean_log2_ratio_F1R2_over_F2R1_50plus'], 
    orientation_data['std_log2_ratio_F1R2_over_F2R1_50plus'], 
    c=tmb_discrepancy,
    cmap='inferno_r',
    s=100*(tmb_discrepancy/np.max(tmb_discrepancy)).flatten(),
)
plt.colorbar()
plt.clim(0,max(tmb_discrepancy))
ax1.set_xlabel('mean_log2_ratio_F1R2_over_F2R1_50plus')
ax1.set_ylabel('std_log2_ratio_F1R2_over_F2R1_50plus')
ax1.set_title('tmb_discrepancy_withOB_absolute')
# for i,j in enumerate(tmb_discrepancy):
#     if abs(j) > 300:
#         ax1.annotate(orientation_data['sample'].values[i],xy=(orientation_data['mean_log2_ratio_F1R2_over_F2R1_50plus'].values[i],orientation_data['std_log2_ratio_F1R2_over_F2R1_50plus'].values[i]))
figs1.savefig(OUTPUTPATH + 'tmb_discrepancy_withOB_absolute.png')


# plot the tmb change
tmb_norm = pd.read_csv(DATAPATH+'tmb_norm.csv',header=None)
tmb_norm['tmb']=tmb_norm[1]/(tmb_norm[2]/1000000)
tmb_norm.columns = ['tumor','snps','coverage','tmb']

tmb_no_ob = pd.read_csv(DATAPATH+'tmb_no-ob.csv',header=None)
tmb_no_ob['tmb']=tmb_no_ob[1]/(tmb_no_ob[2]/1000000)
tmb_no_ob.columns = ['tumor','snps','coverage','tmb']

tmb_discrepancy_all = pd.DataFrame(columns=[
    'tumor',
    'old_tmb',
    'tmb_discrepancy_withOB',
    'tmb_discrepancy_noOB',
    'tmb_discrepancy_OB_threshold',
    'bias'])
for i in biased_pairs_hd['tumor'].tolist():
    tmb_discrepancy_all.loc[len(tmb_discrepancy_all)] = [
        i,hyper_data[hyper_data['TumorName'] == i]['Mut/MB'].values[0],
        tmb_norm[tmb_norm['tumor'] == i]['tmb'].values[0],
        tmb_no_ob[tmb_no_ob['tumor'] == i]['tmb'].values[0],
        tmb_no_ob[tmb_no_ob['tumor'] == i]['tmb'].values[0],
        'biased']
for i in unbiased_pairs_hd['tumor'].tolist():
    tmb_discrepancy_all.loc[len(tmb_discrepancy_all)] = [
        i,hyper_data[hyper_data['TumorName'] == i]['Mut/MB'].values[0],
        tmb_norm[tmb_norm['tumor'] == i]['tmb'].values[0],
        tmb_no_ob[tmb_no_ob['tumor'] == i]['tmb'].values[0],
        tmb_norm[tmb_norm['tumor'] == i]['tmb'].values[0],
        'unbiased']
for i in cross_cluster_hd['tumor'].tolist():
    tmb_discrepancy_all.loc[len(tmb_discrepancy_all)] = [
        i,hyper_data[hyper_data['TumorName'] == i]['Mut/MB'].values[0],
        tmb_norm[tmb_norm['tumor'] == i]['tmb'].values[0],
        tmb_no_ob[tmb_no_ob['tumor'] == i]['tmb'].values[0],
        #turn off ob for cross cluster samples
        tmb_no_ob[tmb_no_ob['tumor'] == i]['tmb'].values[0],
        # tmb_norm[tmb_norm['tumor'] == i]['tmb'].values[0],
        'cross_cluster']


#plot diff between new and old tmb
figs2,axes2 = plt.subplots(figsize=(21,8))
for i,j in zip(range(1,4),['tmb_discrepancy_withOB','tmb_discrepancy_noOB','tmb_discrepancy_OB_threshold']):
    ax2 = plt.subplot(1,3,i)
    for k,l in enumerate(tmb_discrepancy_all['tumor'].tolist()):
        if tmb_discrepancy_all['bias'].values[k] == 'biased': color = '#e4455e'
        elif tmb_discrepancy_all['bias'].values[k] == 'unbiased': color = '#86a1a9'
        elif tmb_discrepancy_all['bias'].values[k] == 'cross_cluster': color = '#b18f6a'
        diff = (tmb_discrepancy_all[j].values[k]-tmb_discrepancy_all['old_tmb'].values[k])
        plt.plot([0,diff],[k,k],color=color)
        # if abs(diff)>100:
        #     if abs(diff)>400:
        #         ax2.annotate(tmb_discrepancy_all['tumor'].values[k],xy=(400,k))
        #     else:
        #         ax2.annotate(tmb_discrepancy_all['tumor'].values[k],xy=(diff,k))
    ax2.set_title(f'{j},midian=%.2f,mean=%.2f'%(np.median(abs(tmb_discrepancy_all[j].values-tmb_discrepancy_all['old_tmb'].values)),np.mean(abs(tmb_discrepancy_all[j].values-tmb_discrepancy_all['old_tmb'].values))))
    ax2.set_xlabel(f'tmb_discrepancy')
    ax2.set_ylabel('tumor ID')
    ax2.set_ylim(0,len(tmb_discrepancy_all))
    ax2.set_xlim(-max(tmb_discrepancy),max(tmb_discrepancy))
    bias_line = mpatches.Patch(color='#e4455e',label=f'biased,,midian=%.2f,mean=%.2f'%(
        np.median(abs(tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased']['old_tmb'].values)),
        np.mean(abs(tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased']['old_tmb'].values))))
    unbiased_line = mpatches.Patch(color='#86a1a9',label=f'unbiased,,midian=%.2f,mean=%.2f'%(
        np.median(abs(tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased']['old_tmb'].values)),
        np.mean(abs(tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased']['old_tmb'].values))))
    cross_cluster_line = mpatches.Patch(color='#b18f6a',label=f'cross_cluster,,midian=%.2f,mean=%.2f'%(
        np.median(abs(tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster']['old_tmb'].values)),
        np.mean(abs(tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster']['old_tmb'].values))))
    ax2.legend(handles=[bias_line,unbiased_line,cross_cluster_line],bbox_to_anchor=(0.9, -0.1))
figs2.savefig(OUTPUTPATH + 'tmb_discrepancy_OB_threshold_diff.png',dpi=300,bbox_inches="tight")

# a percentage version
figs3,axes3 = plt.subplots(figsize=(21,8))
for i,j in zip(range(1,4),['tmb_discrepancy_withOB','tmb_discrepancy_noOB','tmb_discrepancy_OB_threshold']):
    ax3 = plt.subplot(1,3,i)
    for k,l in enumerate(tmb_discrepancy_all['tumor'].tolist()):
        if tmb_discrepancy_all['bias'].values[k] == 'biased': color = '#e4455e'
        elif tmb_discrepancy_all['bias'].values[k] == 'unbiased': color = '#86a1a9'
        elif tmb_discrepancy_all['bias'].values[k] == 'cross_cluster': color = '#b18f6a'
        diff = (tmb_discrepancy_all[j].values[k]-tmb_discrepancy_all['old_tmb'].values[k])/tmb_discrepancy_all['old_tmb'].values[k]
        #check if FFPE
        if 'FFPE' in hyper_data[hyper_data['TumorName']==l]['TissueFormat'].values[0]:
            style = '--'
        else:
            style = '-'
        plt.plot([0,diff],[k,k],color=color,linestyle=style)
        # if abs(diff)>0.1:
        #     if abs(diff)>0.4:
        #         ax3.annotate(tmb_discrepancy_all['tumor'].values[k],xy=(0.4,k))
        #     else:
        #         ax3.annotate(tmb_discrepancy_all['tumor'].values[k],xy=(diff,k))
    ax3.set_title(f'{j},midian=%.2f,mean=%.2f'%(np.median(abs(tmb_discrepancy_all[j].values-tmb_discrepancy_all['old_tmb'].values))/tmb_discrepancy_all['old_tmb'].values.mean(),np.mean(abs(tmb_discrepancy_all[j].values-tmb_discrepancy_all['old_tmb'].values))/tmb_discrepancy_all['old_tmb'].values.mean()))
    ax3.set_xlabel(f'tmb_discrepancy_percentage')
    ax3.set_ylabel('tumor ID')
    ax3.set_ylim(0,len(tmb_discrepancy_all))
    ax3.set_xlim(-0.4,0.4)
    bias_line = mpatches.Patch(color='#e4455e',label=f'biased,,midian=%.2f,mean=%.2f'%(
        np.median(
            abs(
                (tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased']['old_tmb'].values)/tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased']['old_tmb'].values
                )),
        np.mean(
            abs(
                (tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased']['old_tmb'].values)/tmb_discrepancy_all[tmb_discrepancy_all['bias']=='biased']['old_tmb'].values
                ))))
    unbiased_line = mpatches.Patch(color='#86a1a9',label=f'unbiased,,midian=%.2f,mean=%.2f'%(
        np.median(
            abs(
                (tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased']['old_tmb'].values)/tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased']['old_tmb'].values
                )),
        np.mean(
            abs(
                (tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased']['old_tmb'].values)/tmb_discrepancy_all[tmb_discrepancy_all['bias']=='unbiased']['old_tmb'].values
                ))))
    cross_cluster_line = mpatches.Patch(color='#b18f6a',label=f'cross_cluster,,midian=%.2f,mean=%.2f'%(
        np.median(
            abs(
                (tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster']['old_tmb'].values)/tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster']['old_tmb'].values
                )),
        np.mean(
            abs(
                (tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster'][j].values-tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster']['old_tmb'].values)/tmb_discrepancy_all[tmb_discrepancy_all['bias']=='cross_cluster']['old_tmb'].values
                ))))
    dashline = mlines.Line2D([], [], color='#000000', linestyle='--', label='FFPE')
    ax3.legend(handles=[bias_line,unbiased_line,cross_cluster_line,dashline],bbox_to_anchor=(0.9, -0.1))
figs3.savefig(OUTPUTPATH+'tmb_discrepancy_percentage.png',dpi=300,bbox_inches='tight')

# tmb discrepencies between old and new (norm) after thresholding
tumors_normals = pd.read_csv(DATAPATH+'tumors_normals_finished_completely.csv')
tmb_discrepancy_threshold=[]
for i in orientation_data['sample'].values:
    if i in tmb_discrepancy_all['tumor'].values:
        tmb_discrepancy_threshold.append(
            tmb_discrepancy_all[tmb_discrepancy_all['tumor']==i]['tmb_discrepancy_OB_threshold'].values[0]-tmb_discrepancy_all[tmb_discrepancy_all['tumor']==i]['old_tmb'].values[0]
            )
    else:
        tumor=tumors_normals[tumors_normals['normal']==i]['tumor'].values[0]
        tmb_discrepancy_threshold.append(
            tmb_discrepancy_all[tmb_discrepancy_all['tumor']==tumor]['tmb_discrepancy_OB_threshold'].values[0]-tmb_discrepancy_all[tmb_discrepancy_all['tumor']==tumor]['old_tmb'].values[0]
            )
tmb_discrepancy_threshold = np.abs(tmb_discrepancy_threshold)

figs4,axes4 = plt.subplots(figsize=(10,9))
ax4 = plt.subplot(111)
plt.scatter(
    orientation_data['mean_log2_ratio_F1R2_over_F2R1_50plus'], 
    orientation_data['std_log2_ratio_F1R2_over_F2R1_50plus'], 
    c=tmb_discrepancy_threshold,
    cmap='inferno_r',
    s=100*(tmb_discrepancy_threshold/np.max(tmb_discrepancy_threshold)).flatten(),
)
# for i,j in enumerate(tmb_discrepancy_threshold):
#     if abs(j) > 300:
#         anno = f"{orientation_data['sample'].values[i]},   %.1f"%j
#         ax4.annotate(anno,xy=(orientation_data['mean_log2_ratio_F1R2_over_F2R1_50plus'].values[i],orientation_data['std_log2_ratio_F1R2_over_F2R1_50plus'].values[i]))
plt.colorbar()
plt.clim(0,max(tmb_discrepancy))
ax4.set_xlabel('mean_log2_ratio_F1R2_over_F2R1_50plus')
ax4.set_ylabel('std_log2_ratio_F1R2_over_F2R1_50plus')
ax4.set_title('tmb_discrepancy_after_threshold_OB_absolute')
figs4.savefig(OUTPUTPATH + 'tmb_discrepancy_after_threshold_absolute.png',dpi=300,bbox_inches='tight')

tmb_discrepancy_all.to_csv(DATAPATH+'tmb_discrepancy_all.csv')