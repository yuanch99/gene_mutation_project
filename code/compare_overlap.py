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

#load cluster labels
biased_pairs_hd = pd.read_csv(DATAPATH + 'biased_pairs_hd.csv',index_col=0)
unbiased_pairs_hd = pd.read_csv(DATAPATH + 'unbiased_pairs_hd.csv',index_col=0)
cross_cluster_hd = pd.read_csv(DATAPATH + 'cross_cluster_hd.csv',index_col=0)

orientation_data = pd.read_csv(DATAPATH + 'mean_std_log2_all_samples_complete.csv',index_col=0)
tmb_discrepancy_all = pd.read_csv(DATAPATH + 'tmb_discrepancy_all.csv',index_col=0)

overlap_data = pd.read_csv(DATAPATH + 'overlap_counts_norm.csv')
overlap_data_noob = pd.read_csv(DATAPATH + 'overlap_counts_no-ob.csv')
tmb_recal = pd.read_csv(DATAPATH + 'tmb_re-cal_norm.csv')
tmb_recal_noob = pd.read_csv(DATAPATH + 'tmb_re-cal_no-ob.csv')

# compare old and old recalulated tmb
# old_recal - old
figs1,axes1 = plt.subplots(figsize=(16,8))
ax1 = plt.subplot(1,2,1)
diff_total = []

diff_biased1 = []
diff_unbiased1 = []
diff_cross1 = []
for k,l in enumerate(tmb_recal['tumor'].tolist()):
    diff = (tmb_recal['old_tmb_recal'].values[k]-tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['old_tmb'].values[0])
    diff_total.append(diff)
    if tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] == 'biased': 
        color = '#e4455e'
        diff_biased1.append(diff)
    elif tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] == 'unbiased': 
        color = '#86a1a9'
        diff_unbiased1.append(diff)
    elif tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] == 'cross_cluster': 
        color = '#b18f6a'
        diff_cross1.append(diff)
    if abs(diff) > 100:
        ax1.annotate(l,xy=(min(abs(diff),300),k))
    plt.plot([0,diff],[k,k],color=color)
ax1.set_title(f'tmb_diff_oldrecal_old, midian=%.2f,mean=%.2f'%(np.median(np.abs(diff_total)),np.mean(np.abs(diff_total))))
ax1.set_xlabel(f'tmb(old_re-cal - old)')
ax1.set_ylabel('tumor')
ax1.set_ylim(0,len(tmb_recal['tumor'].tolist()))
ax1.set_xlim(-400,400)
bias_line = mpatches.Patch(color='#e4455e',label=f'biased,,midian=%.2f,mean=%.2f'%(
    np.median(np.abs(diff_biased1)),np.mean(np.abs(diff_biased1))
))
unbiased_line = mpatches.Patch(color='#86a1a9',label=f'unbiased,,midian=%.2f,mean=%.2f'%(
    np.median(np.abs(diff_unbiased1)),np.mean(np.abs(diff_unbiased1))
))

cross_cluster_line = mpatches.Patch(color='#b18f6a',label=f'cross_cluster,,midian=%.2f,mean=%.2f'%(
    np.median(np.abs(diff_cross1)),np.mean(np.abs(diff_cross1))
))
ax1.legend(handles=[bias_line,unbiased_line,cross_cluster_line],bbox_to_anchor=(0.9, -0.1))


ax2=plt.subplot(1,2,2)
diff_total = []
diff_biased2 = []
diff_unbiased2 = []
diff_cross2 = []
# new - old_recal
for k,l in enumerate(tmb_recal['tumor'].tolist()):
    diff = tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['tmb_discrepancy_OB_threshold'].values[0]-tmb_recal['old_tmb_recal'].values[k]
    diff_total.append(diff)
    if tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] == 'biased': 
        color = '#e4455e'
        diff_biased2.append(diff)
    elif tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] == 'unbiased': 
        color = '#86a1a9'
        diff_unbiased2.append(diff)
    elif tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] == 'cross_cluster': 
        color = '#b18f6a'
        diff_cross2.append(diff)
    if abs(diff) > 100:
        ax2.annotate(l,xy=(min(abs(diff),300),k))
    plt.plot([0,diff],[k,k],color=color)
ax2.set_title(f'tmb_diff_new_oldrecal, midian=%.2f,mean=%.2f,n=n=%d'%(np.median(np.abs(diff_total)),np.mean(np.abs(diff_total)),len(diff_total)))
ax2.set_xlabel(f'tmb(new - old_re-cal)')
ax2.set_ylabel('tumor')
ax2.set_ylim(0,len(tmb_recal['tumor'].tolist()))
ax2.set_xlim(-400,400)

bias_line = mpatches.Patch(color='#e4455e',label=f'biased,,midian=%.2f,mean=%.2f'%(
    np.median(np.abs(diff_biased2)),np.mean(np.abs(diff_biased2))
))
unbiased_line = mpatches.Patch(color='#86a1a9',label=f'unbiased,,midian=%.2f,mean=%.2f'%(
    np.median(np.abs(diff_unbiased2)),np.mean(np.abs(diff_unbiased2))
))
cross_cluster_line = mpatches.Patch(color='#b18f6a',label=f'cross_cluster,,midian=%.2f,mean=%.2f'%(
    np.median(np.abs(diff_cross2)),np.mean(np.abs(diff_cross2))
))
ax2.legend(handles=[bias_line,unbiased_line,cross_cluster_line],bbox_to_anchor=(0.9, -0.1))

figs1.savefig(OUTPUTPATH + 'tmb_diff_new_oldrecal.png',dpi=300,bbox_inches="tight")

def overlap(overlap_data,type,tmb_discrepancy_all=tmb_discrepancy_all):
    # overlap venn3_diagram
    for i,j in enumerate(overlap_data['tumor'].tolist()):
        bias_status = tmb_discrepancy_all[tmb_discrepancy_all['tumor']==j]['bias'].values[0]
        venn3(
            set_labels = ('old_PASS','new_PASS','orientation_filtered'),
            subsets = (
                overlap_data.iloc[i]['uni_in_old'],
                (overlap_data.iloc[i]['uni_in_new']-overlap_data.iloc[i]['ori_in_uni_new']),
                (overlap_data.iloc[i]['overlap']-overlap_data.iloc[i]['ori_in_overlap']),
                0,
                0,
                overlap_data.iloc[i]['ori_in_uni_new'],
                overlap_data.iloc[i]['ori_in_overlap'],
            )
        )
        bia=tmb_discrepancy_all[tmb_discrepancy_all['tumor']==j]['bias'].values[0]
        plt.title(f'{j}_{bia}')
        plt.savefig(OUTPUTPATH + f'{type}/venn_diagram/venn3_{j}.png',dpi=300,bbox_inches="tight")
        plt.close()

    # overlap venn2_diagram
    overlap_coverage = pd.DataFrame(columns=['tumor','normal','old_cover','new_cover','bias'])
    for i,j in enumerate(overlap_data['tumor'].tolist()):
        set_labels = ('old','new')
        #turn off OB for biased and cross_cluster
        if tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] in ['biased','cross_cluster']: 
            set=(
                overlap_data.iloc[i]['uni_in_old'],
                overlap_data.iloc[i]['uni_in_new'],
                overlap_data.iloc[i]['overlap'],
            )
        #turn on OB for unbiased
        elif tmb_discrepancy_all[tmb_discrepancy_all['tumor']==l]['bias'].values[0] == ['unbiased']: 
            set=(
                overlap_data.iloc[i]['uni_in_old'],
                (overlap_data.iloc[i]['uni_in_new']-overlap_data.iloc[i]['ori_in_uni_new']),
                (overlap_data.iloc[i]['overlap']-overlap_data.iloc[i]['ori_in_overlap']),
            )
        venn2(
            set_labels = set_labels,
            subsets = set,
        )
        bia=tmb_discrepancy_all[tmb_discrepancy_all['tumor']==j]['bias'].values[0]
        old_cover=set[2]/(set[0]+set[2])
        new_cover=set[2]/(set[1]+set[2])
        plt.title(f'{j}_{bia},old_cover=%.2f,new_cover=%.2f'%(old_cover,new_cover))
        
        overlap_coverage.loc[i] = [j,overlap_data['tumor'].tolist()[i],old_cover,new_cover,bia]
        
        plt.savefig(OUTPUTPATH + f'{type}/venn_diagram/venn2_{j}.png',dpi=300,bbox_inches="tight")
        plt.close()

    return overlap_coverage

overlap_coverage = overlap(overlap_data=overlap_data,type='norm',tmb_discrepancy_all=tmb_discrepancy_all)
overlap_coverage_noob = overlap(overlap_data=overlap_data_noob,type='no-ob_hd',tmb_discrepancy_all=tmb_discrepancy_all)

#coverage plot
sns.scatterplot(x='old_cover',y='new_cover',hue='bias',data=overlap_coverage_noob,)
k=0
for i in overlap_coverage_noob.index:
    if overlap_coverage_noob.loc[i]['old_cover'] < 0.8 or overlap_coverage_noob.loc[i]['new_cover'] < 0.8:
        plt.annotate(overlap_coverage_noob.loc[i]['tumor'],xy=(overlap_coverage_noob.loc[i]['old_cover'],overlap_coverage_noob.loc[i]['new_cover']),xytext=(overlap_coverage.loc[i]['old_cover']+0.1,overlap_coverage.loc[i]['new_cover']+0.1),arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))
        k += 0.05
plt.vlines(0.8,0.8,1,linestyles='dashed')
plt.hlines(0.8,0.8,1,linestyles='dashed')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('old_cover')
plt.ylabel('new_cover')
plt.title('old_cover vs new_cover, n=%d'%len(overlap_coverage_noob))
plt.legend()
plt.savefig(OUTPUTPATH + 'no-ob_hd/old_cover_vs_new_cover.png',dpi=300,bbox_inches="tight")

##POLD1 and POLE gene overlap
d = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M',
     }
#flip dict
d_flip = dict((v,k) for k,v in d.items())

new_pol_mut = pd.read_csv(DATAPATH + 'pol_anno.csv',header=None)
pol_overlap = pd.DataFrame(columns=['tumor','pole','pold1'])

for j,i in enumerate(new_pol_mut[0].unique()):
    pole_muts_new = new_pol_mut[(new_pol_mut[0]==i) & (new_pol_mut[2]=='POLE')][3].tolist()
    pold1_muts_new = new_pol_mut[(new_pol_mut[0]==i) & (new_pol_mut[2]=='POLD1')][3].tolist()

    pole_muts_old = hyper_data[hyper_data['TumorName']==i]['POLEcell'].fillna('none').values[0]
    if pole_muts_old == 'none':
        pole_muts_old = 'none'
    else:
        pole_muts_old.split()[0]
        pole_muts_old = pole_muts_old.split('.')[1]
        pole_muts_old = pole_muts_old.split(' ')[0]
        #replace the letters in the flipped dictionary
        pole_muts_old = ''.join([d_flip[j] if j in d_flip else j for j in pole_muts_old])

    pold1_muts_old = hyper_data[hyper_data['TumorName']==i]['POLDcell'].fillna('none').values[0]
    if pold1_muts_old == 'none':
        pold1_muts_old = 'none'
    else:
        pold1_muts_old.split()[0]
        pold1_muts_old = pold1_muts_old.split('.')[1]
        pold1_muts_old = pold1_muts_old.split(' ')[0]
        #replace the letters in the flipped dictionary
        pold1_muts_old = ''.join([d_flip[j] if j in d_flip else j for j in pold1_muts_old])
    #get the overlap
    pole_status=0
    pold1_status=0

    if pole_muts_old == 'none':
        if len(pole_muts_new) == 0:
            pole_status = 1 #no old, no new
        else:
            pole_status = 1 #no old, yes new
    else:
        if pole_muts_old in pole_muts_new:
            if len(pole_muts_new) == 1:
                pole_status = 1 #1 old, 1 new
            else:
                pole_status = 1 #1 old, more than 1 new
        else:
            pole_status = 0 #1 old, no new

    if pold1_muts_old == 'none':
        if len(pold1_muts_new) == 0:
            pold1_status = 1 #no old, no new
        else:
            pold1_status = 1 #no old, yes new
    else:
        if pold1_muts_old in pold1_muts_new:
            if len(pold1_muts_new) == 1:
                pold1_status = 1 #1 old, 1 new
            else:
                pold1_status = 1 #1 old, more than 1 new
        else:
            pold1_status = 0 #1 old, no new

    pol_overlap.loc[j] = [i,pole_status,pold1_status]

#heat map
plt.figure(figsize=(5,13))
sns.set_style('darkgrid')
sns.heatmap(np.array(pol_overlap.iloc[:,1:]).astype(int),cbar=False,annot=False,fmt='d',cmap='Blues')
#set y tick labels
plt.yticks(np.arange(0,len(pol_overlap.index),1),pol_overlap.tumor,rotation=0)
plt.xticks(np.arange(0,2,0.5),['','POLE','','POLD1'])
plt.title('POL gene overlap')
plt.savefig(OUTPUTPATH + 'pol_overlap.png',dpi=300,bbox_inches="tight")

#MMR gene overlap
