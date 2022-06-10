from calendar import c
from lib2to3.pytree import HUGE
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os.path
import seaborn as sns
DATAPATH = '/Users/yuanchang/Documents/MBP/TaboriLab/codes/mutation_project/vcf_validation/data/'
OUTPUTPATH = '/Users/yuanchang/Documents/MBP/TaboriLab/codes/mutation_project/vcf_validation/output/'

all_new_filters = pd.read_csv(DATAPATH+'new_all_filters_count.csv',header=None)
samples = np.unique(all_new_filters[0])
good_sample = ["MD1460T1","MD1440T2","MD1515T1","MD1564T1","4","6","MD1456T2","MD1630T05","MD1468T1","MD1630T02","JG"]
bad_sample = ["D1807_1","D1807_2","MD1468T2","D1577","MD1566T2","D132","MD1385T2_2","D1806","MD1454T1","MD1456T3","D1120","D1804_3","D1144","MD1475T1","MD1385T2_1","MD0134T10","D1121","X3","D1243","D1804_2","D1423","X2","MD1293T1","D1410","MD1273T25","MD0134T11","MD190T2","MD101T10","MD1341T3","MD1273T26","D1764","M1273T2","MD0134T12","MD1303T1"]



#######draw stack bar plot with filter PASS and orientation and (orientation and others) and others filters
def passOrientationAndOtherSum(df):
    n_pass = df[df[2]=='PASS'].iloc[0,1]
    n_orientation = df[df[2]=='orientation'][1].sum()
    n_orientation_and_others = df[df[2].str.contains('orientation')][1].sum()-n_orientation
    n_others = df[1].sum()-n_pass-n_orientation-n_orientation_and_others
    return n_pass.astype(int),n_orientation.astype(int),n_orientation_and_others.astype(int),n_others.astype(int)

stack_bar_plot_good = pd.DataFrame(columns=['PASS','orientation','orientation+others','others'])
stack_bar_plot_bad = pd.DataFrame(columns=['PASS','orientation','orientation+others','others'])
for sample in good_sample:
        stack_bar_plot_good.loc[sample] = passOrientationAndOtherSum(all_new_filters[all_new_filters[0]==sample])
for sample in bad_sample:
        stack_bar_plot_bad.loc[sample] = passOrientationAndOtherSum(all_new_filters[all_new_filters[0]==sample])

tmb_good = pd.read_csv(DATAPATH+'tmb_good.csv')
tmb_bad = pd.read_csv(DATAPATH+'tmb_bad.csv')


fig, ax = plt.subplots(figsize=(13,10))
ax1 = plt.subplot(2,1,1)
ax1.set_title('Good samples')
ax1.set_ylabel('Number of variants')
ax1.set_xlabel('Sample')
ax1.set_xticklabels(good_sample,rotation=45)
ax1.set_xticks(np.arange(len(good_sample)))
ax1.bar(np.arange(len(good_sample)),stack_bar_plot_good['PASS'],color='#08cad1',label='PASS')
ax1.bar(np.arange(len(good_sample)),stack_bar_plot_good['orientation'],bottom=stack_bar_plot_good['PASS'],color='#ff6961',label='orientation')
ax1.bar(np.arange(len(good_sample)),stack_bar_plot_good['orientation+others'],bottom=stack_bar_plot_good['PASS']+stack_bar_plot_good['orientation'],color='#ffb480',label='orientation and other filters')
ax1.bar(np.arange(len(good_sample)),stack_bar_plot_good['others'],bottom=stack_bar_plot_good['PASS']+stack_bar_plot_good['orientation']+stack_bar_plot_good['orientation+others'],color='#9d94ff',label='others')
plt.legend(loc='upper right')
ax2 = ax1.twinx()
ax2.set_ylabel('esitimate TMB')
ax2.scatter(np.arange(len(good_sample)),tmb_good['old TMB'],color='black',label='old TMB')
ax2.scatter(np.arange(len(good_sample)),tmb_good['new TMB'],color='g',label='new TMB')
ax2.set_ylim(ax1.get_ylim()[0]/50,ax1.get_ylim()[1]/50)
ax2.legend(loc='upper center')


ax3 = plt.subplot(2,1,2)
ax3.set_title('Bad samples')
ax3.set_ylabel('Number of variants')
ax3.set_xlabel('Sample')
ax3.set_xticklabels(bad_sample,rotation=90)
ax3.set_xticks(np.arange(len(bad_sample)))
ax3.bar(np.arange(len(bad_sample)),stack_bar_plot_bad['PASS'],color='#08cad1',label='PASS')
ax3.bar(np.arange(len(bad_sample)),stack_bar_plot_bad['orientation'],bottom=stack_bar_plot_bad['PASS'],color='#ff6961',label='orientation')
ax3.bar(np.arange(len(bad_sample)),stack_bar_plot_bad['orientation+others'],bottom=stack_bar_plot_bad['PASS']+stack_bar_plot_bad['orientation'],color='#ffb480',label='orientation and other filters')
ax3.bar(np.arange(len(bad_sample)),stack_bar_plot_bad['others'],bottom=stack_bar_plot_bad['PASS']+stack_bar_plot_bad['orientation']+stack_bar_plot_bad['orientation+others'],color='#9d94ff',label='others')
ax3.set_ylim(ax1.get_ylim())
plt.legend(loc='upper right')
ax4 = ax3.twinx()
ax4.set_ylabel('esitimate TMB')
ax4.scatter(np.arange(len(bad_sample)),tmb_bad['old TMB'],color='black',label='old TMB')
ax4.scatter(np.arange(len(bad_sample)),tmb_bad['new TMB'],color='g',label='new TMB')
ax4.set_ylim(ax1.get_ylim()[0]/50,ax1.get_ylim()[1]/50)
ax4.legend(loc='upper center')

fig.savefig(OUTPUTPATH+'stack_bar_plot_good_bad.png',dpi=300)

###########

#######draw pie plot with filter PASS and orientation and (orientation and others) and others filters
example_new = pd.read_csv(DATAPATH+'example_new_filter_count.csv',header=None)
example_old = pd.read_csv(DATAPATH+'example_old_filter_count.csv',header=None)

def piePlotData(df,sample,version):
    if version=='new':
        filters = ['PASS','orientation','orientation and other filter','others']
        values = [
        df[(df[0]==sample) & (df[2]==filters[0])][1].sum(),
        df[(df[0]==sample) & (df[2]==filters[1])][1].sum(),
        df[(df[0]==sample) & (df[2].str.contains(filters[1]))][1].sum(),
        df[(df[0]==sample)][1].sum()-df[(df[0]==sample) & (df[2].str.contains(filters[0]))][1].sum().astype(int)-df[(df[0]==sample) & (df[2].str.contains(filters[1]))][1].sum().sum()
        ]
    else:
        filters = [df[(df[0]==sample)].iloc[0,2],
            df[(df[0]==sample)].iloc[1,2],
            df[(df[0]==sample)].iloc[2,2],
            df[(df[0]==sample)].iloc[3,2],
            df[(df[0]==sample)].iloc[4,2]]
        values = [df[(df[0]==sample)].iloc[0,1],
            df[(df[0]==sample)].iloc[1,1],
            df[(df[0]==sample)].iloc[2,1],
            df[(df[0]==sample)].iloc[3,1],
            df[(df[0]==sample)].iloc[4,1]]
    return filters,values

def pieChart(ax,df,sample,version):
    filters,values = piePlotData(df,sample,version)
    ax.pie(values,labels=filters,autopct='%1.2f%%' ,startangle=90,
        colors=['#42d6a4','#ff6961','#ffb480','#9d94ff','#ffb480','#59adf6'])
    ax.axis('equal')
    ax.set_title(f'{sample} {version}')

fig2,axes2 = plt.subplots(2,2,figsize=(10,10))
axes2[0,0] = pieChart(axes2[0,0],example_new,sample='JG',version='new')
axes2[0,1] = pieChart(axes2[0,1],example_old,sample='JG',version='old')
axes2[1,0] = pieChart(axes2[1,0],example_new,sample='MD1341T3',version='new')
axes2[1,1] = pieChart(axes2[1,1],example_old,sample='MD1341T3',version='old')
fig2.savefig(OUTPUTPATH+'pie_chart_example_new_old.png',dpi=300)

####vcf histogram
JG_vaf = pd.read_csv(DATAPATH+'JG_vaf.csv',header=None,on_bad_lines='skip')
MD1341T3_vaf = pd.read_csv(DATAPATH+'MD1341T3_vaf.csv',header=None,on_bad_lines='skip')
fig3,axes3 = plt.subplots(2,2,figsize=(10,10))
axes3[0,0].hist(JG_vaf[JG_vaf[0]=='PASS'][1],bins=100,color='#42d6a4')
axes3[0,0].set(xlabel='VAF',ylabel='count')
axes3[0,0].set_title('JG_PASS')
axes3[0,1].hist(JG_vaf[JG_vaf[0]!='PASS'][1],bins=100,color='#42d6a4')
axes3[0,1].set(xlabel='VAF',ylabel='count')
axes3[0,1].set_title('JG_filtered')
axes3[1,0].hist(MD1341T3_vaf[MD1341T3_vaf[0].isin(['PASS','orientation'])][1],bins=100,color='#ff6961')
axes3[1,0].set(xlabel='VAF',ylabel='count')
axes3[1,0].set_title('MD1341T3_PASS+orientation')
axes3[1,1].hist(MD1341T3_vaf[~(MD1341T3_vaf[0].isin(['PASS','orientation']))][1],bins=100,color='#ff6961')
axes3[1,1].set(xlabel='VAF',ylabel='count')
axes3[1,1].set_title('MD1341T3_other filtered')
fig3.savefig(OUTPUTPATH+'vcf_histogram.png',dpi=300)

####orientation histogram
vcf_orientation_count = pd.read_csv(DATAPATH+'vcf_orientation_count.csv',on_bad_lines='skip')

fig4,axes4 = plt.subplots(2,2,figsize=(10,10))
axes4[0,0] = sns.scatterplot(
    x=vcf_orientation_count[vcf_orientation_count['sample'].isin(good_sample)]['T_F1R2_REF']+vcf_orientation_count[vcf_orientation_count['sample'].isin(good_sample)]['T_F1R2_ALT'],
    y=vcf_orientation_count[vcf_orientation_count['sample'].isin(good_sample)]['T_F2R1_REF']+vcf_orientation_count[vcf_orientation_count['sample'].isin(good_sample)]['T_F2R1_ALT'],
    marker='.',
    hue=vcf_orientation_count[vcf_orientation_count['sample'].isin(good_sample)]['sample'],
    ax=axes4[0,0],
    legend = False
    )

#draw a 45 degree line to show the orientation
axes4[0,0].plot(axes4[0,0].get_xlim(),axes4[0,0].get_ylim(),'k--')
axes4[0,0].set(xlabel='F1R2',ylabel='F2R1',title='good_sample')
axes4[0,1] = sns.scatterplot(
    x=vcf_orientation_count[vcf_orientation_count['sample'].isin(bad_sample)]['T_F1R2_REF']+vcf_orientation_count[vcf_orientation_count['sample'].isin(bad_sample)]['T_F1R2_ALT'],
    y=vcf_orientation_count[vcf_orientation_count['sample'].isin(bad_sample)]['T_F2R1_REF']+vcf_orientation_count[vcf_orientation_count['sample'].isin(bad_sample)]['T_F2R1_ALT'],
    marker='.',
    hue=vcf_orientation_count[vcf_orientation_count['sample'].isin(bad_sample)]['sample'],
    ax=axes4[0,1],
    legend = False
    )

axes4[0,1].set(xlabel='F1R2',ylabel='F2R1',title='bad_sample')
axes4[0,1].plot(axes4[0,1].get_xlim(),axes4[0,1].get_ylim(),'k--')
axes4[0,1].plot([0,100],[0,100],'k--')
#plot ration
axes4[1,0] = sns.scatterplot(
    x=vcf_orientation_count[(vcf_orientation_count['sample'].isin(good_sample) & (vcf_orientation_count['filter']=='PASS'))]['T_F1R2_REF']/vcf_orientation_count[(vcf_orientation_count['sample'].isin(good_sample) & (vcf_orientation_count['filter']=='PASS'))]['T_F1R2_ALT'],
    y=vcf_orientation_count[(vcf_orientation_count['sample'].isin(good_sample) & (vcf_orientation_count['filter']=='PASS'))]['T_F2R1_REF']/vcf_orientation_count[(vcf_orientation_count['sample'].isin(good_sample) & (vcf_orientation_count['filter']=='PASS'))]['T_F2R1_ALT'],
    marker='.',
    hue=vcf_orientation_count[vcf_orientation_count['sample'].isin(good_sample)]['sample'],
    ax=axes4[1,0],
    legend = False
    )
axes4[1,0].plot(axes4[1,0].get_xlim(),axes4[1,0].get_ylim(),'k--')
axes4[1,0].set(xlabel='F1R2 REF/ALT',ylabel='F2R1 REF/ALT',title='good_sample (PASS)')

axes4[1,1] = sns.scatterplot(
    x=vcf_orientation_count[(vcf_orientation_count['sample'].isin(bad_sample) & (vcf_orientation_count['filter'].isin(['PASS','orientation'])))]['T_F1R2_REF']/vcf_orientation_count[(vcf_orientation_count['sample'].isin(bad_sample) & (vcf_orientation_count['filter'].isin(['PASS','orientation'])))]['T_F1R2_ALT'],
    y=vcf_orientation_count[(vcf_orientation_count['sample'].isin(bad_sample) & (vcf_orientation_count['filter'].isin(['PASS','orientation'])))]['T_F2R1_REF']/vcf_orientation_count[(vcf_orientation_count['sample'].isin(bad_sample) & (vcf_orientation_count['filter'].isin(['PASS','orientation'])))]['T_F2R1_ALT'],
    marker='.',
    hue=vcf_orientation_count[vcf_orientation_count['sample'].isin(bad_sample)]['sample'],
    ax=axes4[1,1],
    legend = False
    )
axes4[1,1].plot(axes4[1,1].get_xlim(),axes4[1,1].get_ylim(),'k--')
axes4[1,1].set(xlabel='F1R2 REF/ALT',ylabel='F2R1 REF/ALT',title='bad_sample (PASS & orientation)',xlim=(-5,100),ylim=(-5,300))

fig4.savefig(OUTPUTPATH+'vcf_orientation_count.png',dpi=300)






