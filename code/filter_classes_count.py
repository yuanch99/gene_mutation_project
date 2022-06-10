import pandas as pd
from pysam import VariantFile
import matplotlib.pyplot as plt
import numpy as np
import os.path
DATAPATH = '/Users/yuanchang/Documents/MBP/TaboriLab/codes/mutation_project/vcf_validation/data/'
OUTPUTPATH = '/Users/yuanchang/Documents/MBP/TaboriLab/codes/mutation_project/vcf_validation/output/'



def getFilters(data,sample_version):
  filter_file_path = DATAPATH+sample_version+'filters.csv'

  if os.path.exists(filter_file_path):
    filter_new=pd.read_csv(filter_file_path)

  else:
    filter_new = pd.DataFrame(columns=['ale','alt','chr','loc','filters'])
    for rec in data.fetch():
      filter_new = filter_new.append(
        {
          'ale':rec.alleles[0],
          'alt':rec.alleles[1],
          'chr':rec.contig,
          'loc':rec.pos,
          'filters':rec.filter.keys()
        },ignore_index=True
      )
      filter_new['filters_combine'] = filter_new['filters'].apply(lambda x: ','.join(map(str, x)))
    filter_new.to_csv((DATAPATH+sample_version+'filters.csv'))

  print('load filters done!')
  return filter_new

def filterStackBarPlot(filter_new,filter_old,sample_name):
  bar_plot_data = pd.DataFrame(columns=['sample','pass','orientation','others'])
  bar_plot_data = bar_plot_data.append(
      {
      'sample':'new',
      'pass':np.count_nonzero(filter_new['filters_combine']=='PASS'),
      'orientation':np.count_nonzero(filter_new['filters_combine']=='orientation'),
      'others':(len(filter_new)-np.count_nonzero(filter_new['filters_combine']=='PASS')-np.count_nonzero(filter_new['filters_combine']=='orientation'))
      },ignore_index=True
    )
  bar_plot_data = bar_plot_data.append(
      {
      'sample':'old',
      'pass':np.count_nonzero(filter_old['filters_combine']=='PASS'),
      'orientation':np.count_nonzero(filter_old['filters_combine']=='orientation'),
      'others':(len(filter_old)-np.count_nonzero(filter_old['filters_combine']=='PASS')-np.count_nonzero(filter_old['filters_combine']=='orientation'))
      },ignore_index=True
    )
  fig, ax = plt.subplots()
  ax.bar(bar_plot_data['sample'], bar_plot_data['pass'], label='PASS',color='#42d6a4')
  ax.bar(bar_plot_data['sample'], bar_plot_data['orientation'], bottom=bar_plot_data['pass'],label='Orientation',color='#ffb480')
  ax.bar(bar_plot_data['sample'], bar_plot_data['others'], bottom=bar_plot_data['pass']+bar_plot_data['orientation'],label='Others',color='#9d94ff')
  plt.legend(loc=9)
  ax.set_ylabel('Filters')
  ax.set_title(sample_name)
  fig.savefig((OUTPUTPATH+sample_name),dpi=300)
  return fig,ax


D1243_old = VariantFile((DATAPATH+'mutect2.D1243.MMR100B1.b37.merged.vcf.gz'))
D1243_new = VariantFile((DATAPATH+'D1243__MMR100B1.mutect2.filtered.wes.vcf.gz'))
D1243_new_filter = getFilters(D1243_new,'D1243_new')
D1243_old_filter = getFilters(D1243_old,'D1243_old')
filterStackBarPlot(D1243_new_filter,D1243_old_filter,'D1243')

JG_old = VariantFile((DATAPATH+'mutect2.JG.MMR152B1.b37.merged.vcf.gz'))
JG_new = VariantFile((DATAPATH+'JG__MMR152B1.mutect2.filtered.wes.vcf.gz'))
JG_new_filter = getFilters(JG_new,'JG_new')
JG_old_filter = getFilters(JG_old,'JG_old')
filterStackBarPlot(JG_new_filter,JG_old_filter,'JG')

MD1341T3_old = VariantFile((DATAPATH+'mutect2.MD1341T3.MD1341B1.with_PON.merged.20171110.vcf.gz'))
MD1341T3_new = VariantFile((DATAPATH+'MD1341T3__MD1341B1.mutect2.filtered.wes.vcf.gz'))
MD1341T3_new_filter = getFilters(MD1341T3_new,'MD1341T3_new')
MD1341T3_old_filter = getFilters(MD1341T3_old,'MD1341T3_old')
filterStackBarPlot(MD1341T3_new_filter,MD1341T3_old_filter,'MD1341T3')