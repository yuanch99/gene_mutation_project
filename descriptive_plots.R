library(tidyverse)
library(biomaRt)
library(ggplot2)
WRITE <- TRUE

#Load Data
dataPath <- '~/tabori/data/mutation/raw_data/santiago_pipeline/vcf_annotation/'
outPutPath <- "/Users/yuanchang/Documents/MBP/TaboriLab/codes/mutation_project/gene_mutation_project/output/"
initial <- read.csv(paste(dataPath,'all_18_Somatic_snpEff_annotations_20220301.csv',
                          sep=''),nrows=3000)
classes <- sapply(initial,class)
rm(initial)

somMutData <- read.csv(
  paste(dataPath,'all_18_Somatic_snpEff_annotations_20220301.csv',sep=''),
  colClasses=classes)
germMutData <- read.csv(
  paste(dataPath,'all_18_Germline_snpEff_annotations_20220301.csv',sep=''),
  colClasses=classes)

rm(classes)
gc()

#compare germline mutation counts in MMR and PR genes to total mutation counts

germMMRPDMutCount <- read.csv('/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/vanessa_analyses/analyses_Feb_26_2022/analyses/germ_MMR_PD_muts_count.csv',
                              header = FALSE,col.names = c('MMR_PR_count','sample'))

count_sample_muts <- function(germMMRPDMutCount,mutData){
  mutCount <- c()
  for (i in germMMRPDMutCount$sample){
    mutCount <- append(mutCount,length(mutData[mutData$Sample==i,1]))
  }
  germMMRPDMutCount$total_mut <- mutCount
  return(germMMRPDMutCount)
}

germSampleMutCount <- count_sample_muts(germMMRPDMutCount,germMutData)
somSampleMutCount <- count_sample_muts(germMMRPDMutCount,somMutData)

germMMRPDvsTotalPlot <- qplot(data=germSampleMutCount,
                              x=MMR_PR_count,y=total_mut,
                              color=sample,show.legend=FALSE)+
  xlab('Somatic_MMR_PR_Mut_count')+
  ggtitle('Germline MMRPD mut counts vs Total mut counts')
somMMRPDvsTotalPlot <- qplot(data=somSampleMutCount,
                              x=MMR_PR_count,y=total_mut,
                              color=sample,show.legend=FALSE)+
  xlab('Somatic_MMR_PR_Mut_count')+
  ggtitle('Somatic MMRPD mut counts vs Total mut counts')

if (WRITE==TRUE){
  ggsave(filename = paste(outPutPath,'figs/','somMMRPDvsTotalPlot_',Sys.Date(),'.png',sep=''),
         plot=somMMRPDvsTotalPlot,dpi=300)
  ggsave(filename = paste(outPutPath,'figs/','germMMRPDvsTotalPlot_',Sys.Date(),'.png',sep=''),
         plot=germMMRPDvsTotalPlot,dpi=300)
}











