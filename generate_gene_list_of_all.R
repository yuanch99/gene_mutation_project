library(biomaRt)

#Load Data
dataPath <- '~/tabori/data/mutation/raw_data/santiago_pipeline/vcf_annotation/'

initial <- read.csv(paste(dataPath,'all_18_Somatic_snpEff_annotations_20220301.csv',sep=''),nrows=100)
classes <- sapply(initial,class)
rm(initial)

somaticMutData <- read.csv(
  paste(dataPath,'all_18_Somatic_snpEff_annotations_20220301.csv',sep=''))
germlineMutData <- read.csv(
  paste(dataPath,'all_18_Germline_snpEff_annotations_20220301.csv',sep=''))

rm(classes)
gc()

