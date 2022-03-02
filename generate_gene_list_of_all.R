library(tidyverse)

#Load Data
dataPath <- '~/tabori/data/mutation/raw_data/santiago_pipeline/vcf_annotation/'

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

#Annotations
somAnnoList <- unique(somMutData$Annotation)
germAnnoList <- unique(germMutData$Annotation)
sepAnnos <- function(annolist){
  list <- c()
  for (i in annoList){
    for (j in strsplit(i,"&")){
      print(j)
      list <- append(list,j)
    }
  }
  return(list)
}
somAnnoList <- sepAnnos(somAnnoList)
germAnnoList <- sepAnnos(germAnnoList)
unique(somAnnoList,germAnnoList)


#select only coding and exon region mutation
interestedAnnos <- c(
  "synonymous_variant","missense_variant","splice_region_variant",
  "stop_retained_variant","non_coding_exon_variant",
  "5_prime_UTR_premature_start_codon_gain_variant","intergenic_region",
  "intragenic_variant","splice_donor_variant","splice_acceptor_variant",
  "splice_region_variant","frameshift_variant","initiator_codon_variant",
  "start_lost","stop_gained","stop_lost","inframe_insertion","disruptive_inframe_deletion",
  "inframe_deletion","disruptive_inframe_insertion")

somCodeMut <- somMutData %>% 
  filter(grepl(paste(interestedAnnos,sep = "|"),strsplit(Annotation,"&"),fixed = TRUE))

