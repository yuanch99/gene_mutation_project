library(tidyverse)
library(biomaRt)
library(ggplot2)
WRITE <- FALSE

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

#Annotations
somAnnoList <- unique(somMutData$Annotation)
germAnnoList <- unique(germMutData$Annotation)
sepAnnos <- function(annoList){
  list <- c()
  for (i in annoList){
    for (j in strsplit(i,"&")){
      for (k in j){list <- append(list,j)}
    }
  }
  return(list)
}

annoList <- sort(unique(append(sepAnnos(somAnnoList),sepAnnos(germAnnoList))))
rm(somAnnoList,germAnnoList)


#select only coding and exon region mutation
importantAnnos <- paste(
  "5_prime_UTR_premature_start_codon_gain_variant","disruptive_inframe_deletion",
  "disruptive_inframe_insertion","disruptive_inframe_insertion","exon_loss_variant",
  "frameshift_variant","inframe_deletion","inframe_insertion","initiator_codon_variant",
  "missense_variant","splice_acceptor_variant","splice_donor_variant",
  "splice_region_variant","start_lost","stop_gained","stop_lost","stop_retained_variant",
  "synonymous_variant",
  sep = "|")

somImpMut <- somMutData %>% filter(grepl(importantAnnos,Annotation))
germImpMut <- germMutData %>% filter(grepl(importantAnnos,Annotation))

#gene list info

somGene <- c(unique(somMutData$Gene_Name))
germGene <- c(unique(germMutData$Gene_Name))
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
somGeneInfo <- getBM(
  attributes=c("hgnc_symbol",'transcript_length','start_position','end_position',
               'ensembl_gene_id','chromosome_name'),
  filters = "hgnc_symbol", values = somGene, mart = ensembl)
germGeneInfo <- getBM(
  attributes=c("hgnc_symbol",'transcript_length','start_position','end_position',
               'ensembl_gene_id','chromosome_name'),
  filters = "hgnc_symbol", values = germGene, mart = ensembl)

#remove duplicates by selecting the one with largest transcript_length
somGeneInfo <- arrange(somGeneInfo,hgnc_symbol,desc(transcript_length))
germGeneInfo <- arrange(germGeneInfo,hgnc_symbol,desc(transcript_length))
somGeneInfo <- somGeneInfo[!duplicated(somGeneInfo$hgnc_symbol),]
germGeneInfo <- germGeneInfo[!duplicated(germGeneInfo$hgnc_symbol),]

##should find a way to speed up this step
count_gene_muts <- function(geneInfo,impMut){
  mutCount <- c()
  for (i in geneInfo$hgnc_symbol){
    mutCount <- append(mutCount,length(impMut[impMut$Gene_Name==i,1]))
  }
  return(mutCount)
}

somGeneInfo$mutation_count <- count_gene_muts(somGeneInfo,somImpMut)
germGeneInfo$mutation_count <- count_gene_muts(germGeneInfo,germImpMut)

somGeneInfo$mut_count_by_trans_len <- somGeneInfo$mutation_count/somGeneInfo$transcript_length
germGeneInfo$mut_count_by_trans_len <- germGeneInfo$mutation_count/germGeneInfo$transcript_length


somGeneInfo <- arrange(somGeneInfo,desc(mut_count_by_trans_len))
germGeneInfo <- arrange(germGeneInfo,desc(mut_count_by_trans_len))

if (WRITE==TRUE){
  write_csv(somGeneInfo,
            paste(outPutPath,'data/','somGeneInfo_',Sys.Date(),'.csv',sep = ''))
  write_csv(germGeneInfo,
            paste(outPutPath,'data/','germGeneInfo_',Sys.Date(),'.csv',sep = ''))
  write.table(strsplit(importantAnnos,split = '\\|'),
            paste(outPutPath,'data/','importantAnnos_',Sys.Date(),'.txt',sep = ''),
            sep = ',',col.names = F, row.names = F)
  }

#some basic plots
#mutation count by chromosomes
count_chr_muts <- function(mutData){
  mutCount <- c()
  chr <- c()
  for (i in unique(mutData$Chromosome)){
    chr <- append(chr,i)
    mutCount <- append(mutCount,length(mutData[mutData$Chromosome==i,1]))
  }
  chrMutCount <- data.frame(chromosome=chr,count=mutCount)
  chrMutCount <- chrMutCount[!grepl('*_',chrMutCount$chromosome),]
  chrMutCount$length <- as.numeric(c("248956422","242193529","198295559","190214555","181538259",
                          "170805979","159345973","145138636","138394717","133797422",
                          "135086622","133275309","114364328","107043718","101991189",
                          "90338345","83257441","80373285","58617616","64444167",
                          "46709983","50818468","156040895","57227415"))
  chrMutCount$mut_by_length <- chrMutCount$count/chrMutCount$length
  chrMutCount <- arrange(chrMutCount,desc(mut_by_length))
  return(chrMutCount)
}

somChrMutCount <- count_chr_muts(somMutData)
germChrMutCount <- count_chr_muts(germMutData)

somChrMutPlot <- ggplot(somChrMutCount,aes(x=reorder(chromosome,-mut_by_length)))+
  geom_bar(aes(y=count),stat='identity',fill="#69b3a2")+
  geom_point((aes(y=length/13146)))+
  scale_y_continuous(
    name = "Mutation Count",
    sec.axis = sec_axis(~.*13146,name="Chromosome Length")
  )+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Chromosome Mutation Count (Somatic)")
  
germChrMutPlot <- ggplot(germChrMutCount,aes(x=reorder(chromosome,-mut_by_length)))+
  geom_bar(aes(y=count),stat='identity',fill="#69b3a2")+
  geom_point((aes(y=length/13146)))+
  scale_y_continuous(
    name = "Mutation Count",
    sec.axis = sec_axis(~.*13146,name="Chromosome Length")
  )+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Chromosome Mutation Count (Germline)")

if (WRITE==TRUE){
  ggsave(filename = paste(outPutPath,'figs/','somChrMutPlot_',Sys.Date(),'.png',sep=''),
                          plot=somChrMutPlot,dpi=300)
  ggsave(filename = paste(outPutPath,'figs/','germChrMutPlot_',Sys.Date(),'.png',sep=''),
         plot=germChrMutPlot,dpi=300)
}


germChrMutbyLenPlot <- qplot(data=germChrMutCount,
                         x=reorder(chromosome,-mut_by_length),y=mut_by_length)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Mutation by Length (Germline)")

somChrMutbyLenPlot <- qplot(data=somChrMutCount,
                        x=reorder(chromosome,-mut_by_length),y=mut_by_length)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Mutation by Length (Somatic)")

if (WRITE==TRUE){
  ggsave(filename = paste(outPutPath,'figs/','somChrMutbyLenPlot_',Sys.Date(),'.png',sep=''),
         plot=somChrMutbyLenPlot,dpi=300)
  ggsave(filename = paste(outPutPath,'figs/','germChrMutbyLenPlot_',Sys.Date(),'.png',sep=''),
         plot=germChrMutbyLenPlot,dpi=300)
}
