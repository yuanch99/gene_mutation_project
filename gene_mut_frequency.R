library(tidyverse)
library(biomaRt)
library(ggplot2)
WRITE <- TRUE

dataPath <- '/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/raw_data/santiago_pipeline/vcf_annotated/'
outPutPath <- "/Users/yuanchang/Documents/MBP/TaboriLab/codes/mutation_project/gene_mutation_project/output/"

poleSom <- read.csv(paste(dataPath,'all_snpEff_annotations_somatic_pole_22-03-20.csv',sep = ''))
pold1Som <- read.csv(paste(dataPath,'all_snpEff_annotations_somatic_pold1_22-03-20.csv',sep = ''))
# poleGerm <- read.csv(paste(dataPath,'all_snpEff_annotations_germline_pole_22-03-20.csv',sep = ''))
# pold1Germ <- read.csv(paste(dataPath,'all_snpEff_annotations_germline_pold1_22-03-20.csv',sep = ''))
som <- rbind(poleSom,pold1Som)
rm(poleSom,pold1Som)

#select missense variants
importantAnnos <- paste("missense_variant",sep = "|")
somImpMut <- som %>% filter(grepl(importantAnnos,Annotation))

mutGeneCount <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(mutGeneCount) <- c('gene','length','count')
mutGeneCount[, c(2:3)] <- sapply(mutGeneCount[, c(2:3)], as.numeric)
#mutGeneCount[, c(3)] <- sapply(mutGeneCount[, c(3)], as.numeric)
mutGeneCount$gene <- as.character(mutGeneCount$gene)
mutGenes <- c(unique(somImpMut$Gene_ID))

#get length of gene
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
geneLength <- getBM(
  attributes=c("hgnc_symbol",'transcript_length'),
  filters = "hgnc_symbol", values = mutGenes, mart = ensembl)
geneLength <- arrange(geneLength,hgnc_symbol,desc(transcript_length))
geneLength <- geneLength[!duplicated(geneLength$hgnc_symbol),]

for (i in 1:dim(geneLength)[1]){
  genename <- geneLength[i,1]
  genecount <- length(unique(somImpMut[somImpMut$Gene_ID==genename,]$Sample))
  genelength <- geneLength[i,2]
  mutGeneCount <- rbind(mutGeneCount,
                        data.frame('gene' = genename, 'length' = genelength, 'count' = genecount))
}

mutGeneCount$rate <- mutGeneCount$length/mutGeneCount$count
mutGeneCount <- mutGeneCount[order(mutGeneCount$rate),]

#select genes in mini panel
mini_panel <- read_csv('/Users/yuanchang/Documents/MBP/TaboriLab/data/mutation/raw_data/mini_panel.csv')
mutGeneCount_mini_panel <- mutGeneCount[mutGeneCount$gene%in%mini_panel$Associated_Gene_Name,]
`%!in%` <- Negate(`%in%`)
not_in_mini_panel <- mini_panel[mini_panel$Associated_Gene_Name%!in%mutGeneCount$gene,]

#add percentage of each gene
mutGeneCount_mini_panel$percentage <- 100*percent_rank(mutGeneCount_mini_panel$rate)
mutGeneCount$percentage <- 100*percent_rank(mutGeneCount$rate)

if (WRITE==TRUE){
  write_csv(mutGeneCount,
            paste(outPutPath,'data/','mutGeneCount_somatic_',Sys.Date(),'.csv',sep = ''))
  write_csv(mutGeneCount_mini_panel,
            paste(outPutPath,'data/','mutGeneCount_mini_panel_somatic_',Sys.Date(),'.csv',sep = ''))
}

mutGeneCount2 <- mutGeneCount[mutGeneCount$count>1,]
