library(biomaRt)

#get gene start and end locations
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

chr1_genes <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id',
                                 'hgnc_symbol','chromosome_name','start_position',
                                 'end_position'), filters ='chromosome_name',
                    values ="1", mart = ensembl)