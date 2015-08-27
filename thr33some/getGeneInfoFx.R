#This function will get biological info for a given set of gene symbols
#Function uses multiple Bio databases
#Author: Kirill Bessonov


#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#biocLite("org.Hs.eg.db")
#biocLite("KEGGprofile")

library(biomaRt)
library(org.Hs.eg.db)
library(KEGGprofile)


ensemblCxn <- useMart("ensembl"); #connect to Ensemble DB
ensemblDB <- useDataset("hsapiens_gene_ensembl", mart=ensemblCxn); #create object using human ensembl dataset


EnrtrezID2KeggID <- function(x){
	if(length(x) == 1){stop("More than 1 EntrezID required");}
	entrezIDs <- mappedkeys(org.Hs.egPATH)
	entrezIDtoKeggList <-as.list(org.Hs.egPATH[entrezIDs])
	
	return(entrezIDtoKeggList[x])	
}



#annotate gene to pathways
genes = c("A2M","A4GALT","PGK2","PGK1", "PDHB", "PDHA1","PDHA2","PGM2","TPI1","ACSS1")
          

queryGene <- getBM(
  attributes=c("ensembl_gene_id", "wikigene_name", "entrezgene","wikigene_description","chromosome_name","gene_biotype"),
  filters="wikigene_name",values=genes, mart=ensemblDB)


#perform KEGG pathway enrichment analysis				
KEGGenrichAnal <- find_enriched_pathway(queryGene[,"entrezgene"], specis="hsa", returned_pvalue=0.05)[[1]]




