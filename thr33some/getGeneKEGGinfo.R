# Date: 9 Sept 2014
# This file is part of the Thr33s0m3 project 
# Author: Kyrylo Bessonov (c) 2014

# function: getGeneKEGGinfo() - annotate gene symbols and run KEGG enrichmnet analysis
# AIM:this function annotates genes and runs KEGG pathway enrichment analysis
# input:  character vector of gene symbols
# output: annotated dataframes on genes and KEGG pathways enrichment

# ----------------------------------------------------------------------------------------
#example using genes from sugar processing pathway - glycolysis 
#genes = c("A2M","A4GALT","PGK2","PGK1", "PDHB", "PDHA1","PDHA2","PGM2","TPI1","ACSS1");
#print(genes)
#ans <- getGeneKEGGinfo(genes)
#print(ans)
# ----------------------------------------------------------------------------------------



getGeneKEGGinfo <- function(x){
	library(biomaRt); 
  library(KEGGprofile)

	ensemblCxn <- useMart("ensembl"); #connect to Ensemble DB
	ensemblDB <- useDataset("hsapiens_gene_ensembl", mart=ensemblCxn); #create object using human ensembl dataset
	
	#annotate genes
	queryGene <- getBM(attributes=c("ensembl_gene_id","entrezgene","wikigene_description", 
					"chromosome_name", "gene_biotype"), 
					filters="wikigene_name", value=x, mart=ensemblDB)
					
	#remove duplicates
	dupl_idx <- duplicated(queryGene[,2])
	queryGene <- queryGene[!dupl_idx,]
	
	#perform KEGG pathway enrichment analysis				
	KEGGenrichAnal <- find_enriched_pathway(queryGene[,"entrezgene"], specis="hsa", returned_pvalue=0.05)[[1]]

	return(list(queryGene, KEGGenrichAnal))
}


