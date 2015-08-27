####################################################################
# data extraction utilities
####################################################################

assert <- function (expr, error) {
  if (! expr) stop(error, call. = FALSE)
}



# one-to-one mapping methyl id into gene space
# usage: out = extract_gene_from_methyl(namesmethyl)
extract_gene_from_methyl <- function(namesmethyl)
{
  end = length(namesmethyl)
  genespace = {}
  for(i in 1:end) {
    ## get methyl_id
    methyl_id = namesmethyl[i]
    methylgene = strsplit(methyl_id, "_")[[1]][1]
    genespace = c(genespace, methylgene)
  }
  return(genespace)
}