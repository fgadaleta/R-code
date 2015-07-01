################################################
##  Manylasso Survey                          ##
##  (c) 2014 Francesco Gadaleta               ##    
##                                            ##  
################################################

# import golden standard in matrix form

importGold <- function(filename, labels) {
  goldstd <- read.table(filename, header=F, sep="\t", fileEncoding="windows-1252")
  goldadj <- matrix(0, ncol=p, nrow=p)
  colnames(goldadj) <- labels
  rownames(goldadj) <- labels
  
  for(i in 1:nrow(goldstd)) {
    if(goldstd[i,3] == 1)
    {
      from = as.character(goldstd[i,1])
      to = as.character(goldstd[i,2])
      goldadj[from, to] = 1
    }
  }
  return(goldadj)
}
