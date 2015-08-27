####################################################################
## thr33some prj. expr-expr net construction                      ##
## Copyright 2014 Francesco Gadaleta                              ##
## University of Liege - Belgium                                  ##
####################################################################

# merge all results into one adj matrix
# convert adj matrix into graph
# path: 
#Return value: generate xml graph to be visualised in Cytoscape

source("makegraph.R")

mergeResults <- function(path ="./results/", ee=T, em=T, plot=F) {
  #TODO merge tmp files into one network
  if(ee) {
    ee_path = paste(path,"expr_expr", sep = "")
    filenames <- list.files(ee_path)
    numfiles <- length(filenames)  
    files = paste(ee_path, filenames, sep="/")
    expr_expr = matrix(0, nrow=500, ncol=500)
  
    for(i in 1:numfiles)
    {
      tmp <- read.csv(files[i], header=T)
      if(i==1) genenames <- colnames(tmp)     # gene names equal for all 
      #t1 = as.matrix(tmp)
      #mode(t1) = "numeric"
      # TODO normalise matrix
      expr_expr = expr_expr + as.matrix(tmp)
    }
    

    # visualise it and save graph
    dbg <- sprintf("Saving graph XML visualisation file")
    print(dbg)
    makegraph(expr_expr, nodesize=6, labels=genenames, 
              saveto = paste(path, "ee_graph.xml", sep = ""), plot = plot)
  
    
    # save merged file to disk
    dbg <- sprintf("Saving adj matrix or merged results (expr-expr)")
    print(dbg)
    write.matrix(format(as.data.frame(expr_expr), scientific=F), 
                 file = paste(path, "expr_expr_merged.cvs", sep=""))
    
  }

## TODO em_path
if(em) {
  em_path = paste(path,"expr_meth", sep = "")
  filenames <- list.files(em_path)
  numfiles <- length(filenames)  
  files = paste(em_path, filenames, sep="/")
  expr_meth = matrix(0, nrow=500, ncol=500)
  
  for(i in 1:numfiles)
  {
    #tmp <- read.csv(files[i], header=T)
    tmp = read.table(files[i], header=T, quote="\"")
    
    if(i==1) genenames <- colnames(tmp)     # gene names equal for all 
    #t1 = as.matrix(tmp)
    #mode(t1) = "numeric"
    # TODO normalise matrix
    expr_meth = expr_meth + as.matrix(tmp)
  }
  
  # visualise it and save graph
  # visualise it and save graph
  dbg <- sprintf("Saving expr-methyl graph XML visualisation file")
  print(dbg)
  makegraph(expr_meth, nodesize=6, labels=genenames, 
            saveto = paste(path, "em_graph.xml", sep = ""), plot = plot)
  
  dbg <- sprintf("Saving adj matrix or merged results (expr-methyl)")
  print(dbg)
  write.matrix(format(as.data.frame(expr_expr), scientific=F), 
               file = paste(path, "expr_methy_merged.cvs", sep=""))
  
}


# l = length(which(netw>0, arr.ind=T))
# sub = which(netw>0, arr.ind=T)
# rownames(sub) = NULL
# nzidx = which(netw>0, arr.ind=T)[1:l]
# 
# 
# g <- graph.adjacency(newt, weighted=T, mode="undirected", diag=F)
# V(g)$name = labels   # fix this when generalising
# 
# if(outfile == "gml")
#   write.graph(graph=g, file="hierarchical.xml", format="graphml")
# 
# 
}
