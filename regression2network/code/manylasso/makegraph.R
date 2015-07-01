####################################################################
## thr33some prj. expr-expr net construction                      ##
## Copyright 2014 Francesco Gadaleta                              ##
####################################################################

library(igraph)

makegraph <- function(adjmat, nodesize=10, labels="", saveto="", outfile="gml", 
                      layout = layout.fruchterman.reingold, plot=F)
{
  if(!is.matrix(adjmat)) 
    adjmat = as.matrix(adjmat)
  
  g <- graph.adjacency(adjmat, weighted=T, mode="undirected", diag=F)
  V(g)$name = labels
  
  if(outfile == "gml")
    write.graph(graph=g, file=saveto, format="graphml")
  
  if(plot)
    plot.igraph(g, vertex.label=V(g)$name, 
                edge.color="black",
                vertex.size=nodesize,
                vertex.label.dist=0,
                layout=layout.fruchterman.reingold)
}
