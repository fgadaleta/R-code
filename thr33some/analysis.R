# analysis code 
# this is part of the Thr33s0m3 project
# Copyright 2014 Francesco Gadaleta



library(sna)
library(network)
library(igraph)
library(statnet)
library(ggplot2)

# read adj matrix networks
expr_expr <- read.table("./results/expr_expr_12042_merged.cvs", 
                                     header=TRUE, quote="\"")

expr_meth <- read.table("./results/expr_methyl_12042_merged.cvs", 
                        header=TRUE, quote="\"")

genenames = colnames(expr_expr)


# AND matrix 
maxnum = as.integer(dim(expr_expr)[1])
a = expr_expr[1:maxnum, 1:maxnum]
b = expr_meth[1:maxnum, 1:maxnum]
c = `&`(a,b)
mode(c) = "integer"

#expr_and_meth =  `&`(expr_expr, expr_meth)


# build network (sna)
a = symmetrize(a, rule = "weak")
ee_net = graph.adjacency(a, mode = "undirected")

b = symmetrize(b, rule = "weak")
em_net = graph.adjacency(b, mode ="undirected")
and_net = graph.adjacency(c, mode = "undirected")

# compute vectors of degrees
ee_degree = igraph::degree(ee_net)
em_degree = igraph::degree(em_net)
and_degree = igraph::degree(and_net)

#purge non significant predictors
qee = quantile(ee_degree, probs = c(.25,.5,.75,.95))
qem = quantile(em_degree, probs = c(.25,.5,.75,.95))

genes = genenames[which(em_degree>qem[4], arr.ind = T)]
genes = genenames[which(ee_degree>qee[4], arr.ind = T)]
genes = genenames[which(and_degree>0, arr.ind = T)]

source("getGeneKEGGinfo.R")
res <- getGeneKEGGinfo(genes)
query <- res[[1]]
pathway <- res[[2]]
#just order it by significance
pathway = pathway[order(pathway$Percentage, decreasing = T),]

# prepare and plot data 
dat <- data.frame(
  datatype = rep(c("ee", "em"), each=maxnum),
  x = rep(1:maxnum, 2),
  y = c(ee_degree, -em_degree)
)

ggplot(dat, aes(x=x, y=y, fill=datatype)) + 
  geom_bar(stat="identity", position="identity")
