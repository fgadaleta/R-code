source("load_data.R")
source("build_expr_expr_net.R")
source("build_expr_methyl_net.R")
source("merge_results.R")

exprnet = build_expr_expr_net(exprdata, start=1, end=100,ncpu=4)
methnet = build_expr_methyl_net(exprdata=exprdata, methyldata=methyldata, start=1, end=500)

#tmp = read.csv("./results/expr_expr_subnet_500_start_3_end_115.csv",header=T)

# merge results and create graph file 
mergeResults(ee = T, em=T, plot=F)

makegraph(exprnet, outfile = "nofile", plot = T)
makegraph(methnet, outfile = "nofile", plot = T)

source("makegraph.R")
