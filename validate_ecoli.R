source("regulon_tf_adj.R")
true_net = regulon_tf_adj(file = "delete.txt")
pred_net <- as.matrix(read.table("~/Desktop/pred_net_ecoli_checklist", 
                                 header=TRUE, quote="\""))
pred_genes = colnames(pred_net) = tolower(colnames(pred_net))
rownames(pred_net) = pred_genes
true_genes = colnames(true_net)
idx = na.omit(match(true_genes, pred_genes))   # index of pred_genes in true_genes

# symmetrise 
library(sna)
tmp = symmetrize(pred_net)
rownames(tmp) = colnames(tmp) = pred_genes

# TODO create matrix of pred with only existent genes (submatrix of pred_net)
sub_pred = matrix(0, nrow = 1303, ncol=1303)
rownames(sub_pred) = colnames(sub_pred) = pred_genes[idx]

for(i in 1:1303) {
  from = pred_genes[idx[i]]
  to = names(which(tmp[, from] >0, arr.ind = T))
  if(length(to)) 
    if(to %in% pred_genes[idx]) {
      to = to[to %in% pred_genes[idx]]
      sub_pred[from, to] = sub_pred[to, from] = 1  
    }
}



source("makegraph.R")
makegraph(adjmat = pred_net, nodesize = 3, plot = T, saveto = "deleteme")
makegraph(adjmat = sub_pred, nodesize = 3, plot = T, saveto = "deleteme")

sub_genes= colnames(sub_pred)
which(sub_pred[sub_genes[34], ]>0)

tp = 0

check = which(sub_pred>0, arr.ind = T)
for(t in 1:688) {
  i = check[t,1]
  j = check[t,2]
  g = rownames(check)[t]

  connected = c(which(true_net[g,]>0), which(true_net[,g]>0))
  if(g %in% connected)
    tp = tp+1
}
## 

small_true = true_net[sub_genes, sub_genes]

p = 1303
res = bench(gold = small_true, predicted = sub_pred)
names(res) = c("gold_edges", "pred_edges", "tp","fp","tn","fn","mcc", "tpr", "fpr", "acc")
#benchNetwork(gold = small_true, predicted = sub_pred, file = "compareEcoli.csv")

components(small_true)
components(sub_pred)
cl1 = component.largest(small_true)
cl2 = component.largest(sub_pred)
sub_genes[cl1]
sub_genes[cl2]

dp = sna::degree(sub_pred)
dg = sna::degree(small_true)
cor(dp,dg) 

bp = sna::betweenness(sub_pred)
bg = sna::betweenness(small_true)
cor(bp, bg)
which(bp>0, useNames = T, arr.ind = T)
which(bg >0)




connected_to <- function(genename, network) 
{
  genelist = colnames(network)
  if(as.character(genename) %in% genelist == FALSE) {
    dbg <- sprintf("Gene %s is not listed in this network", genename)
    print(dbg)
    return(0)
  }
  
else {
    genes1 = which(network[, as.character(genename) ] > 0, arr.ind = T)
    genes2 = which(network[as.character(genename), ] > 0, arr.ind = T)
    genes = unique(c(genes1, genes2))
    return(genelist[genes])
  }
  
  return(NULL)
}

genename = pred_genes[idx[2]]

connected_to(genename, true_net)
connected_to(genename, pred_net)
