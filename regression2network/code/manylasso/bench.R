################################################
##  Manylasso Survey                          ##
##  (c) 2014 Francesco Gadaleta               ##    
##                                            ##  
################################################

library(sna)

bench<- function(gold, predicted) {  
  
  #goldadj = symmetrize(sign(gold), rule="weak")  
  goldadj = gold                          # gold should be binary
  pred = symmetrize(predicted, rule="weak")
  pred = abs(sign(pred))
  
  tp = 0  # true positives 
  tn = 0  # true negatives
  fp = 0  # false positives
  fn = 0  # false negatives
  tpr = fpr = acc = 0
  
  gold_edges = sum(goldadj)   # edges from the directed net 
  pred_edges = sum(pred)
  
  for(i in 1:p){
    for(j in 1:p) {
      if(goldadj[i,j] == pred[i,j]) {
        if(pred[i,j] == 1 )
          tp = tp + 1
        if(pred[i,j] == 0)  
          tn = tn + 1
      }
      
      if(goldadj[i,j] != pred[i,j]){
        if(pred[i,j] == 1)
          fp = fp+1
        if(pred[i,j] == 0)
          fn=fn+1
      }
    }
  }
  
  # compute mcc (Matthew's Correlation Coefficient)
  mcc = (tp*tn - fp*fn)/ sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  tpr = tp/(tp+fn)
  fpr = fp/(fp+tn)
  acc = (tp+tn)/(tp+tn+fp+fn)
  output = c(gold_edges, pred_edges, tp,fp,tn,fn,mcc, tpr, fpr, acc)
  return(output)
}

benchNetwork <- function(gold, predicted, file) {  
  result = matrix(0,ncol=16, nrow=dim(gold)[1])
  colnames(x=result) = c("gold degree", "pred degree", "cor degree", "gold between", 
                         "pred between", "cor between", "gold num comp", "pred num comp", 
                         "gold connectedness", "pred connectedness", "gold num comp", "pred num comp", 
                         "gold largest comp", "pred largest comp", "true component", "false component")
  rownames(x=result) = colnames(gold)            
  
  #goldadj= symmetrize(gold, rule="weak")
  goldadj = gold
  
  if(sum(predicted > 0)) {
    pred = symmetrize(predicted, rule="weak")
    pred = abs(sign(pred))
    # network properties
    # create sna network objects
    goldnet = network(x=goldadj, directed=F) 
    prednet = network(x=pred, directed=F)     
    
    # local properties
    result[, "pred degree"] = sna::degree(prednet)
    result[, "gold degree"]  = sna::degree(goldnet)
    result[, "gold between"] = sna::betweenness(goldnet)
    result[, "pred between"] = sna::betweenness(prednet)
    result[1, "cor degree"] = cor(result[, "gold degree"], result[, "pred degree"])
    result[1, "cor between"] = cor(result[, "gold between"], result[, "pred between"])                     
    #global properties
    result[1, "gold connectedness"]  = sna::connectedness(goldnet)
    result[1, "pred connectedness"] = sna::connectedness(prednet)
    # number of components
    result[1, "gold num comp"] = sna::components(goldnet)  
    result[1, "pred num comp"] = sna::components(prednet)
    # nodes of the largest component 
    result[, "gold largest comp"] = sna::component.largest(goldnet)
    result[, "pred largest comp"] = sna::component.largest(prednet)
    result[1, "true component"] = sum(result[, "gold largest comp"] ==  result[, "pred largest comp"])                    
    result[1, "false component"] = sum(result[, "gold largest comp"] !=  result[, "pred largest comp"])
  }    
  
  else{
    result[,] = 0 
  }
  
  write.csv(result, file)
  return (result)  
}


