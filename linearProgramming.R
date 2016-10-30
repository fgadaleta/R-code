library(lpSolveAPI)

#used for result visualization
library(ggplot2)
library(reshape)
library(gridExtra)

#define the datasets

train<-data.frame(wagon=c('w1','w2','w3'), weightcapacity=c(10,8,12), spacecapacity=c(5000,4000,8000))
cargo<-data.frame(type=c('c1','c2','c3','c4'), available=c(18,10,5,20), volume=c(400,300,200,500),profit=c(2000,2500,5000,3500))


# /* Objective function */
#   max: +2000 C1 +2500 C2 +5000 C3 +3500 C4 +2000 C5 +2500 C6 +5000 C7 +3500 C8 +2000 C9 +2500 C10 +5000 C11
# +3500 C12;
# 
# /* Constraints */
#   +C1 +C2 +C3 +C4 <= 10;
# +C5 +C6 +C7 +C8 <= 8;
# +C9 +C10 +C11 +C12 <= 12;
# +400 C1 +300 C2 +200 C3 +500 C4 <= 5000;
# +400 C5 +300 C6 +200 C7 +500 C8 <= 4000;
# +400 C9 +300 C10 +200 C11 +500 C12 <= 8000;
# +C1 +C5 +C9 <= 18;
# +C2 +C6 +C10 <= 10;
# +C3 +C7 +C11 <= 5;
# +C4 +C8 +C12 <= 20;



#create model with 10 constraints and 12 decision vars
lpmodel <- make.lp(2*nrow(train)+2*nrow(cargo), 12)

#I used this to keep count within the loops, I admit that this could be done a lot neater
column<-0
row<-0

#build the model column per column
for(wg in train$wagon){
  row<-row+1
  for(type in seq(1,NROW(cargo$type))){
    column<-column+1
    
    #this takes the arguments 'column','values' & 'indices' (as in where these values should be placed in the column)
    set.column(lpmodel,column,c(1, cargo[type,'volume'],1), indices=c(row,NROW(train)+row, NROW(train)*2+type))
  }}

#set rhs weight constraints
set.constr.value(lpmodel, rhs=train$weightcapacity, constraints=seq(1,NROW(train)))

#set rhs volume constraints
set.constr.value(lpmodel, rhs=train$spacecapacity, constraints=seq(NROW(train)+1,NROW(train)*2))


#set rhs volume constraints
set.constr.value(lpmodel, rhs=cargo$available, constraints=seq(NROW(train)*2+1,NROW(train)*2+NROW(cargo)))

#set objective coefficients
set.objfn(lpmodel, rep(cargo$profit,NROW(train)))

#set objective direction
lp.control(lpmodel,sense='max')

#I in order to be able to visually check the model, I find it useful to write the model to a text file
write.lp(lpmodel,'model.lp',type='lp')


# solve the model (if returns 0 an optimal solution found)
solve(lpmodel)

#this return the proposed solution
get.objective(lpmodel)

ggplot(lpmodel)


#the code below generates the plots
results<-data.frame(cargo=rep(cargo$type, 3), wagon=as.vector(sapply(train$wagon, FUN=function(x) rep(x, NROW(cargo)))), solution=get.variables(lpmodel))

r1<-ggplot(results, aes(x=cargo, y=solution, fill=wagon)) + geom_bar(color='black', position='dodge') + geom_text(aes(label=solution), size=2.5, position=position_dodge(width=1), vjust=-.4) + scale_fill_brewer(palette='Set1') + facet_grid(.~wagon) + opts(title='Planning result', legend.position='none') + ylab('Solution (tonnes)')

financialresult<-data.frame(cargo=rep(cargo$type, 3), wagon=as.vector(sapply(train$wagon, FUN=function(x) rep(x, NROW(cargo)))), solution=get.variables(lpmodel), profit_unit=rep(cargo$profit, 3))
financialresult$profit<-financialresult$profit_unit*financialresult$solution

r2<-ggplot(financialresult, aes(x=cargo, y=profit, fill=wagon)) + geom_bar(color='black', position='dodge') + geom_text(aes(label=profit), size=2.5, position=position_dodge(width=1), vjust=-.4) + scale_fill_brewer(palette='Set1') + facet_grid(.~wagon) + opts(title='Financial result', legend.position='none') + ylab('Solution ($)')

spacecapacity<-data.frame(wagon=train$wagon, capacity=train$spacecapacity, used=get.constraints(lpmodel)[4:6])
weightcapacity<-data.frame(wagon=train$wagon, capacity=train$weightcapacity, used=get.constraints(lpmodel)[1:3])

c1<-ggplot(melt(spacecapacity, id='wagon'), aes(x=variable, y=value, fill=variable)) + geom_bar(color='black') + facet_grid(.~wagon) + opts(legend.position='none', title='Volume capacity vs. used') + scale_fill_brewer(palette='Set1') + geom_text(aes(label=value), size=2.5, vjust=-.4) + ylab('Volume (m2)')

c2<-ggplot(melt(weightcapacity, id='wagon'), aes(x=variable, y=value, fill=variable)) + geom_bar(color='black') + facet_grid(.~wagon) + opts(legend.position='none', title='Weight capacity vs. used') + scale_fill_brewer(palette='Set1') + geom_text(aes(label=value), size=2.5, vjust=-.4) + ylab('Weight (tonnes)')

png('model_results.png', width=8, height=8, units='in', res=90)
grid.arrange(r1, r2, c1, c2, ncol=2)
dev.off()


