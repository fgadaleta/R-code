#library(ggplot2)
library(reshape2)
#library(compiler)

#library(Rcpp)
#library(inline)

# rotating functions
rot <- function(x) (1:x %% x) +1
rotvec <- function(vec) vec[rot(length(vec))] 


# parking hell
set.seed(1)
minutes <- 60   # observation time in minutes
duration <- 10  # parking duration in minutes
numlots <- totalcars*2   # maximum number of parking lots

# number of cars per time unit and total
cars <- rpois(minutes, lambda=8)
cars[cars==0] <- 1
totalcars <- sum(cars)

prefs <- matrix(0, nr=0, nc=numlots)
for(cc in 1:totalcars) {
  carPrefs <- sample(1:numlots, replace=F, size=numlots)
  prefs <- rbind(prefs, carPrefs)
}

#logical variables 
served <- numeric(totalcars)    # cars that have been served
busy <- numeric(numlots)         # parking lot availability and duration
j <- 0

price <- matrix(2, nr=numlots, nc=minutes)      # price of parking lots during observation time
attempt <- matrix(0, nr=totalcars, nc=numlots)  # number of attempts of driver for parking lot
tts <- numeric(totalcars)                       # time to be served

for (t in 1:minutes) {              # for each minute (observed time)
  tts[which(served==0)] <- tts[which(served==0)]+1          # clock is ticking for those who are not served
  
  for(car in 1:cars[t] ) {          # for each car in that minute
    j<-j+1
    
  if(served[j] == 0) {                    # current driver has not been served yet
    space <- prefs[j, 1]                  # decide where to park 
    # TODO chose according to price too 
    
    if(busy[space] == 0) {              # if available
      busy[space] = duration            # set duration
      served[j] <- 1                    # car has been served
      break
    }
    else {                                 # not available
      attempt[j, space] <- attempt[j, space] + 1
      # update preferences and break
      rotvec(prefs[j,])
      break
    }
  
  }
    
    # TODO update price of parking lots  
    for(cc in 1:numlots) {
        price[cc, t] <- price[cc, t] + 0.1*sum(attempt[, numlots]) 
        
    }
    
    # update duration 
    busy[busy>0] <- busy[busy>0]-1
    
}     


################################## visualisation ##################################

palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)
plot(palette)
#attempt <- as.matrix(scale(attempt))
heatmap(attempt, scale='none', col=palette)

which.max(preference[1,])   # select parking lot with maximum preference 

hc.rows <- hclust(dist(attempt))
hc.cols <- hclust(dist(t(attempt)))
plot(hc.rows)
heatmap(attempt[cutree(hc.rows, k=2)==1,],Colv=as.dendrogram(hc.cols), scale='none')
