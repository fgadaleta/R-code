

#load data
event_attendees <- read.csv("event_attendees.csv", header=TRUE)
events <- read.csv("events.csv", header=TRUE)
user_friends <- read.csv("user_friends.csv", header=TRUE)
users <- read.csv("users.csv", header=TRUE)
train <- read.csv("train.csv", header=TRUE)
test <- read.csv("test.csv", header=TRUE)


#find unique user ids
colnames(users)
allusers <- unique(users[,1])

#find unique event ids
colnames(events)
allevents <- unique(events[,1])

#find unique user ids in user_friends
colnames(user_friends)
allfriends <- unique(user_friends[,1])

colnames(event_attendees)



########################################################
# pick up a user
#tmpuser <- users[1]  
tmpuser <- 3044012

# load events of user 
userevents <- train[train$user == tmpuser, ]
numevents <- length(userevents)
i <- 1
tmpev <- userevents[i,2]

#search for info about event in the event table
eventInfo <- events[events$event_id == tmpev, ]


# check time of event and joining time
userInfo <- users[users$user_id == tmpuser,]
#eventDateTime <- unlist(strsplit(as.character(eventInfo$start_time), "T"))
#userDateTime <- unlist(strsplit(as.character(userInfo$joinedAt), "T"))
#as.Date(eventDateTime[1]) - as.Date(userDateTime[1])

as.Date(eventInfo$start_time) - as.Date(userInfo$joinedAt)   # + means that user can still partecipate



########################################################
connections <- user_friends[user_friends$user == tmpuser,]

friendList <- unlist(strsplit(as.character(connections$friends), " "))
length(friendList)  # number of friends




# assign a score(E) for all events of user U 




#methods: scoring metric, elastic net, lasso, bayesian network



