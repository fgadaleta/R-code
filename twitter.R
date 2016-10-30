library(twitteR)

fragadaletaTweets <- userTimeline("worldofpiggy", n=1000)
(nDocs <- length(fragadaletaTweets))
fragadaletaTweets[11:15]


#convert to data frame
df <- do.call("rbind", lapply(fragadaletaTweets, as.data.frame))
dim(df)


#transform text
library(tm)
myCorpus <- Corpus(VectorSource(df$text))
myCorpus <- tm_map(myCorpus,tolower)
#no punctuation
myCorpus <- tm_map(myCorpus, removePunctuation)
#no numbers
myCorpus <- tm_map(myCorpus, removeNumbers)
#no urls
removeURL <- function(x) gsub("http[[:alnum:]]*", "", x)
myCorpus <- tm_map(myCorpus,removeURL)
#add extra stop words
myStopwords <- c(stopwords('english'), "available", "via")
#remove r and big
idx <- which(myStopwords %in% c("r","big"))
myStopwords <- myStopwords[-idx]
#remove stopwords from corpus
myCorpus <- tm_map(myCorpus, removeWords, myStopwords)

#stemming words
#library(Snowball)
myCorpusCopy <- myCorpus
#stem
#myCorpus<- tm_map(myCorpus, stemDocument)
myCorpus <- tm_map(myCorpus, stemCompletion, dictionary=myCorpusCopy)
inspect(myCorpus[11:15])

#term-document matrix
myTdm <- TermDocumentMatrix(myCorpus, control=list(wordLengths=c(1,Inf)))
myTdm

#check the first six terms starting with 'r'
idx <- which(dimnames(myTdm)$Terms == "r")
#inspect(myTdm[idx+(0:5),101:110])

myTdm <- TermDocumentMatrix(myCorpus, control = list(minWordLength=1))
rownames(myTdm)   # list of terms

# frequent terms
findFreqTerms(myTdm, lowfreq=10)
termFrequency <- rowSums(as.matrix(myTdm))
termFrequency <- subset(termFrequency, termFrequency>=10)
library(ggplot2)
qplot(names(termFrequency), termFrequency, geom="bar") + coord_flip()

