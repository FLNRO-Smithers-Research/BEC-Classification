#Coding by K.Przeczek, M. Sakals, J. Edinger, and K. Pitman, Feb. 7, 2013

#Kara Pitman's soil moisture - species data
#Output frequency of occurrence of a species within a soil moisture zone, weighted occurrence for each 
#zone, and percentage of observations of each species within a zone.

#Read in data
sm <- read.csv("C:/Users/kpitman/R/All.csv")
result=read.csv("C:/Users/kpitman/R/result.csv") 

#Check data
str(sm)
#set WNA as a factor
sm$WNA <- as.factor(sm$WNA)

#prepare a table with weighting factors
weight <- read.csv("C:/Users/kpitman/R/weighting.csv") 
str(weight)

#The guts:
p=1
for(i in 1:nlevels(sm$Species)){ #run the loop for every factor level in the column "species"
  tmp.spp <- subset(sm, Species == levels(sm$Species)[i]) #subset data set by species level i
    if (nrow(tmp.spp) < 30) {} #if there are fewer than 30 observations of a given species, do not include it.
    else {
    t1 <- table(factor(tmp.spp$WNA,levels=seq(1:9))) #gives freqency of each WNA. COrresponds to histogram results.
    t1 <- as.data.frame(t1) #convert to a data frame
    t1$wt <- weight$wt #add a column with the weight factors in line with appropriate WNA
    t2 <- t1[order(-t1$Freq, t1$Var1, t1$wt),] #sort crosstabulation table by Frequency, from largest to smallest
    rank <- c(1:9) #create column for rank
    species <- rep(levels(sm$Species)[i], 9)
    t3 <- cbind(species, rank, t2) #add that column at the beginning of the table. Otherwise it is added at the end.
    colnames(t3) <- c("Species", "Rank","WNA","Frequency", "Weighting") #rename the columns
    t3$Percent <- round((t3$Frequency*t3$Weighting)/sum(t3$Frequency*t3$Weighting)*100, digits=2) #Calculates the percent of weighting 
    #t3$Percent <- round(t3$Frequency/sum(t3$Frequency)*100, digits=2) #calculate the percent of the total observations within each WNA, doens't account for weighting of plot numbers.
    t3 <- t3[order(-t3$Percent, t3$Weighting, t3$Frequency),] #Puts Percent in the top column in order to extract highest value 
    
    result[i,2]=as.character(t3$Species[2])
    result[i,3]=round(t3$Percent[1],2)
    result[i,4]=t3$WNA[1]
    result[i,5]=round(t3$Percent[1]+t3$Percent[2],2) 
    result[i,6]=t3$WNA[2]
    result[i,7]=round(t3$Percent[1]+t3$Percent[2]+t3$Percent[3],2)
    result[i,8]=t3$WNA[3]
    result[i,9]=round(t3$Percent[1]+t3$Percent[2]+t3$Percent[3]+t3$Percent[4],2)
    result[i,10]=t3$WNA[4]
    result[i,11]=round(t3$Percent[1]+t3$Percent[2]+t3$Percent[3]+t3$Percent[4]+t3$Percent[5],2)
    result[i,12]=t3$WNA[5]
    
    }
     
   write.table(result, "C:/Users/kpitman/R/result.csv", sep = ",", row.names=FALSE, col.names=TRUE)
 }