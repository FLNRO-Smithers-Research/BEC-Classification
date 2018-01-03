# Codes Runs a CART analysis on a Vegetation matrix predicting an SU.
# Data sets are exported from Vpro (3column R for veg and a SU saved as an XLSX file) ll valid alternate variables and thresholds
# Ensure that you have working Datasets  in your working directory:
library("reshape2",lib.loc="/path/to/R-packages/")

#############################################
############# PREPARATION ###################
#############################################
install.packages("caret")
install.packages("tcltk")
install.packages("ddply")
install.packages("rattle")
install.packages("rpart")
install.packages("randomForest")
install.packages("rpart.plot")
install.packages("tcltk")
install.packages("plyr")
install.packages("reshape")
install.packages("reshape2")
install.packages("VSURF")
install.packages("ipred")
install.packages("base64")
install.packages("e1071")
install.packages("parallel")
install.packages("foreach")
install.packages("doParallel")
install.packages("Matrix")
install.packages ("labdsv")
install.packages ("vegan")
install.packages("gdata")
install.packages("MASS")
install.packages ("xlsx")
install.packages("C50")
# Required R packages:
library(caret)
library(tcltk)
library(ddply)
library(rattle)
library(rpart)
library(randomForest)
library(rpart.plot)
library(tcltk)
library(plyr)
library(reshape)
library(reshape2)
library(VSURF)
library(ipred)
library(base64)
library(e1071)
library(parallel)
library(foreach)
library(doParallel)
require(doBy)
library(rbenchmark)
library(Matrix)
require(labdsv)
require(bigmemory)
require(bigpca)
require(vegan)
require(gdata)
require(MASS)
require(xlsx)
require (C50)
# clean up workspace
rm(list=ls())
wd=tk_choose.dir(); setwd(wd)

# set file name - what do you want to call the Vegetation Summary file?
fname="VegetationMatrix"

# read in 3-column vegetation data
# user chooses file interactively
# choose the .txt file produced from VPro
a=file.choose()
veg=read.table(a,as.is=TRUE,header=TRUE,stringsAsFactors=FALSE)
setwd(gsub(basename(a),"",a))
colnames(veg)[1:3]=c("PlotNumber", "Species", "Cover")

#choose the CSV matrix of Site Units
b=file.choose()
su <- read.xlsx(b, sheetName = "Sheet1",as.is=TRUE, header=TRUE, row.names=1,  stringsAsFactors=FALSE)
su.order <- as.matrix(su[order(row.names(su)),])

####Count number of columns (species) required
Nspp = length(unique(veg$Species))
Listspp <- unique(veg$Species)
Plots <- unique(veg$PlotNumber)
####Counts the number of instances of each unique species
Countspp<-ddply(veg,~Species,summarise,sppcount=length(unique(PlotNumber)))

veg1 <- veg
##########Create Mtrix from List from package labdsv######
vegmtrx <- matrify(veg1)
vegmtrx[vegmtrx==0] <- 0.0001
sort (su.order, decreasing=F)
###remove uncommon species
vegmtrx2 <- dropspc(vegmtrx, minocc=2, minabu=.1)
###Merge SU and Veg matrices
allmtrx <- cbind(su.order,vegmtrx)

# save dataset as .Rda
save(allmtrx,file=paste(fname,".Rda",sep=""))


#############################################
############# START ANALYSIS ################
#############################################

X1 <- vegmtrx2
X1$SU <- as.factor (allmtrx$SiteUnit)

################run a CART on matrix##########
###Change cp value to show more or less detail
vegtree <- rpart(SU ~ ., data = X1, method ='class', control = rpart.control (minsplit = 5,  
              cp = 0.005, maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10))
##rpart.plot(vegtree, main= "type=1", under=T, extra=1, fallen.leaves = F)
prp(vegtree)
save(vegtree,file = "vegtree.RData")
printcp(vegtree)
plotcp(vegtree)
# return summary output to text file
sink("CARTsummary.txt", append=FALSE, split=FALSE)
summary(vegtree)
sink()

pred = predict(vegtree, type = "class")
confuse <- table(pred, X1$SU)
write.csv(confuse, file= "CART_ConfusionMatrix.csv")

################run a C50 on matrix##########
###Change cp value to show more or less detail
vegtreeC50 <- C5.0(SU ~ ., data = X1, trials=1, rules=F)
##rpart.plot(vegtree, main= "type=1", under=T, extra=1, fallen.leaves = F)
save(vegtreeC50,file = "vegtreeC50.RData")
plot(vegtreeC50, subtree=1)
# return summary output to text file
sink("C5.0summary.txt", append=FALSE, split=FALSE)
summary(vegtreeC50)
sink()

pred = predict(vegtreeC50, newdata = X1, type = "class")
confuse <- table(pred, X1$SU)
write.csv(confuse, file= "CART_ConfusionMatrix.csv")







####run Random Forests on matrix###############
set.seed(71)
vegRF <- randomForest(SU ~ ., data=X1, method= 'class')
print (vegRF$confusion, digits=2)
write.csv(vegRF$confusion, file= "RF_Vegunit_ConfusionMatrix.csv")


# plot tree
plot(vegtree, uniform=TRUE,
     main="Classification Tree for Vegtree")
text(vegtree, use.n=F, all=F, cex=.8)


write.csv(longCart,"longCARTresults.csv",row.names=FALSE)
write.csv(PairT,"Paired.BGC.Table.csv",row.names=FALSE)
# create attractive postscript plot of tree
pdf("tree.pdf",
     title = "Classification Tree for Kyphosis")
dev.off()
################an NMDS of the data############
veg.dis <-dist(X1, method="euclidian")
nms <- isoMDS(veg.dis, y = cmdscale(veg.dis, 2))
plot(nms$points, type = "n")
text(nms$points, labels = as.character(nms$row.names))#as.character(1:nrow(nms)))
nms.sh <- Shepard(veg.dis, nms$points)
plot(nms.sh, pch = ".")
lines(nms.sh$x, nms.sh$yf, type = "S")
####run Random Forests on matrix###############
set.seed(71)
vegRF <- randomForest(SU ~ ., data=X1, method= 'class')
print (vegRF$confusion, digits=2)
write.csv(vegRF$confusion, file= "F:/RF_Vegunit_ConfusionMatrix.csv")


###############predict new data##########################
#####will need to build a function to align the species columns in so the model will predict model will run########
c=file.choose()
vegnew=read.table(c,as.is=TRUE,header=TRUE,stringsAsFactors=FALSE)
setwd(gsub(basename(c),"",c))
colnames(veg)[1:3]=c("PlotNumber", "Species", "Cover")
newvegmtrx <- matrify(vegnew)
Class.pred <- newvegmtrx [, c(1:2)]
colnames(Class.pred)[1:2] = c("PlotNumber", "SU.pred")
Class.pred$SU.pred <- predict(vegtree, newvegmtrx, type="class", na.action=na.omit)
write.csv(Class.pred, "Vegunit_Predicted.csv")

#Output list plot list with Zone prediction
BGC.pred <- Y1
BGC.pred = BGC.pred [ , c (1:2)]
colnames(BGC.pred)[2] = c("Subz")
write.csv(BGC.pred, "Zone_Prediction_by_Plot.csv")






