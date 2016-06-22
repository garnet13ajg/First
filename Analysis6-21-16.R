#load packages with command library() with the name of the package between the parentheses
#packages I may need: pastecs, vegan, simba, ecodist, cluster, MASS

#Import file speabu4-5.csv
spedata <- read.csv('presabs.csv',header=TRUE,row.names=1)

#Examine structure of data set
str(spedata)

#Summary statistics # nbr.null is number of objects in which variable equals zero
stat.desc(spedata)

#In how many sites do species occur and in what abundance to species occur?
#ECDF or histogram seems to display data best and is easiest to interpret
foa.plots(spedata)

#Empirical cumulative distribution function (ECDF) - rank order distribution of increasing values of variable
ecdf.plots(spedata)

#Or you can display information as a histogram
hist.plots(spedata)

#Or as a box plot
box.plots(spedata)

#Or normal quantile-quantil (QQ) plot compares ECDF to expected ECDF for normal distributed plot
qqnorm.plots(spedata)

#Species Accumulation curve (vegan)
specaccum(spedata, method="exact", permutations=100,conditioned=TRUE, gamma="jack1")

#Species Accumulation curve for each basin
adm <- read.csv('admiralty.csv',header=TRUE,row.names=1)
cen <- read.csv('central.csv',header=TRUE,row.names=1)
sou <- read.csv('south.csv',header=TRUE,row.names=1)
whi <- read.csv('whidbey.csv',header=TRUE,row.names=1)
ros <- read.csv('rosario.csv',header=TRUE,row.names=1)
hoo <- read.csv('hood.csv',header=TRUE,row.names=1)

#Method "exact" gives expected SA,"collector" gives observed richness
specaccum(ros, method="collector", permutations=100,conditioned=TRUE, gamma="jack1")
rossa<-specaccum(ros)
plot(rossa, add = FALSE, ci = 2, ci.type = c("bar", "line", "polygon"), 
     col = par("fg"), ci.col = col, ci.lty = 1, main = "Rosario",ylab=NA, xlab=NA,ylim=c(0,50), cex.main=1.5, cex.axis=1.5)


specaccum(whi, method="collector", permutations=100, gamma="jack1")
whisa<-specaccum(whi)
plot(whisa, add = FALSE, ci = 2, ci.type = c("bar", "line", "polygon"), 
     col = par("fg"), ci.col = col, ci.lty = 1, main = "Whidbey",ylab=NA, xlab=NA, ylim=c(0,50), cex.main=1.5, cex.axis=1.5)

specaccum(adm, method="collector", permutations=100,conditioned=TRUE, gamma="jack1")
admsa<-specaccum(adm)
plot(admsa, add = FALSE, ci = 2, ci.type = c("bar", "line", "polygon"), 
     col = par("fg"), ci.col = col, ci.lty = 1, main = "Admiralty",ylab=NA, xlab=NA, ylim=c(0,50), cex.main=1.5, cex.axis=1.5)


specaccum(hoo, method="collector", permutations=100,conditioned=TRUE, gamma="jack1")
hoosa<-specaccum(hoo)
plot(hoosa, add = FALSE, ci = 2, ci.type = c("bar", "line", "polygon"), 
     col = par("fg"), ci.col = col, ci.lty = 1, main = "Hood Canal",ylab=NA, xlab=NA, ylim=c(0,50), cex.main=1.5, cex.axis=1.5)

specaccum(cen, method="collector", permutations=100,conditioned=TRUE, gamma="jack1")
censa<-specaccum(cen)
plot(censa, add = FALSE, ci = 2, ci.type = c("bar", "line", "polygon"), 
     col = par("fg"), ci.col = col, ci.lty = 1, main = "Central",ylab=NA, xlab=NA, ylim=c(0,50), cex.main=1.5, cex.axis=1.5)

specaccum(sou, method="collector", permutations=100,conditioned=TRUE, gamma="jack1")
sousa<-specaccum(sou)
plot(sousa, add = FALSE, ci = 2, ci.type = c("bar", "line", "polygon"), 
     col = par("fg"), ci.col = col, ci.lty = 1, main = "South",ylab=NA, xlab=NA, ylim=c(0,50), cex.main=1.5, cex.axis=1.5)



#Combine multiple plots, oma allows edits of outer margin, mtext adds labels to outer margin, cex allows change in text size
par(mfrow=c(2,3), oma=c(2,2,2,2))
mtext("Samples", side = 1, outer=TRUE, cex=2) 
mtext("Richness", side = 2, outer=TRUE, cex=2)
mtext("Species Accumulation Curve", side = 3, outer=TRUE, cex=2)


#Log transform and save data
logdata <- data.trans(spedata,method='log')

#Transform to pres/abs data, plot=F so you don't hav eto press enter for each plot
presabs <- data.trans(spedata,method='power',exp=0, plot=F)
save(presabs, file="presabs")
write.csv(presabs, file = "presabs.csv")
presabs <- read.csv('presabs.csv',header=TRUE,row.names=1)

#Calculate similarity among sites according to fish community composition based on pres/abs
#Sim function in simba library to calculate series of (dis)similarity coefficients
#Calculate jaccard's similarity matrix (assymetrical association coefficient)
jacsim <- sim(presabs, method="jaccard", plot=F)

#Species specific univariate analysis
#Logistic regression
Amodel1 <- glm(Amhex~1, family=binomial, data=presabs)
Amodel2 <- glm(Amhex~Basin, family=binomial, data=presabs)
Amodel3 <- glm(Amhex~Month, family=binomial, data=presabs)
Amodel4 <- glm(Amhex~Basin + Month, family=binomial, data=presabs)
Amodel5 <- glm(Amhex~Basin + Month + Basin:Month, family=binomial, data=presabs)

Amodel1.AIC <- AIC(Amodel1)
Amodel2.AIC <- AIC(Amodel2)
Amodel3.AIC <- AIC(Amodel3)
Amodel4.AIC <- AIC(Amodel4)
Amodel5.AIC <- AIC(Amodel5)

Pmodel1 <- glm(Platstel~1, family=binomial, data=presabs)
Pmodel2 <- glm(Platstel~Basin, family=binomial, data=presabs)
Pmodel3 <- glm(Platstel~Month, family=binomial, data=presabs)
Pmodel4 <- glm(Platstel~Basin + Month, family=binomial, data=presabs)
Pmodel5 <- glm(Platstel~Basin + Month + Basin:Month, family=binomial, data=presabs )

Pmodel1.AIC <- AIC(Pmodel1)
Pmodel2.AIC <- AIC(Pmodel2)
Pmodel3.AIC <- AIC(Pmodel3)
Pmodel4.AIC <- AIC(Pmodel4)
Pmodel5.AIC <- AIC(Pmodel5)

Cmodel1 <- glm(Clupal~1, family=binomial, data=presabs)
Cmodel2 <- glm(Clupal~Basin, family=binomial, data=presabs)
Cmodel3 <- glm(Clupal~Month, family=binomial, data=presabs)
Cmodel4 <- glm(Clupal~Basin + Month, family=binomial, data=presabs)
Cmodel5 <- glm(Clupal~Basin + Month +Basin:Month, family=binomial, data=presabs)

Cmodel1.AIC <- AIC(Cmodel1)
Cmodel2.AIC <- AIC(Cmodel2)
Cmodel3.AIC <- AIC(Cmodel3)
Cmodel4.AIC <- AIC(Cmodel4)
Cmodel5.AIC <- AIC(Cmodel5)
