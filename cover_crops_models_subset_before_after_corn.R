### Egg card data analaysis


covercrop<-read.csv(file="GLBRC_egg_cards_14.csv", header=TRUE, na.strings="NA")

#remove sampling date that was unreplicated in Wisconsin
covercrop<-covercrop[which(covercrop$Sample!="1.5"),]

#subset data to only look at what's going on in open cages
open.only<-covercrop[which(covercrop$Treatment=="Open"),]
attach(open.only)

#binomial model 

#must have an NA free dataset for this to work
#cut NAs out of data
open.only.complete<-na.omit(open.only)

binomial.model.open<-glm(cbind(Number_eggs_removed, Initial_egg_number)~State*Year+Block+Crop*as.factor(Sample), family=binomial(logit), na.action="na.fail", data=open.only.complete)
summary(binomial.model.open)

anova(binomial.model.open, test="Chisq")
summary(anova(binomial.model.open, test="Chisq"))


#Compute a posthoc t-test for the differences between sample-crop interactions
#start by computing a residual dataset with these factors removed from the model
binomial.model.nosample<-glm(cbind(Number_eggs_removed, Initial_egg_number)~State*Year+Block+Crop, family=binomial(logit), na.action="na.fail", data=open.only.complete)
residual.vector.nosample<-binomial.model.nosample$resid
binomial.model.nocrop<-glm(cbind(Number_eggs_removed, Initial_egg_number)~State*Year+Block+as.factor(Sample), family=binomial(logit), na.action="na.fail", data=open.only.complete)
residual.vector.nocrop<-binomial.model.nocrop$resid
binomial.model.noint<-glm(cbind(Number_eggs_removed, Initial_egg_number)~State*Year+Block+Crop+as.factor(Sample), family=binomial(logit), na.action="na.fail", data=open.only.complete)
residual.vector.noint<-binomial.model.noint$resid

#create a vector for the interaction we're trying to get at
crop.sample.interaction<-interaction(open.only.complete[c(5, 3)], drop=TRUE)
attach(open.only.complete)
#run holm-adjusted t-tests
pairwise.t.test(residual.vector.nosample, as.factor(Sample), p.adj = "holm")
pairwise.t.test(residual.vector.nocrop, Crop, p.adj = "holm")
pairwise.t.test(residual.vector.noint, crop.sample.interaction, p.adj = "holm")

#create variable for proportion of eggs eaten
open.only.complete$prop.eaten<-open.only.complete$Number_eggs_removed/open.only.complete$Initial_egg_number
#get means and SEs out
library(plyr)
#year
ddply(open.only.complete, c("Year"), summarise,
      mean = mean(prop.eaten), sd = sd(prop.eaten),
      sem = sd(prop.eaten)/sqrt(length(prop.eaten)))
#state
ddply(open.only.complete, c("State"), summarise,
      mean = mean(prop.eaten), sd = sd(prop.eaten),
      sem = sd(prop.eaten)/sqrt(length(prop.eaten)))
#sample date, where 1= before corn planting, 2= after corn planting, 3= midseason
ddply(open.only.complete, c("Sample"), summarise,
      mean = mean(prop.eaten), sd = sd(prop.eaten),
      sem = sd(prop.eaten)/sqrt(length(prop.eaten)))
#crop
ddply(open.only.complete, c("Crop"), summarise,
      mean = mean(prop.eaten), sd = sd(prop.eaten),
      sem = sd(prop.eaten)/sqrt(length(prop.eaten)))

#sample date by crop
ddply(open.only.complete, c("Crop", "Sample"), summarise,
      mean = mean(prop.eaten), sd = sd(prop.eaten),
      sem = sd(prop.eaten)/sqrt(length(prop.eaten)))

## Vacuum sample analysis

vacuum<-read.csv(file="GLBRC_vacuum_samples_14.csv", header=TRUE, na.strings="NA")
#dataset is roughly the same as egg card data, with no open/closed treatment
#2 response variables Abundance and Simpson diversity "Diversity"
#almost identical approach to above, but we'll model Abundance assuming a Poission
#error distribution and Diversity assuming a normal error distribution
#must have an NA free dataset for this to work
#cut NAs out of data
vacuum.complete<-na.omit(vacuum)
attach(vacuum.complete)

#abundance analysis first

#abundance model 
library(pscl)
abundance.model<-glm.nb(Abundance~State*Year+Block+Crop*as.factor(Sample), na.action="na.fail", data=vacuum.complete)
summary(abundance.model)

anova(abundance.model, test="Chisq")
summary(anova(abundance.model, test="Chisq"))


#Compute a posthoc t-test for the differences between sample-crop interactions, as above
#start by computing a residual dataset with these factors removed from the model
abundance.model.nosample<-glm.nb(Abundance~Year+Block+Crop,na.action="na.fail", data=vacuum.complete)
a.residual.vector.nosample<-abundance.model.nosample$resid
abundance.model.nocrop<-glm.nb(Abundance~Year+Block+as.factor(Sample), na.action="na.fail", data=vacuum.complete)
a.residual.vector.nocrop<-abundance.model.nocrop$resid
abundance.model.noint<-glm.nb(Abundance~Year+Block+Crop+as.factor(Sample), na.action="na.fail", data=vacuum.complete)
a.residual.vector.noint<-abundance.model.noint$resid

#create a vector for the interaction we're trying to get at
a.crop.sample.interaction<-interaction(vacuum.complete[c(5, 3)], drop=TRUE)
attach(vacuum.complete)
#run holm-adjusted t-tests
pairwise.t.test(a.residual.vector.nosample, as.factor(Sample), p.adj = "holm")
pairwise.t.test(a.residual.vector.nocrop, Crop, p.adj = "holm")
pairwise.t.test(a.residual.vector.noint, a.crop.sample.interaction, p.adj = "holm")

#get means and SEs out
library(plyr)
#year
ddply(vacuum.complete, c("Year"), summarise,
      mean = mean(Abundance), sd = sd(Abundance),
      sem = sd(Abundance)/sqrt(length(Abundance)))
#state
ddply(vacuum.complete, c("State"), summarise,
      mean = mean(Abundance), sd = sd(Abundance),
      sem = sd(Abundance)/sqrt(length(Abundance)))
#sample date
ddply(vacuum.complete, c("Sample"), summarise,
      mean = mean(Abundance), sd = sd(Abundance),
      sem = sd(Abundance)/sqrt(length(Abundance)))
#crop
ddply(vacuum.complete, c("Crop"), summarise,
      mean = mean(Abundance), sd = sd(Abundance),
      sem = sd(Abundance)/sqrt(length(Abundance)))

#sample date by crop
ddply(vacuum.complete, c("Crop", "Sample"), summarise,
      mean = mean(Abundance), sd = sd(Abundance),
      sem = sd(Abundance)/sqrt(length(Abundance)))

#Diversity model next

#diversity model 
diversity.model<-glm(Diversity~State*Year+Block+Crop*as.factor(Sample), na.action="na.fail", data=vacuum.complete)
summary(diversity.model)

anova(diversity.model, test="Chisq")
summary(anova(diversity.model, test="Chisq"))


#Compute a posthoc t-test for the differences between sample-crop interactions, as above
#start by computing a residual dataset with these factors removed from the model
diversity.model.nosample<-glm(Diversity~Year+State+Crop, na.action="na.fail", data=vacuum.complete)
a.residual.vector.nosample<-diversity.model.nosample$resid
diversity.model.nocrop<-glm(Diversity~Year+State+as.factor(Sample), na.action="na.fail", data=vacuum.complete)
a.residual.vector.nocrop<-diversity.model.nocrop$resid
diversity.model.noint<-glm(Diversity~Year+State+Crop+as.factor(Sample), na.action="na.fail", data=vacuum.complete)
a.residual.vector.noint<-diversity.model.noint$resid


#run holm-adjusted t-tests
pairwise.t.test(a.residual.vector.nosample, as.factor(Sample), p.adj = "holm")
pairwise.t.test(a.residual.vector.nocrop, Crop, p.adj = "holm")
pairwise.t.test(a.residual.vector.noint, a.crop.sample.interaction, p.adj = "holm")

#get means and SEs out
library(plyr)
#year
ddply(vacuum.complete, c("Year"), summarise,
      mean = mean(Diversity), sd = sd(Diversity),
      sem = sd(Diversity)/sqrt(length(Diversity)))
#state
ddply(vacuum.complete, c("State"), summarise,
      mean = mean(Diversity), sd = sd(Diversity),
      sem = sd(Diversity)/sqrt(length(Diversity)))
#sample date
ddply(vacuum.complete, c("Sample"), summarise,
      mean = mean(Diversity), sd = sd(Diversity),
      sem = sd(Diversity)/sqrt(length(Diversity)))
#crop
ddply(vacuum.complete, c("Crop"), summarise,
      mean = mean(Diversity), sd = sd(Diversity),
      sem = sd(Diversity)/sqrt(length(Diversity)))

#sample date by crop
ddply(vacuum.complete, c("Crop", "Sample"), summarise,
      mean = mean(Diversity), sd = sd(Diversity),
      sem = sd(Diversity)/sqrt(length(Diversity)))


#Analysis of community data
community<-read.csv(file="GLBRC_community_matrix_14.csv", header=TRUE, na.strings="NA")

#create dataset without NAs for analysis
community.complete<-community[complete.cases(community),]

#need to total captures up by treatmentX sampling date to take care of zero bias in individual samples
library(reshape2)
melted<-melt(community.complete, id=1:6, na.rm=TRUE)
community.by.block<-dcast(melted, State+Year+Sample+Crop~variable, sum)

#also remove samples that captured one or fewer insects
com.by.block.1<-community.by.block[rowSums(community.by.block[5:34])>1,]

#rename rows to provide 
row.names(com.by.block.1) <- NULL 

#Create matrix of environemntal variables
env.matrix<-com.by.block.1[c(1:4)]
#create matrix of community variables
com.matrix<-com.by.block.1[c(5:34)]
#delete columns of community matrix where non-predators were recorded
com.matrix<-subset(com.matrix, select=-c(Crickets,Slugs))
#delete columns from community matrix that have 2 or fewer observations
com.matrix.1<-com.matrix[,colSums(com.matrix)>2]

#transform data upfront so all tests will be done using same data
#com.matrix.1<-wisconsin(com.matrix.1)

library(vegan)
#NMDS of communitiy
ord<-metaMDS(com.matrix.1, autotransform=FALSE)
ord
most_abund<-colSums(com.matrix.1)>100
plot(ord, disp='sites', type="n")

#display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="1"), pch=19, col="black")
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="2"), pch=15, col="black")
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="1"), pch=21, col="black")
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="2"), pch=22, col="black")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="lightblue", kind="se", conf=0.95, label=TRUE)
legend("bottomright", title=NULL, pch=c(19,15,21,22,19), col=c("black","black","black","black", "lightblue"), cex=0.75, legend=c("WI 2013", "WI 2014", "MI 2013", "MI 2014", "Crop"))
#points(ord, display="species", select=which(most_abund==FALSE), pch=21, cex=1, col="red")
#ordilabel(ord, display="species", select=which(most_abund==TRUE), cex=0.75, col="red", fill="white")



#overlay environmental variables
ordfit<-envfit(ord~State+as.factor(Year)+as.factor(Sample)+Crop, data=env.matrix, perm=1000)
summary(ordfit)
ordfit

#try a permanona
mod1<-adonis(com.matrix.1~State+as.factor(Year)+as.factor(Sample)*Crop, data=env.matrix, method="bray")
mod1

#because adonis doesn't have functionality for posthoc tests, we will need to do pairwise comparisons 
#by subsetting the data into each combination of crops and re-running the analyis.

#corn-cover
com.matrix.corn.cover<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="COVER"),]
env.matrix.corn.cover<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="COVER"),]

mod1<-adonis(com.matrix.corn.cover~State+as.factor(Year)+as.factor(Sample)*Crop, data=env.matrix.corn.cover, method="bray")
mod1

#corn-switch
com.matrix.corn.switch<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="SWITCH"),]
env.matrix.corn.switch<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="SWITCH"),]

mod1<-adonis(com.matrix.corn.switch~State+as.factor(Year)+as.factor(Sample)*Crop, data=env.matrix.corn.switch, method="bray")
mod1

#corn-prairie
com.matrix.corn.prairie<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.corn.prairie<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.corn.prairie~State+as.factor(Year)+as.factor(Sample)*Crop, data=env.matrix.corn.prairie, method="bray")
mod1

#cover-switch
com.matrix.cover.switch<-com.matrix.1[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="SWITCH"),]
env.matrix.cover.switch<-env.matrix[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="SWITCH"),]

mod1<-adonis(com.matrix.cover.switch~State+as.factor(Year)+as.factor(Sample)*Crop, data=env.matrix.cover.switch, method="bray")
mod1

#cover-prairie
com.matrix.cover.prairie<-com.matrix.1[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.cover.prairie<-env.matrix[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.cover.prairie~State+as.factor(Year)+as.factor(Sample)*Crop, data=env.matrix.cover.prairie, method="bray")
mod1

#switch-prairie
com.matrix.switch.prairie<-com.matrix.1[which(env.matrix$Crop=="SWITCH"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.switch.prairie<-env.matrix[which(env.matrix$Crop=="SWITCH"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.switch.prairie~State+as.factor(Year)+as.factor(Sample)*Crop, data=env.matrix.switch.prairie, method="bray")
mod1

#repeat analysis for before and after corn planting

#create data frames for subsets of data taken before (sample 1-3) or after (sample 4,5) corn planting
com.by.block.before<-com.by.block.1[which(as.numeric(com.by.block.1$Sample)<=3),]
com.by.block.after<-com.by.block.1[which(as.numeric(com.by.block.1$Sample)>=4),]


#before corn planting
#Create matrix of environemntal variables
env.matrix<-com.by.block.before[c(1:4)]
#create matrix of community variables
com.matrix<- com.by.block.before [c(5:34)]
#delete columns of community matrix where non-predators were recorded
com.matrix<-subset(com.matrix, select=-c(Crickets,Slugs))
#delete columns from community matrix that have 2 or fewer observations
com.matrix.1<-com.matrix[,colSums(com.matrix)>2]

#transform data upfront so all tests will be done using same data
#com.matrix.1<-wisconsin(com.matrix.1)

#NMDS of community before planting corn
ord<-metaMDS(com.matrix.1, autotransform=FALSE)
ord
most_abund<-colSums(com.matrix.1)>100
plot(ord, disp='sites', type="n")
title(main="A", cex.main=2, adj=0)
#ordiellipse(ord, env.matrix$Sample, draw="polygon",  col="brown", kind="se", conf=0.95)
#display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="1"), pch=19, col="black")
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="2"), pch=15, col="black")
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="1"), pch=21, col="black")
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="2"), pch=22, col="black")
#levels(env.matrix$Crop)=c("MAIZE","COVER","PRAIRIE","SWITCH")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="lightblue", kind="se", conf=0.95, label=TRUE)
legend("bottomright", title=NULL, pch=c(19,15,21,22,19), col=c("black","black","black","black", "lightblue"), cex=0.75, legend=c("WI 2013", "WI 2014", "MI 2013", "MI 2014", "Crop"))
#points(ord, display="species", select=which(most_abund==FALSE), pch=21, cex=1, col="red")
#ordilabel(ord, display="species", select=which(most_abund==TRUE), cex=0.75, col="red", fill="white")



#try a permanona
mod1<-adonis(com.matrix.1~State+as.factor(Year)+Crop, data=env.matrix, method="bray")
mod1

#because adonis doesn't have functionality for posthoc tests, we will need to do pairwise comparisons 
#by subsetting the data into each combination of crops and re-running the analyis.

#corn-cover
com.matrix.corn.cover<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="COVER"),]
env.matrix.corn.cover<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="COVER"),]

mod1<-adonis(com.matrix.corn.cover~State+as.factor(Year)+Crop, data=env.matrix.corn.cover, method="bray")
mod1

#corn-switch
com.matrix.corn.switch<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="SWITCH"),]
env.matrix.corn.switch<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="SWITCH"),]

mod1<-adonis(com.matrix.corn.switch~State+as.factor(Year)+Crop, data=env.matrix.corn.switch, method="bray")
mod1

#corn-prairie
com.matrix.corn.prairie<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.corn.prairie<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.corn.prairie~State+as.factor(Year)+Crop, data=env.matrix.corn.prairie, method="bray")
mod1

#cover-switch
com.matrix.cover.switch<-com.matrix.1[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="SWITCH"),]
env.matrix.cover.switch<-env.matrix[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="SWITCH"),]

mod1<-adonis(com.matrix.cover.switch~State+as.factor(Year)+Crop, data=env.matrix.cover.switch, method="bray")
mod1

#cover-prairie
com.matrix.cover.prairie<-com.matrix.1[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.cover.prairie<-env.matrix[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.cover.prairie~State+as.factor(Year)+Crop, data=env.matrix.cover.prairie, method="bray")
mod1

#switch-prairie
com.matrix.switch.prairie<-com.matrix.1[which(env.matrix$Crop=="SWITCH"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.switch.prairie<-env.matrix[which(env.matrix$Crop=="SWITCH"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.switch.prairie~State+as.factor(Year)+Crop, data=env.matrix.switch.prairie, method="bray")
mod1


#after corn planting
#Create matrix of environemntal variables
env.matrix<-com.by.block.after[c(1:4)]
#create matrix of community variables
com.matrix<- com.by.block.after [c(5:34)]
#delete columns of community matrix where non-predators were recorded
com.matrix<-subset(com.matrix, select=-c(Crickets,Slugs))
#delete columns from community matrix that have 2 or fewer observations
com.matrix.1<-com.matrix[,colSums(com.matrix)>2]

#transform data upfront so all tests will be done using same data
#com.matrix.1<-wisconsin(com.matrix.1)
#NMDS of community after planting corn
ord<-metaMDS(com.matrix.1, autotransform=FALSE)
ord
most_abund<-colSums(com.matrix.1)>100
plot(ord, disp='sites', type="n")
title(main="B", cex.main=2, adj=0)
#ordiellipse(ord, env.matrix$Sample, draw="polygon",  col="brown", kind="se", conf=0.95)
#display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="1"), pch=19, col="black")
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="2"), pch=15, col="black")
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="1"), pch=21, col="black")
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="2"), pch=22, col="black")
#levels(env.matrix$Crop)=c("MAIZE","COVER","PRAIRIE","SWITCH")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="lightblue", kind="se", conf=0.95, label=TRUE)
legend("bottomright", title=NULL, pch=c(19,15,21,22,19), col=c("black","black","black","black", "lightblue"), cex=0.75, legend=c("WI 2013", "WI 2014", "MI 2013", "MI 2014", "Crop"))
#points(ord, display="species", select=which(most_abund==FALSE), pch=21, cex=1, col="red")
#ordilabel(ord, display="species", select=which(most_abund==TRUE), cex=0.75, col="red", fill="white")

#try a permanona
mod1<-adonis(com.matrix.1~State+as.factor(Year)+Crop, data=env.matrix, method="bray")
mod1


#because adonis doesn't have functionality for posthoc tests, we will need to do pairwise comparisons 
#by subsetting the data into each combination of crops and re-running the analyis.

#corn-cover
com.matrix.corn.cover<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="COVER"),]
env.matrix.corn.cover<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="COVER"),]

mod1<-adonis(com.matrix.corn.cover~State+as.factor(Year)+Crop, data=env.matrix.corn.cover, method="bray")
mod1

#corn-switch
com.matrix.corn.switch<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="SWITCH"),]
env.matrix.corn.switch<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="SWITCH"),]

mod1<-adonis(com.matrix.corn.switch~State+as.factor(Year)+Crop, data=env.matrix.corn.switch, method="bray")
mod1

#corn-prairie
com.matrix.corn.prairie<-com.matrix.1[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.corn.prairie<-env.matrix[which(env.matrix$Crop=="CORN"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.corn.prairie~State+as.factor(Year)+Crop, data=env.matrix.corn.prairie, method="bray")
mod1

#cover-switch
com.matrix.cover.switch<-com.matrix.1[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="SWITCH"),]
env.matrix.cover.switch<-env.matrix[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="SWITCH"),]

mod1<-adonis(com.matrix.cover.switch~State+as.factor(Year)+Crop, data=env.matrix.cover.switch, method="bray")
mod1

#cover-prairie
com.matrix.cover.prairie<-com.matrix.1[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.cover.prairie<-env.matrix[which(env.matrix$Crop=="COVER"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.cover.prairie~State+as.factor(Year)+Crop, data=env.matrix.cover.prairie, method="bray")
mod1

#switch-prairie
com.matrix.switch.prairie<-com.matrix.1[which(env.matrix$Crop=="SWITCH"|env.matrix$Crop=="PRAIRIE"),]
env.matrix.switch.prairie<-env.matrix[which(env.matrix$Crop=="SWITCH"|env.matrix$Crop=="PRAIRIE"),]

mod1<-adonis(com.matrix.switch.prairie~State+as.factor(Year)+Crop, data=env.matrix.switch.prairie, method="bray")
mod1
## aphid predator exclusion

aphids<-read.csv(file="GLBRC_aphids_14.csv", header=TRUE, na.strings="NA")

library(pscl)

#assume negative binomial distribution for aphid data
aphid.model<-glm.nb(aphids~Day*Treatment+as.factor(Block)+offset(tillers), na.action="na.fail", data=aphids)
summary(aphid.model)

anova(aphid.model, test="Chisq")
summary(anova(aphid.model, test="Chisq"))

#Compute a posthoc t-test for the differences between day-treatment interactions

aphid.interaction<-interaction(aphids[c(1, 4)], drop=TRUE)

#run holm-adjusted t-test- don't know if this is useful- we just really need the treatment differences and the ANOVA does that
pairwise.t.test(aphids$aphids, aphid.interaction, p.adj = "holm")


#get means and SEs out
library(plyr)
#day
ddply(aphids, c("Day"), summarise,
      mean = mean(aphids), sd = sd(aphids),
      sem = sd(aphids)/sqrt(length(aphids)))
#treatment
ddply(aphids, c("Treatment"), summarise,
      mean = mean(aphids), sd = sd(aphids),
      sem = sd(aphids)/sqrt(length(aphids)))
#day by treatment
ddply(aphids, c("Day", "Treatment"), summarise,
      mean = mean(aphids), sd = sd(aphids),
      sem = sd(aphids)/sqrt(length(aphids)))



#####################

#figures

#####################


library(reshape2)
library(ggplot2)
library(gridExtra)


###
#figure 1 NMDS
###

#set up two panel vertical stack design

par(mfrow=c(2,1))

#panel A
#before corn planting


#Create matrix of environemntal variables
env.matrix<-com.by.block.before[c(1:4)]
#create matrix of community variables
com.matrix<-com.by.block.before [c(5:34)]
#delete columns of community matrix where non-predators were recorded
com.matrix<-subset(com.matrix, select=-c(Crickets,Slugs))
#delete columns from community matrix that have 2 or fewer observations
com.matrix.1<-com.matrix[,colSums(com.matrix)>2]

#NMDS of community before planting corn
ord<-metaMDS(com.matrix.1, autotransform=FALSE)
ord
most_abund<-colSums(com.matrix.1)>100
plot(ord, disp='sites', type="n")
title(main="A", cex.main=1.5, adj=0)
#display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="1"), pch=19, col="black", cex=0.75)
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="2"), pch=15, col="black", cex=0.75)
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="1"), pch=21, col="black", cex=0.75)
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="2"), pch=22, col="black", cex=0.75)
levels(env.matrix$Crop)=c("MAIZE","COVER","PRAIRIE","SWITCH")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#a6cee3", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="MAIZE")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#1f78b4", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="COVER")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#b2df8a", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="SWITCH")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#33a02c", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="PRAIRIE")
legend(1.65,0.1, title=NULL, pch=c(19,15,21,22,19,19,19,19), ncol=1, text.width=0.6,
      col=c("black","black","black","black", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"), 
      cex=0.5, legend=c("WI 2013", "WI 2014", "MI 2013", "MI 2014", "Maize", "Cover", "Switch", "Prairie"))

#panel B
#after corn planting
#Create matrix of environemntal variables
env.matrix<-com.by.block.after[c(1:4)]
#create matrix of community variables
com.matrix<- com.by.block.after [c(5:34)]
#delete columns of community matrix where non-predators were recorded
com.matrix<-subset(com.matrix, select=-c(Crickets,Slugs))
#delete columns from community matrix that have 2 or fewer observations
com.matrix.1<-com.matrix[,colSums(com.matrix)>2]

#NMDS of community after planting corn
ord<-metaMDS(com.matrix.1, autotransform=FALSE)
ord
most_abund<-colSums(com.matrix.1)>100
plot(ord, disp='sites', type="n", cex=0.5)
title(main="B", cex.main=1.5, adj=0)
#display WI data as solid shapes, MI as outlines, and 2013 as circles, 2014 as squares
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="1"), pch=19, col="black", cex=0.75)
points(ord, display="sites", select=which(env.matrix$State=="WI"& env.matrix$Year=="2"), pch=15, col="black", cex=0.75)
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="1"), pch=21, col="black", cex=0.75)
points(ord, display="sites", select=which(env.matrix$State=="MI"& env.matrix$Year=="2"), pch=22, col="black", cex=0.75)
levels(env.matrix$Crop)=c("MAIZE","COVER","PRAIRIE","SWITCH")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#a6cee3", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="MAIZE")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#1f78b4", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="COVER")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#b2df8a", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="SWITCH")
ordiellipse(ord, env.matrix$Crop, draw="polygon", col="#33a02c", kind="se", conf=0.95, label=FALSE, cex=0.75, show.groups="PRAIRIE")


###
# Bar graphs
###

#set common graphical parameters

#labels
#original colours:
#col.vec<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")

#greyscale
col.vec<-c("white", "gray83", "gray68", "gray43")
sample.vec<-c("Pre Maize Week 1", "Pre Maize Week 2", "Pre Maize Week 3", "Post Maize Planting", "Mid Season")
egg.sample.vec<-c("Pre Maize planting",  "Post Maize Planting", "Mid Season")

#assign ordinal values to categories so they all plot in the same order
vacuum.complete$Crop<-factor(vacuum.complete$Crop, levels=c("CORN","COVER","SWITCH","PRAIRIE"), labels=c("Continuous Maize", "Cover Crop System", "Switchgrass", "Prairie"))
vacuum.complete$Sample<-factor(vacuum.complete$Sample, levels=c("1","2","3","4","5"), labels=sample.vec)
#egg cards
open.only.complete$Crop<-factor(open.only.complete$Crop, levels=c("CORN","COVER","SWITCH","PRAIRIE"), labels=c("Continuous Maize", "Cover Crop System", "Switchgrass", "Prairie"))
open.only.complete$Sample<-factor(open.only.complete$Sample, levels=c("1","2","3"), labels=egg.sample.vec)


###
#figure 2 Abundance
###



#Panel A

# create summary data frame
abundance.by.crop<- ddply(vacuum.complete, c("Crop"), summarise,
                      mean = mean(Abundance), sd = sd(Abundance),
                      sem = sd(Abundance)/sqrt(length(Abundance)))

fig2A<-ggplot(abundance.by.crop, aes(y=mean, x=Crop, colour=Crop))+
  geom_bar(stat="identity", fill=col.vec, colour="black", # Use black outlines,
           size=0.5) +
  scale_fill_manual(values=Crop)+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                size=0.3,    # Thinner lines
                width=0.2,
                position=position_dodge(.9),
                colour="black") +
  xlab("Crop") +
  ylab("Mean abundance per sample") +
  ggtitle("A\n")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0, size=22),
        panel.background=element_rect(colour="black"), 
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))+
  annotate("text", x=1, y=4, label="C")+ #maize grouping
  annotate("text", x=2, y=5, label="C")+ #cover grouping
  annotate("text", x=3, y=9, label="B")+ #switch grouping
  annotate("text", x=4, y=18, label="A   ") #prairie grouping

fig2A

#Panel B

# create summary data frame
abundance.by.crop.date<- ddply(vacuum.complete, c("Crop", "Sample"), summarise,
                          mean = mean(Abundance), sd = sd(Abundance),
                          sem = sd(Abundance)/sqrt(length(Abundance)))

fig2B<-ggplot(abundance.by.crop.date, aes(y=mean, x=Sample, fill=Crop))+
  geom_bar(position=position_dodge(), stat="identity", colour="black", aes(fill=Crop), # Use black outlines,
           size=0.5) +
  scale_fill_manual(values=col.vec)+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                size=0.3,    # Thinner lines
                position="dodge",
                colour="black") +
  xlab("Date") +
  ylab(NULL) +
  ggtitle("B\n")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0, size=22),
        panel.background=element_rect(colour="black"), 
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")

fig2B

#stack figures together

grid.arrange(arrangeGrob(fig2A,fig2B, ncol=2, widths=c(1,2)))

###
#figure 3 Diversity
###



#Panel A

# create summary data frame
diversity.by.crop<- ddply(vacuum.complete, c("Crop"), summarise,
                          mean = mean(Diversity), sd = sd(Diversity),
                          sem = sd(Diversity)/sqrt(length(Diversity)))

fig3A<-ggplot(diversity.by.crop, aes(y=mean, x=Crop, colour=Crop))+
  geom_bar(stat="identity", fill=col.vec, colour="black", # Use black outlines,
           size=0.5) +
  scale_fill_manual(values=Crop)+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                size=0.3,    # Thinner lines
                width=0.2,
                position=position_dodge(.9),
                colour="black") +
  xlab("Crop") +
  ylab("Mean diversity per sample") +
  ggtitle("A\n")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0, size=22),
        panel.background=element_rect(colour="black"), 
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))+
  annotate("text", x=1, y=0.25, label="B")+ #maize grouping
  annotate("text", x=2, y=0.3, label="B")+ #cover grouping
  annotate("text", x=3, y=0.5, label="A")+ #switch grouping
  annotate("text", x=4, y=0.5, label="A") #prairie grouping

fig3A

#Panel B

# create summary data frame
diversity.by.crop.date<- ddply(vacuum.complete, c("Crop", "Sample"), summarise,
                               mean = mean(Diversity), sd = sd(Diversity),
                               sem = sd(Diversity)/sqrt(length(Diversity)))

fig3B<-ggplot(diversity.by.crop.date, aes(y=mean, x=Sample, fill=Crop))+
  geom_bar(position=position_dodge(), stat="identity", colour="black", aes(fill=Crop), # Use black outlines,
           size=0.5) +
  scale_fill_manual(values=col.vec)+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                size=0.3,    # Thinner lines
                position="dodge",
                colour="black") +
  xlab("Date") +
  ylab(NULL) +
  ggtitle("B\n")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0, size=22),
        panel.background=element_rect(colour="black"), 
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")

fig3B

#stack figures together

grid.arrange(arrangeGrob(fig3A,fig3B, ncol=2, widths=c(1,2)))

###
#figure 4 Egg cards
###



#Panel A

# create summary data frame
biocontrol.by.crop<- ddply(open.only.complete, c("Crop"), summarise,
                          mean = mean(prop.eaten), sd = sd(prop.eaten),
                          sem = sd(prop.eaten)/sqrt(length(prop.eaten)))

fig4A<-ggplot(biocontrol.by.crop, aes(y=mean, x=Crop, colour=Crop))+
  geom_bar(stat="identity", fill=col.vec, colour="black", # Use black outlines,
           size=0.5) +
  scale_fill_manual(values=Crop)+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                size=0.3,    # Thinner lines
                width=0.2,
                position=position_dodge(.9),
                colour="black") +
  xlab("Crop") +
  ylab("Mean Percent Removal in 48h") +
  ggtitle("A\n")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0, size=22),
        panel.background=element_rect(colour="black"), 
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))+
  annotate("text", x=1, y=0.27, label="B")+ #maize grouping
  annotate("text", x=2, y=0.27, label="B")+ #cover grouping
  annotate("text", x=3, y=0.59, label="A")+ #switch grouping
  annotate("text", x=4, y=0.54, label="A") #prairie grouping

fig4A

#Panel B

# create summary data frame
biocontrol.by.crop.date<- ddply(open.only.complete, c("Crop", "Sample"), summarise,
                               mean = mean(prop.eaten), sd = sd(prop.eaten),
                               sem = sd(prop.eaten)/sqrt(length(prop.eaten)))

fig4B<-ggplot(biocontrol.by.crop.date, aes(y=mean, x=Sample, fill=Crop))+
  geom_bar(position=position_dodge(), stat="identity", colour="black", aes(fill=Crop), # Use black outlines,
           size=0.5) +
  scale_fill_manual(values=col.vec)+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                size=0.3,    # Thinner lines
                position="dodge",
                colour="black") +
  xlab("Date") +
  ylab(NULL) +
  ggtitle("B\n")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0, size=22),
        panel.background=element_rect(colour="black"), 
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")

fig4B

#stack figures together

grid.arrange(arrangeGrob(fig4A,fig4B, ncol=2, widths=c(1,2)))

###
#figure 5 time series
###

timeseries<-ddply(aphids, c("Day", "Treatment"), summarise,
        mean = mean(aphids), sd = sd(aphids),
        sem = sd(aphids)/sqrt(length(aphids)))

fig5<-ggplot(timeseries, aes(y=mean, x=Day, group=Treatment, shape=Treatment))+
      geom_line(aes(linetype=Treatment), size=1)+
      geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
              size=0.3, 
              width=0.2,
              colour="black") +
      geom_point(size=5, fill="white")+
      scale_shape_manual(values=c(22,21))+
      xlab("Day") +
      ylab("Mean Aphids per Cage") +
      theme_bw()+
      theme(panel.background=element_rect(colour="black"), 
        panel.grid=element_blank(),
        axis.text.x=element_text(hjust=1),
        legend.position=c(0.1, 0.8),
        legend.key=element_blank())

fig5
