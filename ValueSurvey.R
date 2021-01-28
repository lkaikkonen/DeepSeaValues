##### VALUES FOR REMOTE ENVIRONMENTS ####

#### Survey Results 2019-10 ####

install.packages("psych")
install.packages(c("foreign", "survey", "knitr"))
install.packages("dplyr")
install.packages("tidyverse","cluster","factoextra","dendextend")

library(psych)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(dplyr)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend)

source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "RBGL", "Rgraphviz"))


data1<- read.table("\FINALresults_added_4s.txt", header=TRUE,fill=TRUE, sep='\t')


## Prune & subset data

df <- data1[,4:99]
PVQ<-df[,65:81]
NEP<-df[,82:96]

# NEP score

df$NEPSc<-rowMeans(df[,c("NEP1", "NEP2", "NEP3", "NEP4", "NEP5", "NEP6", "NEP7", "NEP8", "NEP9", "NEP10", "NEP11", "NEP12", "NEP13", "NEP14", "NEP15")], na.rm = TRUE)

## PVQ

Bio<-df[,c("PVQ_Bio1", "PVQ_Bio2","PVQ_Bio3","PVQ_Bio4" )] #Biospheric values
df$BioV<-rowMeans(df[,c("PVQ_Bio1", "PVQ_Bio2","PVQ_Bio3","PVQ_Bio4")], na.rm = TRUE)

Alt<-df[,c( "PVQ_Alt1", "PVQ_Alt2","PVQ_Alt3","PVQ_Alt4", "PVQ_Alt5" )]#Altruistic
df$AltV<-rowMeans(df[,c("PVQ_Alt1", "PVQ_Alt2","PVQ_Alt3","PVQ_Alt4", "PVQ_Alt5" )], na.rm = TRUE)

Hed<-df[,c( "PVQ_Hed1","PVQ_Hed2","PVQ_Hed3" )]#Hedonic
df$HedV<-rowMeans(df[,c( "PVQ_Hed1","PVQ_Hed2","PVQ_Hed3" )], na.rm = TRUE)

Ego<-df[,c( "PVQ_Ego1", "PVQ_Ego2","PVQ_Ego3", "PVQ_Ego4","PVQ_Ego5")] #Egoistic
df$EgoV<-rowMeans(df[,c("PVQ_Ego1", "PVQ_Ego2","PVQ_Ego3", "PVQ_Ego4","PVQ_Ego5" )], na.rm = TRUE)

### SYMBOLIC VALUES

# Convert 1-7 scale to -3 to 3

df[,5:12]<-df[,5:12]-4
df[,16:23]<-df[,16:23]-4
df[,27:34]<-df[,27:34]-4
df[,38:45]<-df[,38:45]-4

## Symb value scores

df$DS_symb<-rowMeans(df[,c("DS_Beauty", "DS_Mystical", "DS_Importance", "DS_Abundance", "DS_Exciting", "DS_Inviting", "DS_Relaxing", "DS_Calm")])
df$Ant_symb<-rowMeans(df[,c("Ant_Beauty", "Ant_Mystical", "Ant_Importance", "Ant_Abundance", "Ant_Exciting", "Ant_Inviting", "Ant_Relaxing", "Ant_Calm")])
df$RT_symb<-rowMeans(df[,c("RT_Beauty", "RT_Mystical", "RT_Importance", "RT_Abundance", "RT_Exciting", "RT_Inviting", "RT_Relaxing", "RT_Calm")])
df$Mo_symb<-rowMeans(df[,c("Mo_Beauty", "Mo_Mystical", "Mo_Importance", "Mo_Abundance", "Mo_Exciting", "Mo_Inviting", "Mo_Relaxing", "Mo_Calm")])

# Mean scores for environments

mean(df$DS_symb, na.rm=TRUE)
mean(df$Ant_symb, na.rm=TRUE)
mean(df$RT_symb, na.rm=TRUE)
mean(df$Mo_symb, na.rm=TRUE)

########################################
## EFFECT OF DEMOGRAPHICS ON THE RESPOSES 
########################################

# Calculate mean values for the responses over the four environments

meandf<-df
meandf$mcare<-rowMeans(df[,c("Ant_care","DS_care", "RT_care", "Mo_care")])
meandf$msymb<-rowMeans(df[,c("Ant_symb","DS_symb", "RT_symb", "Mo_symb")])
meandf$merisk<-rowMeans(df[,c("EnvRisk_DS","EnvRisk_Ant", "EnvRisk_RT", "EnvRisk_Mo")])
meandf$msrisk<-rowMeans(df[,c("SocRisk_DS","SocRisk_Ant", "SocRisk_RT", "SocRisk_Mo")])
meandf$mknowe<-rowMeans(df[,c("DS_knowENV","Ant_knowENV", "RT_knowENV", "Mo_knowENV")])
meandf$mknowh<-rowMeans(df[,c("DS_knowHA","Ant_knowHA", "RT_knowHA", "Mo_knowHA")])

# Test differences between means of multiple groups with a non-parametric test
            kruskal.test(msymb~Age,data=meandf)
            kruskal.test(msymb~Gender,data=meandf) 
            kruskal.test(msymb~Education,data=meandf)

            kruskal.test(mcare~Age,data=meandf)
            kruskal.test(mcare~Gender,data=meandf) 
            kruskal.test(mcare~Education,data=meandf)

            kruskal.test(merisk~Age,data=meandf)
            kruskal.test(merisk~Gender,data=meandf) 
            kruskal.test(merisk~Education,data=meandf)

            kruskal.test(mknowe~Age,data=meandf)
            kruskal.test(mknowe~Gender,data=meandf) 
            kruskal.test(mknowe~Education,data=meandf)

            kruskal.test(mknowh~Age,data=meandf)
            kruskal.test(mknowh~Gender,data=meandf) 
            kruskal.test(mknowh~Education,data=meandf)

#########################################
#### DIFFERENCES BETWEEN ENVIRONMENTS ###
#########################################

# Mann-Whitney U test for checking statistical significance between distributions

library(plyr)

data2<-na.omit(data1[,c(3:59)])
data1<-data2[ ,!(colnames(data2) %in% c("LivingArea","Education","Country", "Gender","Age"))]

combos <- combn(ncol(data2),2)

MWcombo<-adply(combos, 2, function(x) {
  test <- wilcox.test(data2[, x[1]], data2[, x[2]])
  
  out <- data.frame("Row" = colnames(data2)[x[1]]
                    , "Column" = colnames(data2[x[2]])
                    ,  "df"= test$parameter
                    ,  "p.value" = round(test$p.value, 3)
  )
  
  return(out)
}) 


# test for siginificant difference between symbolic scores

DSsymb<-as.data.frame(df$DS_symb)
  DSsymb$Group <- paste0("DS ", DSsymb$Group) # set labels to groups
  colnames(DSsymb)<-c("Value", "Group")

Ansymb<-as.data.frame(df$Ant_symb)
  Ansymb$Group <- paste0("Ant ", Ansymb$Group) # set labels to groups
  colnames(Ansymb)<-c("Value", "Group")

rtsymb<-as.data.frame(df$RT_symb)
  rtsymb$Group <- paste0("RT ", rtsymb$Group) # set labels to groups
  colnames(rtsymb)<-c("Value", "Group")

mosymb<-as.data.frame(df$Mo_symb)
  mosymb$Group <- paste0("Mo ", mosymb$Group)
  colnames(mosymb)<-c("Value", "Group")

symb<-rbind(DSsymb,Ansymb,rtsymb,mosymb)
symb<-rbind(DSsymb,Ansymb)
            
fligner.test(Value ~ Group, data=symb) # test for equal variances, H0: variances equal
leveneTest(Value ~ Group, data = symb)

# doenst pass the test of homogeneity of vaiances, will use non-prametric test
  kruskal.test(Value ~ Group, data = symb)

  antest<-aov(Value ~ Group, data = symb)

#Visualize means
boxplot(Value ~ Group, data = symb, ylab="Symbolic value score", xlab="Environment")


###############################################
### CHECK INTERNAL CONSISTENCY OF RESPONSES ####
###############################################

# Cronbachs alpha
  psych::alpha(df) # check std.alpha, which is "the standardised alpha based upon the correlations"

###############################
##### PLOTTING THE RESULTS ####
###############################

## Deep sea values 5:12

            ds1<-ggplot(df, aes(x=df[,5]))+ xlab("Ugly-Beautiful")+geom_histogram(fill="navy")+theme_bw()
            ds2<-ggplot(df, aes(x=df[,6]))+ xlab("Ordinary-Mystical")+geom_histogram(fill="navy")+theme_bw()
            ds3<-ggplot(df, aes(x=df[,7]))+ xlab("Insignificant-Important ")+geom_histogram(fill="navy")+theme_bw()
            ds4<-ggplot(df, aes(x=df[,8]))+ xlab(" Empty - Abundant ")+geom_histogram(fill="navy")+theme_bw()
            ds5<-ggplot(df, aes(x=df[,9]))+ xlab(" Boring -Exciting")+geom_histogram(fill="navy")+theme_bw()
            ds6<-ggplot(df, aes(x=df[,10]))+ xlab(" Repelling-Inviting")+geom_histogram(fill="navy")+theme_bw()
            ds7<-ggplot(df, aes(x=df[,11]))+ xlab("Scary-Relaxing")+geom_histogram(fill="navy")+theme_bw()
            ds8<-ggplot(df, aes(x=df[,12]))+ xlab("Stressful-Calm")+geom_histogram(fill="navy")+theme_bw()

png("DeepseaValues_prelim.png")
DS<-grid.arrange(ds1,ds2,
ds3,
ds4,
ds5,
ds6,
ds7,
ds8, nrow=2, top="Deep sea")

dev.off()


## Antarctica values  16:23

            an1<-ggplot(df, aes(x=df[,16]))+ xlab("Ugly-Beautiful")+geom_histogram(fill="lightblue")+theme_bw()
            an2<-ggplot(df, aes(x=df[,17]))+ xlab("Ordinary-Mystical")+geom_histogram(fill="lightblue")+theme_bw()
            an3<-ggplot(df, aes(x=df[,18]))+ xlab("Insignificant-Important ")+geom_histogram(fill="lightblue")+theme_bw()
            an4<-ggplot(df, aes(x=df[,19]))+ xlab(" Empty - Abundant ")+geom_histogram(fill="lightblue")+theme_bw()
            an5<-ggplot(df, aes(x=df[,20]))+ xlab(" Boring -Exciting")+geom_histogram(fill="lightblue")+theme_bw()
            an6<-ggplot(df, aes(x=df[,21]))+ xlab(" Repelling-Inviting")+geom_histogram(fill="lightblue")+theme_bw()
            an7<-ggplot(df, aes(x=df[,22]))+ xlab("Scary-Relaxing")+geom_histogram(fill="lightblue")+theme_bw()
            an8<-ggplot(df, aes(x=df[,23]))+ xlab("Stressful-Calm")+geom_histogram(fill="lightblue")+theme_bw()

Ant<-grid.arrange(an1,
an2,
an3,
an4,
an5,
an6,
an7,
an8,nrow=2, top="Antarctica")


## Remote terrestrial values 27:34

            RT1<-ggplot(df, aes(x=df[,27]))+ xlab("Ugly-Beautiful")+geom_histogram(fill="forestgreen")+theme_bw()
            RT2<-ggplot(df, aes(x=df[,28]))+ xlab("Ordinary-Mystical")+geom_histogram(fill="forestgreen")+theme_bw()
            RT3<-ggplot(df, aes(x=df[,29]))+ xlab("Insignificant-Important ")+geom_histogram(fill="forestgreen")+theme_bw()
            RT4<-ggplot(df, aes(x=df[,30]))+ xlab(" Empty - Abundant ")+geom_histogram(fill="forestgreen")+theme_bw()
            RT5<-ggplot(df, aes(x=df[,31]))+ xlab(" Boring -Exciting")+geom_histogram(fill="forestgreen")+theme_bw()
            RT6<-ggplot(df, aes(x=df[,32]))+ xlab(" Repelling-Inviting")+geom_histogram(fill="forestgreen")+theme_bw()
            RT7<-ggplot(df, aes(x=df[,33]))+ xlab("Scary-Relaxing")+geom_histogram(fill="forestgreen")+theme_bw()
            RT8<-ggplot(df, aes(x=df[,34]))+ xlab("Stressful-Calm")+geom_histogram(fill="forestgreen")+theme_bw()

RT<-grid.arrange(RT1,
RT2,
RT3,
RT4,
RT5,
RT6,
RT7,
RT8, nrow=2, top="Remote terrestrial environments")


# Moon values

            Mo1<-ggplot(df, aes(x=df[,38]))+ xlab("Ugly-Beautiful")+geom_histogram(fill="darkorange")+theme_bw()
            Mo2<-ggplot(df, aes(x=df[,39]))+ xlab("Ordinary-Mystical")+geom_histogram(fill="darkorange")+theme_bw()
            Mo3<-ggplot(df, aes(x=df[,40]))+ xlab("Insignificant-Important ")+geom_histogram(fill="darkorange")+theme_bw()
            Mo4<-ggplot(df, aes(x=df[,41]))+ xlab(" Empty - Abundant ")+geom_histogram(fill="darkorange")+theme_bw()
            Mo5<-ggplot(df, aes(x=df[,42]))+ xlab(" Boring -Exciting")+geom_histogram(fill="darkorange")+theme_bw()
            Mo6<-ggplot(df, aes(x=df[,43]))+ xlab(" Repelling-Inviting")+geom_histogram(fill="darkorange")+theme_bw()
            Mo7<-ggplot(df, aes(x=df[,44]))+ xlab("Scary-Relaxing")+geom_histogram(fill="darkorange")+theme_bw()
            Mo8<-ggplot(df, aes(x=df[,45]))+ xlab("Stressful-Calm")+geom_histogram(fill="darkorange")+theme_bw()

Mo<-grid.arrange(Mo1,
Mo2,
Mo3,
Mo4,
Mo5,
Mo6,
Mo7,
Mo8, nrow=2, top="The Moon")

### KNOWLEDGE OF DIFFERENT ENVIRONMENTS

# Knowledge of environments
envknodf<- melt(data[,c("DS_knowENV","Ant_knowENV", "RT_knowENV", "Mo_knowENV")])

ggplot(envknodf, aes(value, fill=variable))+ 
  geom_bar(position = "dodge")+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.spacing.y=unit(0.2,"cm"))+
  theme(legend.key.size = unit(1, 'lines'))+
  ggtitle("Knowledge of the environmental conditions of different environments")+
  scale_fill_manual(values = c("navy",
                               "lightblue",
                               "forestgreen",
                               "darkorange"),labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  theme(legend.text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),legend.position=c(0.1, 0.9), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## SYMBOLIC VALUE SCORES

symbdf<- melt(df[,c("Mo_symb", "RT_symb","Ant_symb", "DS_symb")])

ggplot(symbdf, aes(value, fill=variable), alpha=0.8, lty="blank")+ 
  geom_bar(position = "dodge")+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.spacing.y=unit(0.2,"cm"))+
  theme(legend.key.size = unit(1, 'lines'))+
  ggtitle("Mean score of symbolic values for different environments")+
  scale_fill_manual(values = c("navy",
                               "lightblue",
                               "forestgreen",
                               "darkorange"),labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),legend.position=c(0.1, 0.9), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Distribution plot by groups

library(ggridges)# control overlap w scale, transparency w alpha

symbv<-
  ggplot(symbdf, aes(x = value, y = variable))+
  stat_density_ridges(aes(fill = variable), size=0.05, quantile_lines = FALSE,alpha=0.8, lty="blank") +
  xlab("Symbolic value score")+
  theme_ridges()+
  scale_fill_manual(values = c("darkorange" , "forestgreen", "lightblue","#00AFBB"))+
  theme(legend.position = "none")+
  scale_y_discrete(labels = c("Moon","Remote Terrestrial","Antarctica","Deep Sea"))+
  ylab("")+
theme(plot.title = element_text(hjust = 0.5),legend.position="none", panel.background = element_blank())


### ATTITUDES TO EXTRACTION

# Care for different environments

caredf<- melt(data[,c("DS_care","Ant_care", "RT_care", "Mo_care")])

ggplot(caredf, aes(value, fill=variable))+ 
  geom_bar(position = "dodge")+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.spacing.y=unit(0.2,"cm"))+
  theme(legend.key.size = unit(1, 'lines'))+
  ggtitle("Care for environments")+
  scale_fill_manual(values = c("navy",
                               "lightblue",
                               "forestgreen",
                               "darkorange"),labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  theme(legend.text=element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),legend.position=c(0.1, 0.9), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Histograms as ridgeplot

envrisk<-envriskdf %>%
  ggplot( aes(y=variable, x=value,  fill=variable, lty="blank")) +
  geom_density_ridges(alpha=0.8, stat="binline", bins=15, scale=0.9,lty="blank") +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(1.5, "lines"),
    strip.text.x = element_text(size = 8),
    text=element_text(size=12,  family="Arial")
  ) +
  
  scale_fill_manual(values = c("navy",
                               "lightblue",
                               "forestgreen",
                               "darkorange"),labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  scale_y_discrete(labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  xlab("Perception of environmental risk of mining") +
  theme(text=element_text(family="Lato",size=12))+
  ylab("")



# Cumulative plot

ggplot(socriskdf, aes(value, fill=variable)) + 
  stat_ecdf(aes(colour=variable),size = 0.8)+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.spacing.y=unit(0.2,"cm"))+
  ggtitle("Perception of societal risk of mineral extraction")+
  scale_colour_manual(values = c("navy",
                                 "lightblue",
                                 "forestgreen",
                                 "darkorange"),labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),legend.position=c(0.11, 0.89), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Histograms as ridgeplot

socrisk<-socriskdf %>%
  ggplot( aes(y=variable, x=value,  fill=variable, lty="blank")) +
  geom_density_ridges(alpha=0.8, stat="binline", bins=9, scale=0.9,lty="blank") +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(1.5, "lines"),
    strip.text.x = element_text(size = 8),
    text=element_text(size=12,  family="Arial")
  ) +
  
  scale_fill_manual(values = c("navy",
                               "lightblue",
                               "forestgreen",
                               "darkorange"),labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  scale_y_discrete(labels = c("Deep sea","Antarctica","Remote terrestrial", "Moon"))+
  xlab("Perception of societal risk of mining") +
  theme(text=element_text(family="Lato",size=12))+
  ylab("")

######################
### DATA ANALYSIS ####
######################

#### Correlations

table(data$Demo2, data$DS_care)
prop.table(table(data$Age, data$DS_care))  # get the percentages/proportions

prop.table(table(df$NEPSc, df$DS_care))

prop.table(table(data$BioV, data$DS_care))


## CHECKING PROPORTIONS for contingency tables

prop.table(table(df$Country))

### CHI SQUARE TEST ####

## Loop chi square test for all combinations of variables

library(plyr)

data2<-data[,3:length(data)]

combos <- combn(ncol(data2),2)

chicombo<-adply(combos, 2, function(x) {
  test <- chisq.test(data2[, x[1]], data2[, x[2]])
  
  out <- data.frame("Row" = colnames(data2)[x[1]]
                    , "Column" = colnames(data2[x[2]])
                    , "Chi.Square" = round(test$statistic,3)
                    ,  "df"= test$parameter
                    ,  "p.value" = round(test$p.value, 3)
  )
  return(out)
})  

write.table(chicombo, "\\\\ad.helsinki.fi/home/l/lmkaikko/Documents/chicombo_finaldata.txt", sep="\t") 


# ## Loop SPEARMAN&kendall CORRELATION for all combinations of variables

data2<-df[,c(3:length(df))]
data1<-data2[ ,!(colnames(data2) %in% c("LivingArea","Education","Country", "Gender","Age"))]

combos <- combn(ncol(data1),2)

corcombo<-adply(combos, 2, function(x) {
  test <- cor.test(data1[, x[1]], data1[, x[2]],method = c("kendall"))
  
  out <- data.frame("Row" = colnames(data1)[x[1]]
                    , "Column" = colnames(data1[x[2]])
                    , "R2" = round(test$estimate,3)
                    ,  "p.value" = round(test$p.value, 3)
  )
  
  return(out)
}) 

write.table(corcombo, "\\\\ad.helsinki.fi/home/l/lmkaikko/Documents/PhD tutkimus/corcombo_finaldata_kendall.txt", sep="\t") 

## ## Loop PEARSON CORRELATION for all combinations of variables

data2<-df[,c(3:length(df))]
data1<-data2[ ,!(colnames(data2) %in% c("LivingArea","Education","Country", "Gender","Age"))]

combos <- combn(ncol(data1),2)

corcombo<-adply(combos, 2, function(x) {
  test <- cor.test(data1[, x[1]], data1[, x[2]])
  
  out <- data.frame("Row" = colnames(data1)[x[1]]
                    , "Column" = colnames(data1[x[2]])
                    , "R2" = round(test$estimate,3)
                    ,  "df"= test$parameter
                    ,  "p.value" = round(test$p.value, 3)
  )

  return(out)
}) 

write.table(corcombo, "\\\\ad.helsinki.fi/home/l/lmkaikko/Documents/PhD tutkimus/corcombo_finaldata_new.txt", sep="\t") 


#####################################
##### BAYESIAN NETWORK APPROACH #####
#####################################

library(bnlearn)
data<-na.omit(df)

## ANTARCTIC BN ###

## Create antarctic dataframe 

#Selected parameters from the theoretical framework
Antdf<-df[,c("Ant_care","EnvRisk_Ant","SocRisk_Ant","Ant_ExtLH","NEPSc","Ant_knowENV","Ant_symb","Ant_knowHA", "BioV", "AltV", "EgoV", "HedV")]
Antcat<-na.omit(Antdf)# Transform data to categorical

# Discretise symbolic value and NEP scores

for (i in 1:(length(Antcat$NEPSc))){
if (Antcat$NEPSc[i]<2.1){
  Antcat$NEPSc[i]="Low"
} else if  (Antcat$NEPSc[i]>=2.1 & Antcat$NEPSc[i]<=3.5 ){
  Antcat$NEPSc[i]="Medium"
  
}  else if  (Antcat$NEPSc[i]>3.5 ){
  Antcat$NEPSc[i]="High"
  
}
}

for (i in 1:(length(Antcat$Ant_symb))){
  if (Antcat$Ant_symb[i]< 0){
    Antcat$Ant_symb[i]="Negative"
  } else if  (Antcat$Ant_symb[i]>=0 & Antcat$Ant_symb[i]<=1.5 ){
    Antcat$Ant_symb[i]="Neutral"
    
  }  else if  (Antcat$Ant_symb[i]>1.5 ){
    Antcat$Ant_symb[i]="Positive"
    
  }
}

Antcat[] <- lapply(Antcat, as.factor) #factorize for BN

# learn BN
dagL=hc(Antcat, score="bic", start = random.graph(names(Antcat)))#random start

dagAnt2<-mmhc(Antcat)

fitted = bn.fit(dagAnt2, data =Antcat)#learn CPTs

plot(dagL)
plot(dagAnt2)

## Bootstrapping for consensus network

# learn a set of 2500 network structures (takes ca. 3 min)

bootA <- boot.strength(Antcat, R = 2500, algorithm = "hc",
                       algorithm.args = list(score="bde",iss=10))

# select arcs which present in at least 70% of the DAGS

bootA[(bootA$strength > 0.434) & (bootA$direction >= 0.5), ]

# Having computed the significance for all possible arcs, we can now build averaged network using averaged.network and the appropriate threshold.
averaged.network(bootA)#check threshold
plot(averaged.network(bootA, threshold= 0.434)) # get all connections
avg.bootA <- averaged.network(bootA, threshold = 0.434)
plot(avg.bootA)

# check BIC score
score(avg.bootA, data = Antcat, type = "bic")

arc.strength(avg.bootA, data = Antcat, criterion = "bic")

## Alternative visualization
library(bnviewer)
bn.to.igraph(dagL)


viewer(dagL,
       bayesianNetwork.width = "100%",
       bayesianNetwork.height = "80vh",
       bayesianNetwork.layout = "layout_with_sugiyama",
       bayesianNetwork.title="Discrete Bayesian Network - Antarctica",
       bayesianNetwork.subtitle = "Monitoring of emergency care patients",
       bayesianNetwork.footer = "Fig. 1 - Layout with Sugiyama"
)

## DS BN ###

## Create Deep sea dataframe 

#Selected parameters form the theoretical framework
DSdf<-df[,c("DS_care","EnvRisk_DS","SocRisk_DS","DS_ExtLH","NEPSc","DS_knowENV","DS_symb","DS_knowHA", "BioV", "AltV", "EgoV", "HedV")]
DScat<-na.omit(DSdf)# Transform data to categorical

# modify symbolic value and NEP scores

for (i in 1:(length(DScat$NEPSc))){
  if (DScat$NEPSc[i]<2.1){
    DScat$NEPSc[i]="Low"
  } else if  (DScat$NEPSc[i]>=2.1 & DScat$NEPSc[i]<=3.5 ){
    DScat$NEPSc[i]="Medium"
    
  }  else if  (DScat$NEPSc[i]>3.5 ){
    DScat$NEPSc[i]="High"
    
  }
}


for (i in 1:(length(DScat$DS_symb))){
  if (DScat$DS_symb[i]< -1){
    DScat$DS_symb[i]="Negative"
  } else if  (DScat$DS_symb[i]>=-1 & DScat$DS_symb[i]<=1.5 ){
    DScat$DS_symb[i]="Neutral"
    
  }  else if  (DScat$DS_symb[i]>1.5 ){
    DScat$DS_symb[i]="Positive"
    
  }
}

DScat[] <- lapply(DScat, as.factor) #factorize for BN

# learn BN
dagDS = hc(DScat) # learnt structure
dagDS=hc(DScat, score="bic", start = random.graph(names(DScat))) # start at random

dagDS2<-mmhc(DScat)

fitted = bn.fit(dagDS, data =DScat)#learn CPTs

graphviz.plot(dagDS)
plot(dagDS)

score(dagDS, data = DScat, type = "bic")

arc.strength(dagDS, data = DScat, criterion = "x2")

#### Bootstrapping for consensus network ####

alpha.star(dagDS, DScat, debug = FALSE) # check ideal BDe value

# learn a set of 2500 network structures.

bootD <- boot.strength(DScat, R = 2500, algorithm = "hc", algorithm.args = list(score = "bde", iss = 10))

# select arcs which present in at least 70% of the DAGS

bootD[(bootD$strength >= 0.478 ) & (bootD$direction >= 0.5), ]

# Having computed the significance for all possible arcs, we can now build averaged network using averaged.network and the 85% threshold.

averaged.network(bootD) # check optimal threshold
avg.bootD <- averaged.network(bootD, threshold = 0.47)
plot(avg.bootD)

# check BIC score
score(avg.bootD, data = DScat, type = "bic")

arc.strength(avg.bootD, data = DScat, criterion = "bic")


#############
## RT BN ###
#############

## Create RT dataframe 

#Selected parameters from the theoretical framework
RTdf<-df[,c("RT_care","EnvRisk_RT","SocRisk_RT","RT_ExtLH","NEPSc","RT_knowENV","RT_symb","RT_knowHA", "BioV", "AltV", "EgoV", "HedV")]
RTcat<-na.omit(RTdf)# Transform data to categorical

# modify symbolic value and NEP scores

for (i in 1:(length(RTcat$NEPSc))){
  if (RTcat$NEPSc[i]<2.1){
    RTcat$NEPSc[i]="Low"
  } else if  (RTcat$NEPSc[i]>=2.1 & RTcat$NEPSc[i]<=3.5 ){
    RTcat$NEPSc[i]="Medium"
    
  }  else if  (RTcat$NEPSc[i]>3.5 ){
    RTcat$NEPSc[i]="High"
    
  }
}


for (i in 1:(length(RTcat$RT_symb))){
  if (RTcat$RT_symb[i]< -1){
    RTcat$RT_symb[i]="Negative"
  } else if  (RTcat$RT_symb[i]>=-1 & RTcat$RT_symb[i]<=1.5 ){
    RTcat$RT_symb[i]="Neutral"
    
  }  else if  (RTcat$RT_symb[i]>1.5 ){
    RTcat$RT_symb[i]="Positive"
    
  }
}

RTcat[] <- lapply(RTcat, as.factor) #factorize for BN

# learn BN
dagRT = hc(RTcat, score="bic", start = random.graph(names(RTcat))) # learnt structure w hill climbing with random start position( obtains 3 different results with very similar BIC scores)
# Evaluate
score(dagRT, data = RTcat, type = "bic")

choose.direction(dagRT, arc = c("RT_care", "RT_symb"),RTcat, debug = TRUE)

fitted = bn.fit(dagRT, data =RTcat)#learn CPTs

arc.strength(dagRT, data = RTcat, criterion = "x2")

choose.direction(dagRT, arc = c("RT_care", "SocRisk_RT"),RTcat, debug = TRUE)
choose.direction(dagRT, arc = c("RT_symb", "RT_knowENV"),RTcat, debug = TRUE)

### Bootstrapping for consensus network
# check ideal BDe value

alpha.star(dagRT, RTcat, debug = FALSE)

# learn a set of 2500 network structures.

set.seed(123)

bootR <- boot.strength(RTcat, R = 2500, algorithm = "hc",algorithm.args = list(score = "bde", iss = 10))

# select arcs which present in at least 85% of the DAGS

bootR[(bootR$strength > 0.557) & (bootR$direction >= 0.5), ]

# Having computed the significance for all possible arcs, we can now build averaged network using averaged.network and the 85% threshold.

averaged.network(bootR) # check optimal threshold = 0.558 
avg.bootR <- averaged.network(bootR, threshold = 0.557) 
plot(avg.bootR)


# check BIC score
score(avg.bootR, data = RTcat, type = "bic")

arc.strength(avg.bootR, data = RTcat, criterion = "bic")


## Moon BN ###

## Create Moon dataframe 
#Selected parameters form the theoretical framework
Modf<-df[,c("Mo_care","EnvRisk_Mo","SocRisk_Mo","Mo_ExtLH","NEPSc","Mo_knowENV","Mo_symb","Mo_knowHA", "BioV", "AltV", "EgoV", "HedV")]
Mocat<-na.omit(Modf)# Transform data to categorical

# modify symbolic value and NEP scores

for (i in 1:(length(Mocat$NEPSc))){
  if (Mocat$NEPSc[i]<2.1){
    Mocat$NEPSc[i]="Low"
  } else if  (Mocat$NEPSc[i]>=2.1 & Mocat$NEPSc[i]<=3.5 ){
    Mocat$NEPSc[i]="Medium"
    
  }  else if  (Mocat$NEPSc[i]>3.5 ){
    Mocat$NEPSc[i]="High"
    
  }
}


for (i in 1:(length(Mocat$Mo_symb))){
  if (Mocat$Mo_symb[i]< -1){
    Mocat$Mo_symb[i]="Negative"
  } else if  (Mocat$Mo_symb[i]>=-1 & Mocat$Mo_symb[i]<=1.5 ){
    Mocat$Mo_symb[i]="Neutral"
    
  }  else if  (Mocat$Mo_symb[i]>1.5 ){
    Mocat$Mo_symb[i]="Positive"
    
  }
}


Mocat[] <- lapply(Mocat, as.factor) #factorize for BN


# learn BN
dagMo = hc(Mocat) # learnt structure
dagMo= hc(Mocat, score="bic", start = random.graph(names(Mocat)))

dagMo2<-mmhc(Mocat)

fitted = bn.fit(dagMo, data =Mocat)#learn CPTs
plot(dagMo)

score(dagMo, data = Mocat, type = "bic")
arc.strength(dagMo, data = Mocat, criterion = "x2")

### Bootstrapping for consensus network
# check ideal BDe value

alpha.star(dagMo,Mocat, debug = FALSE)

# learn a set of 1000 network structures.
set.seed(123)

bootM <- boot.strength(Mocat, R = 2500, algorithm = "hc",algorithm.args = list(score = "bde", iss = 10))

# select arcs which present in at least 68% % of the DAGS; this threshold chosen because the best BIC score was obtained for this (as opposed to 70% and 63%)

bootM[(bootM$strength > 0.68) & (bootM$direction >= 0.5), ]

# Having computed the significance for all possible arcs, we can now build averaged network using averaged.network and the 85% threshold.
averaged.network(bootM)
avg.bootM <- averaged.network(bootM, threshold = 0.68) # CHECK inclusion threshold= 0.63
plot(avg.bootM)
arc.strength(avg.boot, data = RTcat, criterion = "x2")# arc strength

# check BIC score
score(avg.bootM, data = Mocat, type = "bic")

arc.strength(avg.bootM, data = Mocat, criterion = "bic")

### FURTHER BN STUFF TO PLAY WITH & ANALYTICS ####

# Defined structure

dagAnt<-empty.graph(nodes = c("Ant_care","EnvRisk_Ant","SocRisk_Ant","Ant_ExtLH","NEPSc","Ant_knowENV","Ant_symb","Ant_knowHA"))

arc.set1 = matrix(c("NEPSc", "EnvRisk_Ant",
                   "NEPSc", "Ant_care",
                  "Ant_ExtLH","EnvRisk_Ant",
                   "Ant_symb","EnvRisk_Ant",
                   "Ant_symb","Ant_care",
                  "Ant_knowENV","Ant_care",
                  "Ant_knowHA", "Ant_ExtLH",
                  "Ant_knowENV","Ant_symb",
                  "Ant_care","EnvRisk_Ant",
                  "EnvRisk_Ant", "SocRisk_Ant"
                   ),
                 byrow = TRUE, ncol = 2,
                 dimnames = list(NULL, c("from", "to")))

arcs(dagAnt) = arc.set1

## Refined strcuture based on tests

arc.set2 = matrix(c("NEPSc", "Ant_care",
                    "Ant_ExtLH","EnvRisk_Ant",
                    "Ant_symb","EnvRisk_Ant",
                    "Ant_symb","Ant_care",
                    "Ant_knowENV","Ant_care",
                    "Ant_knowHA", "Ant_ExtLH",
                    "Ant_knowENV","Ant_symb",
                    "Ant_care","EnvRisk_Ant",
                    "EnvRisk_Ant", "SocRisk_Ant"
),
byrow = TRUE, ncol = 2,
dimnames = list(NULL, c("from", "to")))
dagAnt2<-empty.graph(nodes = c("Ant_care","EnvRisk_Ant","SocRisk_Ant","Ant_ExtLH","NEPSc","Ant_knowENV","Ant_symb","Ant_knowHA"))
arcs(dagAnt2) = arc.set2

g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dagAnt))
graph::nodeRenderInfo(g) <- list(fontsize=40)
Rgraphviz::renderGraph(g)

graphviz.plot(dagAnt)

fitted = bn.fit(dagAnt, data =Antcat,method = "mle") # learning CPTs from data


