### Congruence Three Way Interaction ANOVA ###

#clear environment
rm(list=ls())


#settings
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(afex)
library(rstatix)
library(emmeans)
library(moments)

#set working directory
setwd("R:/Typeface_Stroop_Task/Data_Analysis/TaskData/cleandata/collateddata")

#set conditions
condfile1<- "nameColour_descriptivestats"
condfile2<- "nameWord_descriptivestats"

#set file names
filename1 <- paste(condfile1,"_cong.csv",sep="")
filename2 <- paste(condfile1,"_incong.csv",sep="")
filename3 <- paste(condfile2,"_cong.csv",sep="")
filename4 <- paste(condfile2,"_incong.csv",sep="")

#load data
cong_nc<- read.csv(filename1)
incong_nc <- read.csv(filename2)
cong_nw<- read.csv(filename3)
incong_nw <- read.csv(filename4)

### Organise data structures
cong_nc<-as.matrix(cong_nc)
rnames<- cong_nc[,1]
cong_nc<-cong_nc[,-1]
rownames(cong_nc)<- rnames
cong_nc<-t(cong_nc)
cong_nc<-as.data.frame(cong_nc)
cong_nc_sessions<-cong_nc$session
cong_nc<-apply(cong_nc,2,as.numeric)
cong_nc<-as.data.frame(cong_nc)
cong_nc$session<-cong_nc_sessions


cong_nw<-as.matrix(cong_nw)
rnames<- cong_nw[,1]
cong_nw<-cong_nw[,-1]
rownames(cong_nw)<- rnames
cong_nw<-t(cong_nw)
cong_nw<-as.data.frame(cong_nw)
cong_nw_sessions<-cong_nw$session
cong_nw<-apply(cong_nw,2,as.numeric)
cong_nw<-as.data.frame(cong_nw)
cong_nw$session<-cong_nw_sessions


incong_nc<-as.matrix(incong_nc)
rnames<- incong_nc[,1]
incong_nc<-incong_nc[,-1]
rownames(incong_nc)<- rnames
incong_nc<-t(incong_nc)
incong_nc<-as.data.frame(incong_nc)
incong_nc_sessions<-incong_nc$session
incong_nc<-apply(incong_nc,2,as.numeric)
incong_nc<-as.data.frame(incong_nc)
incong_nc$session<-incong_nc_sessions


incong_nw<-as.matrix(incong_nw)
rnames<- incong_nw[,1]
incong_nw<-incong_nw[,-1]
rownames(incong_nw)<- rnames
incong_nw<-t(incong_nw)
incong_nw<-as.data.frame(incong_nw)
incong_nw_sessions<-incong_nw$session
incong_nw<-apply(incong_nw,2,as.numeric)
incong_nw<-as.data.frame(incong_nw)
incong_nw$session<-incong_nw_sessions

#add sessions back in
cong_nc$session<-cong_nc_sessions
cong_nw$session<-cong_nw_sessions
incong_nc$session<-incong_nc_sessions
incong_nw$session<-incong_nw_sessions


#make ID variable and condition
cong_nc$ID<-rownames(cong_nc)
incong_nc$ID<-rownames(incong_nc)
cong_nw$ID<-rownames(cong_nw)
incong_nw$ID<-rownames(incong_nw)

cong_nc$congruency[1:length(cong_nc$mean)]<-"cong"
cong_nw$congruency[1:length(cong_nw$mean)]<-"cong"
incong_nc$congruency[1:length(incong_nc$mean)]<-"incong"
incong_nw$congruency[1:length(incong_nw$mean)]<-"incong"

cong_nc$condition[1:length(cong_nc$mean)]<-"nameColour"
cong_nw$condition[1:length(cong_nw$mean)]<-"nameWord"
incong_nc$condition[1:length(incong_nc$mean)]<-"nameColour"
incong_nw$condition[1:length(incong_nw$mean)]<-"nameWord"

###combine into single data frame
finaldata<-rbind(cong_nc,incong_nc,cong_nw,incong_nw)

### delete participants with large error rate
finaldata$`error_rate_(%)`<-as.numeric(finaldata$`error_rate_(%)`)
delfinaldata<-c()
for(e in 1:length(finaldata$`error_rate_(%)`)){
  if(finaldata$`error_rate_(%)`[e] > 50){
    print(finaldata$ID[e])
    delfinaldata[e]<- e
  }
}

delfinaldata<-na.omit(delfinaldata)
delfinaldata<-as.vector(delfinaldata)

finaldata<-finaldata[-c(delfinaldata),]


### translate data to long form so ANOVA can run
means<-c(finaldata$mean_L1,finaldata$mean_L2,finaldata$mean_L3,finaldata$mean_L4)
ERs<-c(finaldata$ER_L1,finaldata$ER_L2,finaldata$ER_L3,finaldata$ER_L4)
condition<-c(finaldata$condition,finaldata$condition,finaldata$condition,finaldata$condition)
congruency<-c(finaldata$congruency,finaldata$congruency,finaldata$congruency,finaldata$congruency)
ID<-c(finaldata$ID,finaldata$ID,finaldata$ID,finaldata$ID)

#set disfluency levels
ln1<-length(finaldata$mean)
ln2<-ln1+ln1
ln3<-ln1+ln2
ln4<-ln1+ln3
level<-c()
level[1:ln1]<-1
level[(ln1+1):ln2]<-2
level[(ln2+1):ln3]<-3
level[(ln3+1):ln4]<-4


analysis<-data.frame(ID,level,condition,congruency,means,ERs)
analysis$level<-as.factor(analysis$level)
analysis$ID<-as.factor(analysis$ID)

### Visualisation
#Means by level
ggboxplot(data=analysis,x="level", y="means")

#Means by condition
ggboxplot(data=analysis,x="condition", y="means")

#histogram
hist(analysis$means)

#check and remove outliers by level
l1<-filter(.data = analysis, analysis$level == 1)
l2<-filter(.data = analysis, analysis$level == 2)
l3<-filter(.data = analysis, analysis$level == 3)
l4<-filter(.data = analysis, analysis$level == 4)

meanzscore1<- (l1$means[1:length(l1$means)]-mean(l1$means))/sd(l1$means)
outliers1<-c()
for(i in 1:length(l1$means)){
  if(meanzscore1[i] > 3 | meanzscore1[i] < -3  ){
    outliers1[i]<-i
  }
}
outliers1<-na.omit(outliers1)

l1<-l1[-c(outliers1),]


meanzscore2<- (l2$means[1:length(l2$means)]-mean(l2$means))/sd(l2$means)
outliers2<-c()
for(i in 1:length(l2$means)){
  if(meanzscore2[i] > 3 | meanzscore2[i] < -3  ){
    outliers2[i]<-i
  }
}
outliers2<-na.omit(outliers2)

l2<-l2[-c(outliers2),]

meanzscore3<- (l3$means[1:length(l3$means)]-mean(l3$means))/sd(l3$means)
outliers3<-c()
for(i in 1:length(l3$means)){
  if(meanzscore3[i] > 3 | meanzscore3[i] < -3  ){
    outliers3[i]<-i
  }
}
outliers3<-na.omit(outliers3)

l3<-l3[-c(outliers3),]


meanzscore4<- (l4$means[1:length(l4$means)]-mean(l4$means))/sd(l4$means)
outliers4<-c()
for(i in 1:length(l4$means)){
  if(meanzscore4[i] > 3 | meanzscore4[i] < -3  ){
    outliers4[i]<-i
  }
}
outliers4<-na.omit(outliers4)

l4<-l4[-c(outliers4),]



#recombine data
analysis_m<-rbind(l1,l2,l3,l4)

#identify participants that had their data removed
rem<-c()

rem$l1<-as.character(analysis$ID[outliers1])
rem$l2<-as.character(analysis$ID[outliers2])
rem$l3<-as.character(analysis$ID[outliers3])
rem$l4<-as.character(analysis$ID[outliers4])
rem<- c(rem$l1,rem$l2,rem$l3,rem$l4)

unique_r<-unique(rem)
unique_r

###ANOVA
mean_anova<-aov_ez(id = "ID", dv = "means", data = analysis_m, within = c("level","condition","congruency"))
mean_anova

#histogram of residuals
hist(mean_anova$lm$residuals)

#skewness and kurtosis of residuals
mean(skewness(mean_anova$lm$residuals))
mean(kurtosis(mean_anova$lm$residuals))

#### three-way interaction post-hocs
emmeans(mean_anova, "level", contr = "pairwise", adjust = "none")

#dislfuency level comparisons
contrasting <- emmeans(mean_anova, ~ level|congruency*condition)
contrast(contrasting, method =  "pairwise", adjust = "none")
eff_size(contrasting, sigma = mean(sigma(mean_anova$lm)), edf=df.residual(mean_anova$lm))


#congruence comparisons
contrasting <- emmeans(mean_anova, ~ congruency|condition*level)
contrast(contrasting, method =  "pairwise", adjust = "none")
eff_size(contrasting, sigma = mean(sigma(mean_anova$lm)), edf=df.residual(mean_anova$lm))


#condition comparisons
contrasting <- emmeans(mean_anova, ~ condition|congruency*level)
contrast(contrasting, method =  "pairwise", adjust = "none")
eff_size(contrasting, sigma = mean(sigma(mean_anova$lm)), edf=df.residual(mean_anova$lm))

#two-way post-hocs
#dislfuency level x condition comparisons
contrasting <- emmeans(mean_anova, ~ level|condition)
contrast(contrasting, method =  "pairwise", adjust = "none")
eff_size(contrasting, sigma = mean(sigma(mean_anova$lm)), edf=df.residual(mean_anova$lm))

#congruence x condition comparisons
contrasting <- emmeans(mean_anova, ~ congruency|condition)
contrast(contrasting, method =  "pairwise", adjust = "none")
eff_size(contrasting, sigma = mean(sigma(mean_anova$lm)), edf=df.residual(mean_anova$lm))


# CALCULATE MEANS AND SDS USED IN THE ANALYSIS BY EACH FACTOR
level<-unique(mean_anova$data$long$level)
congruence<-unique(mean_anova$data$long$congruency)
condition<-unique(mean_anova$data$long$condition)

#means and sds stratified by all 3 factors
countr<-0
means<-c()
sds<-c()
cols<-c()
for(c in condition){
  for(co in congruence){
    for(l in level){
      d <-filter(
        mean_anova$data$long,
        mean_anova$data$long$level==l & mean_anova$data$long$condition == c & mean_anova$data$long$congruency == co)
      d<-d$means
      m<-mean(d)
      s<-sd(d)
      print(paste(c, co, l, "mean =", m, sep = ' '))
      print(paste(c, co, l, "SD =", s, sep = ' '))
      countr<-countr+1
      means[countr]<-m
      sds[countr]<-s
      cols[countr]<-paste(c, co, l,sep='_')
    }
  }
}

strat_m_sd<-matrix(nrow = 2, ncol = length(means))
colnames(strat_m_sd)<-cols
rownames(strat_m_sd)<-c('means','sds')
strat_m_sd[1,]<-means
strat_m_sd[2,]<-sds

#means and sds for congruence
mc<-filter(
  mean_anova$data$long,
  mean_anova$data$long$congruency == 'cong'
  )$means
print(paste("congruent mean =",mean(mc)))
print(paste("congruent sd =",sd(mc)))

mi<-filter(
  mean_anova$data$long,
  mean_anova$data$long$congruency == 'incong'
)$means
print(paste("incongruent mean =",mean(mi)))
print(paste("incongruent sd =",sd(mi)))


#means and sds for condition x congruence interation
mcc<-filter(
  mean_anova$data$long,
  mean_anova$data$long$congruency == 'cong' & mean_anova$data$long$condition == 'nameColour'
)$means
print(paste("congruent NC mean =",mean(mcc)))
print(paste("congruent NC sd =",sd(mcc)))

mcw<-filter(
  mean_anova$data$long,
  mean_anova$data$long$congruency == 'cong' & mean_anova$data$long$condition == 'nameWord'
)$means
print(paste("congruent NW mean =",mean(mcw)))
print(paste("congruent NW sd =",sd(mcw)))

mic<-filter(
  mean_anova$data$long,
  mean_anova$data$long$congruency == 'incong' & mean_anova$data$long$condition == 'nameColour'
)$means
print(paste("incongruent NC mean =",mean(mic)))
print(paste("incongruent NC sd =",sd(mic)))

miw<-filter(
  mean_anova$data$long,
  mean_anova$data$long$congruency == 'incong' & mean_anova$data$long$condition == 'nameWord'
)$means
print(paste("incongruent NW mean =",mean(miw)))
print(paste("incongruent NW sd =",sd(miw)))


#boxplot data wrangling
analysis_m$interaction <- interaction(analysis_m$level, analysis_m$congruency, analysis_m$condition, sep=":")
analysis_m$interaction2 <- interaction(analysis_m$congruency,analysis_m$level, sep=":")

anno_df<-compare_means(means~interaction, group.by="interaction2", data = analysis_m, p.adjust.method = "holm", method = "t.test")%>%
  mutate(y_pos = c(1300,1500,1350,1550,1400,1600,1450,1650), p.adj = format.pval(p.adj, digits = 2))

anno_df$congruence<-c("cong","incong","cong","incong","cong","incong","cong","incong")
anno_df<-anno_df[order(anno_df$congruence),]
anno_df<-anno_df[-c(1,2,3,7),]
anno_df$p.signif<-c("","","","")
anno_df$y_pos<-anno_df$y_pos-200
positions<-c(1250)
for(i in 2:4){
  positions[i]<-positions[i-1]+50
}
anno_df$y_pos<-positions
rm(congruency)
analysis_m$congruence<-analysis_m$congruency
rm(analysis_m$congruence)

#boxplot
boxplot<-ggplot(aes(y = means, x = interaction, fill=congruence), data = analysis_m) + 
  geom_boxplot()
boxplot<-boxplot + 
  scale_fill_discrete(name = "Congruence", labels = c("Congruent","Incongruent")) +
  labs(title="Three Way Interaction", x = "Colour Identification                                                                                             Word Identification", y = "Reaction Time (ms)")+
  scale_x_discrete(labels= c("D1","D2","D3","D4","D1","D2","D3","D4","D1","D2","D3","D4","D1","D2","D3","D4"))

boxplot

    
#IDs for demographics calculations
#get IDs from ANOVA df
ID<-as.character(unique(mean_anova$data$long$ID))
sessions<-c()
for(i in 1:length(ID)){
  sessions[i]<-filter(cong_nc,cong_nc$ID == ID[i])$session
}

#load demographics table
fn <- "Stroop participants-demographic_data.csv"
demdata <- read.csv(fn)

#remove excluded participants
locs<-c()
for(s in seq_along(sessions)){
  locs[s]<-which(demdata$session_id == sessions[s])
}

demdata<-demdata[locs,]

#mean age
mean(demdata$age, na.rm=T)

#number of males
sum(demdata$Sex == 'Male')

