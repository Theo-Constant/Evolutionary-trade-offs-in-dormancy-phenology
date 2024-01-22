############################################################
#T.Constant 
#workspace prep
############################################################

#set directory
setwd(" ")

###load packages

library(caper)
library(apTreeshape)
library(phytools)
library(ape)
library(caper)
library(sjPlot)
library(devtools)
library(nlme)
library(geiger)
library(AICcmodavg)
library(sciplot)
library(MuMIn)
library(car)
library(ggplot2)
library(Rtools)

#Phylogenetic data


Tree_m_1<-read.nexus("model_1.nex")
Tree_m_2<-read.nexus("model_2.nex")
Tree_m_3<-read.nexus("model_3.nex")

#Biological data 
#all quantitative variables were standardized (using z-scores) in multi-factor models

model_1<-read.csv2("model_1.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
#standardization
model_1$zbody_mass_male=scale(model_1$body_mass_male,center=TRUE,scale = TRUE)
model_1$zbody_mass_change_during_mating=scale(model_1$body_mass_change_during_mating,center=TRUE,scale=TRUE)
model_1$zmin_temper=scale(model_1$min_temper,center=TRUE,scale=TRUE)
model_1$zlate_mating=scale(model_1$late_mating,center=TRUE,scale=TRUE)

model_2<-read.csv2("model_2.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
model_2$zprecipitation=scale(model_2$precipitation,center=TRUE,scale = TRUE)
model_2$zbody_mass_end_mating=scale(model_2$body_mass_end_mating,center=TRUE,scale=TRUE)
model_2$zlog_spe_repro_effort=scale(model_2$log_spe_repro_effort,center=TRUE,scale=TRUE)
model_2$zmaternal_effort=scale(model_2$maternal_effort,center=TRUE,scale=TRUE)

model_3<-read.csv2("model_3.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
#standardization
model_3$zprecipitation=scale(model_3$precipitation,center=TRUE,scale = TRUE)
model_3$zbody_mass_end_mating=scale(model_3$body_mass_end_mating,center=TRUE,scale=TRUE)
model_3$zlog_spe_repro_effort=scale(model_3$log_spe_repro_effort,center=TRUE,scale=TRUE)
model_3$zdimorphisme_imm=scale(model_3$dimorphisme_imm,center=TRUE,scale=TRUE)



#Build consensus tree from multiphylo object

Tree_m_1_cons<-consensus(Tree_m_1,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_1_cons

Tree_m_2_cons<-consensus(Tree_m_2,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_2_cons

Tree_m_3_cons<-consensus(Tree_m_3,p=1,check.labels = TRUE, rooted = TRUE)
Tree_m_3_cons

#Figure for paper

par(mfrow = c(1, 1))

plot(Tree_m_1_cons,cex=1,no.margin = TRUE,label.offset = 1)
plot(Tree_m_2_cons,cex=1,no.margin = TRUE,label.offset = 1)
plot(Tree_m_3_cons,cex=1,no.margin = TRUE,label.offset = 1)


#Determine branch lengths of consensus trees 
#Are the tree ultrametric?  (the answer must be TRUE)


Tree_m_1_comp <- compute.brlen(Tree_m_1_cons, method="Grafen") 
is.ultrametric(Tree_m_1_comp)

Tree_m_2_comp <- compute.brlen(Tree_m_2_cons, method="Grafen") 
is.ultrametric(Tree_m_2_comp)

Tree_m_3_comp <- compute.brlen(Tree_m_3_cons, method="Grafen") 
is.ultrametric(Tree_m_3_comp)
#Are the tree dichotomus? (the answer must be TRUE)

is.binary.multiPhylo(Tree_m_1)
is.binary.tree(Tree_m_1_comp)

is.binary.multiPhylo(Tree_m_2)
is.binary.tree(Tree_m_2_comp)

is.binary.multiPhylo(Tree_m_3)
is.binary.tree(Tree_m_3_comp)
#Combine both datafiles (phylo and Biological data)


Tree_m_1_comp$node.label<-NULL
comb_model_1<-comparative.data(phy=Tree_m_1_comp,data=model_1,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_1)
sort(unique(Tree_m_1_comp$tip.label))
sort(unique(model_1$species))
sort(unique(comb_model_1$phy$tip.label))
comb_model_1$dropped$tips
sort(unique(comb_model_1$phy$tip.label))==sort(unique(model_1$species))

Tree_m_2_comp$node.label<-NULL
comb_model_2<-comparative.data(phy=Tree_m_2_comp,data=model_2,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_2)
sort(unique(Tree_m_2_comp$tip.label))
sort(unique(model_2$species))
sort(unique(comb_model_2$phy$tip.label))
comb_model_2$dropped$tips
sort(unique(comb_model_2$phy$tip.label))==sort(unique(model_2$species))

Tree_m_3_comp$node.label<-NULL
comb_model_3<-comparative.data(phy=Tree_m_3_comp,data=model_3,names.col=species,vcv=TRUE,na.omit=FALSE,warn.dropped=TRUE)
str(comb_model_3)
sort(unique(Tree_m_3_comp$tip.label))
sort(unique(model_3$species))
sort(unique(comb_model_3$phy$tip.label))
comb_model_3$dropped$tips
sort(unique(comb_model_3$phy$tip.label))==sort(unique(model_3$species))

#PGLS models

#pgls model with lambda evaluation by ML (1=high phylogenetic covariance; 0 = no phylogenetic covariance)

#model 1

fit1_1<-pgls(protandry~late_mating+foodstoring+body_mass_male+dimorphism_body_mass_emergence+body_mass_change_during_mating*min_temper+
              body_mass_change_during_mating*precipitation
            ,data=comb_model_1,lambda ="ML")
summary(model.avg(dredge(fit1_1, rank="AICc",m.lim = c(NA, 4)),delta<3))

fit1_2<-pgls(protandry~zbody_mass_male+zbody_mass_change_during_mating+zmin_temper+zlate_mating,data=comb_model_1,lambda ="ML")
summary(fit1_2)

#Normality and homoscedasticity are checked by graphical observation

#1) normality
#Ploting normal qqplot
#If the data are normally distributed, the points on a Q-Q graph 
#will lie on a straight diagonal line.

#2) homoscedasticity 
#Ploting residual value vs fitted value allows to check the homoscedasticity
#The residuals have to be randomly distributed around the 0 line, 
#approximately form a horizontal band around the 0 line 
#with no extreme value from the basic random pattern of residuals. 

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit1_2)


#equation for regression line of model 6 with mean value of late mating
fit1_3<-pgls(protandry~
              body_mass_male+body_mass_change_during_mating+min_temper+late_mating,data=comb_model_1,lambda ="ML")
summary(fit1_3)

eq1 = function(x){ (20.6398273  +(0.0037926*499.9)+(0.712107*-8.31)+(-1.9893747*1.55)+(-0.4381697*x))}

ggplot(model_1, aes(x=body_mass_change_during_mating, y=protandry,colour=min_temper))+
  geom_point(size=4)+
  theme_classic()+
  
  
  labs(x="Body mass change during mating (% body mass)",
       y="Protandry (day)")+
  scale_y_continuous(limits = c(0,55), expand = c(0, 0)) +
  scale_color_gradient2("Min temperature (°C)",midpoint=-8.31,  low="blue4", mid="snow3",
                        high="red4",  space = "Lab")+
  geom_function(fun=eq1, color="black",linewidth=1)



#Model 2


fit2_1<-pgls(sex_diff_imm~log_spe_repro_effort*min_temp+log_spe_repro_effort*precipitation+body_mass_imm+min_temp*body_mass_end_mating+body_mass_end_mating*precipitation+maternal_effort*precipitation+maternal_effort*min_temp+dimorphisme_imm,data=comb_model_2,lambda = "ML")
summary(model.avg(dredge(fit2_1, rank="AICc",m.lim = c(NA, 4)),delta<6))

fit2_2<-pgls(sex_diff_imm~zlog_spe_repro_effort+zbody_mass_end_mating+zprecipitation+zmaternal_effort,data=comb_model_2,lambda = "ML")
summary(fit2_2)


#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit2_2)

#the points on the Q-Q graph 
#don't follow a a straight diagonal line.


#We found an an outlier, Tachyglossus aculeatus, a species for which gestation and lactation for females was an extremely long period (i.e. 167 day). 
#we test a new model (model 3) without this specie

#Model 3
fit3_1<-pgls(sex_diff_imm~log_spe_repro_effort*min_temp+log_spe_repro_effort*precipitation+body_mass_imm+min_temp*body_mass_end_mating+body_mass_end_mating*precipitation+maternal_effort*precipitation+maternal_effort*min_temp+dimorphisme_imm,data=comb_model_3,lambda = "ML")
summary(model.avg(dredge(fit3_1, rank="AICc",m.lim = c(NA, 4)),delta<6))

fit3_2<-pgls(sex_diff_imm~zlog_spe_repro_effort+zbody_mass_end_mating+zprecipitation+zdimorphisme_imm,data=comb_model_3,lambda = "ML")
summary(fit3_2)

#Normality and homoscedasticity are checked by graphical observation

par(mfrow = c(2, 2))
par(mar=c(4,4,2,2))
plot(fit3_2)


#equation for regression line of model 6 with mean value of late mating

fit3_3<-pgls(sex_diff_imm~log_spe_repro_effort+body_mass_end_mating+precipitation+dimorphisme_imm,data=comb_model_3,lambda = "ML")
summary(fit3_3)
eq2 = function(x){ 180.103912+(-57.962616*2.25)+(-0.049655*583.5)+(-23.754942*1.2)+(-0.709451*x)}


ggplot(model_3, aes(x=body_mass_end_mating, y=sex_diff_imm,colour=precipitation))+
  geom_point(size=4)+
  theme_classic()+
  scale_colour_gradient("Precipitation",
    low = "#56B1F7",
    high = "#132B43")+
  ylab("Sex difference in immergence (day)")+
  xlab("Body mass variation through the end of mating (% body mass)")+
  geom_function(fun=eq2, color="black",linewidth=1)


eq3 = function(x){ 180.103912+(-0.709451*-5.22)+(-0.049655*583.5)+(-23.754942*1.2)+(-57.962616*x)}


ggplot(model_3, aes(x=log_spe_repro_effort, y=sex_diff_imm,colour=precipitation))+
  geom_point(size=4)+
  theme_classic()+
  scale_colour_gradient("Precipitation",
                        low = "#56B1F7",
                        high = "#132B43")+
  ylab("Sex difference in immergence (day)")+
  xlab("log(Female specific reproductive effort)")+
  geom_function(fun=eq3, color="black",linewidth=1)


