#----START----
rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/3eme année/isotopic turnover/DEB/AmP_data")

library(DHARMa)
library(ggplot2)

#importing functions
source("functions_lambda_DEB.r")

#----OPENING DATA----
##DEB database##
data <- read.csv2("DEB_paramV6.csv",sep=",")
data<-cbind(data[,1:2],sapply(data[,-c(1,2)],as.numeric))



##opening dataset 1(T&C and VZ)##
observations <- read.csv2("dataset_commun_DEB_add_temp.csv",fileEncoding = "latin1")
#adding t1_2_corr at 20°C for observations
TC_obs=NULL
for (i in 1:nrow(observations)){
  # i=1
  if (observations$model[i]=="no"){
    TC_obs=c(TC_obs,temp_correction(T_A=8000,T_ref=293.15,t=(observations$temp_..C.[i]+273.15)))
  }
  else {
    name=observations$Species[i]
    param=subset(data,data$Species==name)
    TC_obs=c(TC_obs,temp_correction(T_A=param$T_A,T_ref=param$T_ref,t=(observations$temp_..C.[i]+273.15)))
  }
  
}
observations$TC<- TC_obs
observations$t1_2_corr <- observations$t_1.2_.d.*observations$TC
##Splitting dataset between specoes in DEB and thos without
obs_noDEB=subset(observations,observations$model=="no")
obs=subset(observations,observations$model!="no")


###opening dataset 2##
review_new<- read.csv2("half_life_review.csv",fileEncoding = "latin1")[,1:11]
########TAKING OUT Perca_fluviatilis and Tursiops_truncatus for
#we don't have mass for now
review_new_bis <- subset(review_new,review_new$Species!="Perca_fluviatilis"& review_new$Species!="Tursiops_truncatus")
#Adding a TC col
TC=NULL
for (i in 1:nrow(review_new_bis)){
  # i=28
  name=review_new_bis$Species[i]
  DEB_param =subset(data,data$Species==name)
  TC=c(TC,temp_correction(T_A=DEB_param$T_A,T_ref=DEB_param$T_ref
                          ,t=(review_new_bis$Body.temperature...C.[i]+273.15)))
  
}
review_new_bis$TC <- TC
review_new_bis$mass<- review_new_bis$Mass..g.
review_new_bis$t1_2_deg <- review_new_bis$t1.2*review_new_bis$TC #corr at 20°C


#----RELASHIONSHIP LAMBDA /WW ALL SPECIES----


#for every species in the DEB base between f=0.5 and f=1
data_sub=subset(data,data$Model=="std"|data$Model=="abj"|data$Model=="asj"|data$Model=="stf"|data$Model=="stx")
every_species=unique(data_sub$Species)
# rel_every_species= relashionship(num_points=10,f_val=seq(0.5,1,by=0.1)
#                                  ,species_names=every_species,data=data)
# write.csv2(rel_every_species, file = "relationships_all_species.csv", row.names = FALSE)
rel_every_species=read.csv2("relationships_all_species.csv")


#####graph relathionship: FIGURE 5#####
#at f=0.8

rel_every_species_f1=subset(rel_every_species,rel_every_species$f==0.8)
tiff("relathionship_dataset.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(-10,8),ylim=c(-2,4), xlab="log10(weight)",ylab="log10(t1/2)")
for (i in 1:length(every_species)){
  # i=19
  name=every_species[i]
  subs_2=subset(rel_every_species_f1,rel_every_species_f1$species==name)
  points(log10(subs_2$Ww),log10(subs_2$t1_2),type="l",col="grey")
}
points(log10(obs_noDEB$Mass_.g.),log10(obs_noDEB$t1_2_corr),col="green",pch=16)
points(log10(obs$Mass_.g.),log10(obs$t1_2_corr),type="p",col="red",pch=16)
points(log10(review_new_bis$mass),log10(review_new_bis$t1_2_deg),type="p",col="blue",pch=16)

dev.off()

####finding the outlier ####
# high_mass=subset(rel_every_species,rel_every_species$L>0.5)
# 
# min(high_mass$t1_2)
# qui=subset(high_mass,high_mass$t1_2<0.022)
# 
# coupable=subset(rel_every_species,rel_every_species$species=="Nematoflustra_flagellata")
# 
# plot(xlim=c(-10,8),ylim=c(-2,4),log10(coupable$Ww),log10(coupable$t1_2))
#----FIGURE 6----

#legend
tiff("Legend figure 8.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c(" λ-DEB Predictions","Observations 1","Observations 2","C","N","S")
       ,cex=1,col=c("grey","red","blue","black","black","black")
       ,lty=c(1,NA,NA,NA,NA,NA),box.lty=0,pch=c(NA,16,16,16,15,17))
dev.off()

#####observations and predictions dataset 1  #####

rel_lam_ww= relashionship(num_points=10,f_val=seq(0.5,1,by=0.1)
                          ,species_names=unique(obs$Species),data=data)
species_names=unique(obs$Species)
f_val=seq(0.5,1,by=0.1)
for (i in 1:length(species_names)){ #for each species
  # i=1
  name=species_names[i]
  #isolating predictions and observations for the species
  subs_rel=subset(rel_lam_ww, rel_lam_ww$species==name)
  subs_obs=subset(obs,obs$Species==name)
  title=paste0(name,".tiff")
  tiff(title, , width =4 , height = 5, units = "in", res = 300)
  plot(NULL,main=name,xlim=c(-8,6),ylim=c(-1,3), xlab="log10(weight)",ylab="log10(t1/2)")
  
  for (j in 1:length(f_val)){ #for each f 
    # j=1
    # f=f_val[j]
    subs_rel_f= subset(subs_rel,subs_rel$f==f_val[j])
    points(log10(subs_rel_f$Ww),log10(subs_rel_f$t1_2),type="l",col="grey")
  }
  num=NULL
  for (j in 1:nrow(subs_obs)){
    if ( subs_obs$element[j]=="C"){num=c(num,16)}
    if ( subs_obs$element[j]=="N"){num=c(num,15)}
    if ( subs_obs$element[j]=="S"){num=c(num,17)}
  }
  points(log10(subs_obs$Mass_.g.),log10(subs_obs$t1_2_corr),type="p",col="red",pch=num)
  
  dev.off()
}


#####observations and predictions dataset 2  #####

rel_lam_ww= relashionship(num_points=10,f_val=seq(0.5,1,by=0.1)
                          ,species_names=unique(review_new_bis$Species),data=data)
species_names=unique(review_new_bis$Species)
f_val=seq(0.5,1,by=0.1)
for (i in 1:length(species_names)){ #for each species
  # i=1
  name=species_names[i]
  #isolating predictions and observations for the species
  subs_rel=subset(rel_lam_ww, rel_lam_ww$species==name)
  subs_obs=subset(review_new_bis,review_new_bis$Species==name)
  title=paste0(name,".tiff")
  tiff(title , width =4 , height = 5, units = "in", res = 300)
  plot(NULL,main=name,xlim=c(-8,6),ylim=c(-1,3), xlab="log10(weight)",ylab="log10(t1/2)")
  
  for (j in 1:length(f_val)){ #for each f 
    # j=1
    # f=f_val[j]
    subs_rel_f= subset(subs_rel,subs_rel$f==f_val[j])
    points(log10(subs_rel_f$Ww),log10(subs_rel_f$t1_2),type="l",col="grey")
  }
  num=NULL
  for (j in 1:nrow(subs_obs)){
    if ( subs_obs$element[j]=="C"){num=c(num,16)}
    if ( subs_obs$element[j]=="N"){num=c(num,15)}
    if ( subs_obs$element[j]=="S"){num=c(num,17)}
  }
  points(log10(subs_obs$mass),log10(subs_obs$t1_2_deg),type="p",col="blue",pch=num)
  
  dev.off()
}



#----FIGURE 7 (dataset 1)----

##lambda-DEB predictions vs obs####
f=0.8
#for each obs estimating structural length and if longer than Li, L is Li
obs$Mass<-obs$Mass_.g.
L_obs=size_est_corr(f=0.8,obs=obs,data=data)
obs$L_corr<- L_obs$L_corr

#for each obs with L predicting lambda and t1/2
obs$t1_2_pred<- pred_lam_half_life(f=f,obs=obs,data=data)


tiff("ATTENTION ESSAI.tiff",  width =4 , height = 5, units = "in", res = 300)

# tiff("obs_vs_pred_DEB_OBS_f0.8_.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,3),ylim=c(-1,3),xlab="observed Half-life",
     ylab="λ-DEB predicted Half-life")

for (i in 1:nrow(obs)){
  if (obs$tissue[i]=="whole body"& obs$type[i]=="endotherm" ){
    color="red4"
    if ( obs$element[i]=="C"){num=1}
    if ( obs$element[i]=="N"){num=0}
    if ( obs$element[i]=="S"){num=2}
  }
  if (obs$tissue[i]=="muscle"& obs$type[i]=="endotherm" ){
    color="red"
    if ( obs$element[i]=="C"){num=1}
    if ( obs$element[i]=="N"){num=0}
    if ( obs$element[i]=="S"){num=2}
  }
  if (obs$tissue[i]=="whole body"& obs$type[i]=="ectotherm" ){
    color="red4"
    if ( obs$element[i]=="C"){num=16}
    if ( obs$element[i]=="N"){num=15}
    if ( obs$element[i]=="S"){num=17}
  }
  if (obs$tissue[i]=="muscle" & obs$type[i]=="ectotherm"){
    color="red"
    if ( obs$element[i]=="C"){num=16}
    if ( obs$element[i]=="N"){num=15}
    if ( obs$element[i]=="S"){num=17}
  }
  points(log10(obs$t_1.2_.d.[i]),log10(obs$t1_2_pred[i]),col=color,pch=num)
}
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")
curve(0.62313+x*0.39952 , add = TRUE, col = "red")
curve(0.62313+x*0.39952-0.86854 , add = TRUE, col = "red4")
dev.off()

##fitting linear models lambda-DEB####

#NOT KEPT AIC HIGHER
# mod=lm(log10(obs$t1_2_pred)~log10(obs$t_1.2_.d.))
# summary(mod)## intercept no diff of 0
# simulateResiduals(fittedModel = mod, plot = T) ##ok with the conditions
# linearHypothesis(mod, "log10(obs$t_1.2_.d.) = 1") ##on teste si coeff significativement diff de 1, et oui.
# 



#glm with lambda_DEB~t.1.2*tissue
obs$tissue<- as.factor(obs$tissue)
mod_glm=glm(log10(obs$t1_2_pred)~log10(obs$t_1.2_.d.)*obs$tissue)
summary(mod_glm)## intercept no diff of 0
simulateResiduals(fittedModel = mod_glm, plot = T)
linearHypothesis(mod, "log10(obs$t_1.2_.d.) = 1") ##coeff diff of 1

plot(log10(obs$t_1.2_.d.),log10(obs$t1_2_pred),xlim=c(0,3),ylim=c(0,3),xlab="log10(obst1/2)",
     ylab="log10(predt1/2)",pch=16,main="obs vs pred lambda DEB dataset 1")
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")
curve(0.62313+x*0.39952 , add = TRUE, col = "red")
curve(0.56210+x*0.41320-0.86854 , add = TRUE, col = "red4")



##TaC predictions vs obs####
#predicting lambda using thomas and crowther equations
TaC=NULL
for (i in 1:nrow(obs)){
  if (obs$element[i]=="N"){TaC=c(TaC,TaC_t_1_2_N(Ww=obs$Mass_.g.[i],t=obs$temp_..C.[i]))}
  else {TaC=c(TaC,TaC_t_1_2_C(Ww=obs$Mass_.g.[i],t=obs$temp_..C.[i]))}
  
}


tiff("obs_vs_pred_TaC_OBS.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,3),ylim=c(-1,3),xlab="observed Half-life",
     ylab="T&C predicted Half-life")
for (i in 1:nrow(obs)){
  if (obs$tissue[i]=="whole body"& obs$type[i]=="endotherm" ){
    color="red4"
    if ( obs$element[i]=="C"){num=1}
    if ( obs$element[i]=="N"){num=0}
    if ( obs$element[i]=="S"){num=2}
  }
  if (obs$tissue[i]=="muscle"& obs$type[i]=="endotherm" ){
    color="red"
    if ( obs$element[i]=="C"){num=1}
    if ( obs$element[i]=="N"){num=0}
    if ( obs$element[i]=="S"){num=2}
  }
  if (obs$tissue[i]=="whole body"& obs$type[i]=="ectotherm" ){
    color="red4"
    if ( obs$element[i]=="C"){num=16}
    if ( obs$element[i]=="N"){num=15}
    if ( obs$element[i]=="S"){num=17}
  }
  if (obs$tissue[i]=="muscle" & obs$type[i]=="ectotherm"){
    color="red"
    if ( obs$element[i]=="C"){num=16}
    if ( obs$element[i]=="N"){num=15}
    if ( obs$element[i]=="S"){num=17}
  }
  points(log10(obs$t_1.2_.d.[i]),TaC[i],col=color,pch=num)
}
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")
curve( 1.07892 +x*0.32030, add = TRUE, col = "red")
curve( 1.07892 +x*0.32030-0.62580  , add = TRUE, col = "red4")
dev.off()

##fitting linear model T&C####

# linear model, does not work
mod_ThoC=lm(TaC~log10(obs$t_1.2_.d.))
summary(mod_ThoC) ##intercept diff of 0
simulateResiduals(fittedModel = mod_ThoC, plot = T) ##bof with conditions
linearHypothesis(mod_ThoC, "log10(obs$t_1.2_.d.) = 1") ##on teste si coeff significativement diff de 1, et oui.



#glm gaussian log link, Tac~t.1.2*tissue
mod_glm_ThoC = glm(TaC~log10(obs$t_1.2_.d.)*obs$tissue)
summary(mod_glm_ThoC)
simulateResiduals(fittedModel = mod_glm_ThoC, plot = T) #works well

plot(log10(obs$t_1.2_.d.),TaC,xlim=c(0,3),ylim=c(0,3),xlab="log10(obst1/2)",
     ylab="log10(predt1/2)",pch=16,main="obs vs pred T&C dataset 1")
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")
curve( 1.07892 +x*0.32030, add = TRUE, col = "red")
curve( 1.07892 +x*0.32030-0.62580  , add = TRUE, col = "red4")


#legend
tiff("legend_fig6.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c("Ectotherms","Endotherms","Muscle Tissue","Whole Body","C","N","S"),cex=1,col=c("black","black","red","red4","black","black","black")
       ,pch=c(16,1,16,16,16,15,17),box.lty=0)
dev.off()




#----FIGURE 8 (dataset 2)----

###size of individuals
f=0.8
Size<- size_est_corr(f=f,obs=review_new_bis,data=data)
review_new_bis$L_corr <- Size$Lobs


##prediction
review_new_bis$t_1.2_pred=pred_lam_half_life(f=f,obs=review_new_bis,data=data)

#####Lambda-DEB predictions vs obs#####


tiff("obs_vs_pred_DEB_rev.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,3),ylim=c(0,3),xlab="observed Half-life",
     ylab="λ-DEB predicted Half-life",pch=16)
for (i in 1:nrow(review_new_bis)){
  if (review_new_bis$tissue[i]=="whole body"& review_new_bis$type[i]=="endotherm" ){
    color="blue4"
    if ( review_new_bis$element[i]=="C"){num=1}
    if ( review_new_bis$element[i]=="N"){num=0}
    if ( review_new_bis$element[i]=="S"){num=2}
  }
  if (review_new_bis$tissue[i]=="muscle"& review_new_bis$type[i]=="endotherm" ){
    color="blue"
    if ( review_new_bis$element[i]=="C"){num=1}
    if ( review_new_bis$element[i]=="N"){num=0}
    if ( review_new_bis$element[i]=="S"){num=2}
  }
  if (review_new_bis$tissue[i]=="whole body"& review_new_bis$type[i]=="ectotherm" ){
    color="blue4"
    if ( review_new_bis$element[i]=="C"){num=16}
    if ( review_new_bis$element[i]=="N"){num=15}
    if ( review_new_bis$element[i]=="S"){num=17}
  }
  if (review_new_bis$tissue[i]=="muscle" & review_new_bis$type[i]=="ectotherm"){
    color="blue"
    if ( review_new_bis$element[i]=="C"){num=16}
    if ( review_new_bis$element[i]=="N"){num=15}
    if ( review_new_bis$element[i]=="S"){num=17}
  }
  points(log10(review_new_bis$t1.2[i]),log10(review_new_bis$t_1.2_pred[i]),col=color,pch=num)
}
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")
dev.off()

##fitting GLM lambda-DEB####
##glm loi gamma
mod_2=glm(review_new_bis$t_1.2_pred~review_new_bis$t1.2,
          family = Gamma(link="log"))
summary(mod_2)
simulateResiduals(fittedModel = mod_2, plot = T)

plot(log10(review_new_bis$t1.2),log10(review_new_bis$t_1.2_pred),xlim=c(0,3),ylim=c(0,3),xlab="log(obst1/2)",
     ylab="log(predt1/2)",pch=16)
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")
curve( 2.350241/log(10) +x*0.014339/log(10)  , add = TRUE, col = "red")

curve(coef(mod_2)[1]+coef(mod_2)[2]*x,add=TRUE)


mod_3=glm(log10(review_new_bis$t_1.2_pred)~log10(review_new_bis$t1.2),
          family = Gamma)
summary(mod_2)
simulateResiduals(fittedModel = mod_2, plot = T)


#####T&C predictions vs obs#####
Tac_rev=TaC_t_1_2_C(Ww=review_new_bis$Mass..g.,t=review_new_bis$Body.temperature...C.)
# tiff("obs_vs_pred_TaC_rev.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(log10(review_new_bis$t1.2),Tac_rev,xlim=c(0,3),ylim=c(0,3),xlab="log10(obst1/2)",
     ylab="log10(predt1/2)",pch=16)
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")
# dev.off()

tiff("obs_vs_pred_TaC_rev.tiff",  width =4 , height = 5, units = "in", res = 300)
plot(NULL,xlim=c(0,3),ylim=c(0,3),xlab="observed Half-life",
     ylab="T&C predicted Half-life",pch=16)
for (i in 1:nrow(review_new_bis)){
  if (review_new_bis$tissue[i]=="whole body"& review_new_bis$type[i]=="endotherm" ){
    color="blue4"
    if ( review_new_bis$element[i]=="C"){num=1}
    if ( review_new_bis$element[i]=="N"){num=0}
    if ( review_new_bis$element[i]=="S"){num=2}
  }
  if (review_new_bis$tissue[i]=="muscle"& review_new_bis$type[i]=="endotherm" ){
    color="blue"
    if ( review_new_bis$element[i]=="C"){num=1}
    if ( review_new_bis$element[i]=="N"){num=0}
    if ( review_new_bis$element[i]=="S"){num=2}
  }
  if (review_new_bis$tissue[i]=="whole body"& review_new_bis$type[i]=="ectotherm" ){
    color="blue4"
    if ( review_new_bis$element[i]=="C"){num=16}
    if ( review_new_bis$element[i]=="N"){num=15}
    if ( review_new_bis$element[i]=="S"){num=17}
  }
  if (review_new_bis$tissue[i]=="muscle" & review_new_bis$type[i]=="ectotherm"){
    color="blue"
    if ( review_new_bis$element[i]=="C"){num=16}
    if ( review_new_bis$element[i]=="N"){num=15}
    if ( review_new_bis$element[i]=="S"){num=17}
  }
  points(log10(review_new_bis$t1.2[i]),Tac_rev[i],col=color,pch=num)
}
points(seq(0,3,by=0.1),seq(0,3,by=0.1),type="l")

dev.off()

##fitting GLM T&C####
#glm
mod_ThoC_2=glm(10^Tac_rev~review_new_bis$t1.2,
               family = Gamma(link = "log"))
summary(mod_ThoC_2)
simulateResiduals(fittedModel = mod_ThoC_2, plot = T)



#legend
tiff("legend_fig7.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c("Ectotherms","Endotherms","Muscle Tissue","Whole Body","C","N","S"),cex=1,col=c("black","black","blue","blue4","black","black","black")
       ,pch=c(16,1,16,16,16,15,17),box.lty=0)
dev.off()


#----FIGURE 9----
#respiration Oncorhynchus nerka
H_d_p=subset(data,data$Species=="Oncorhynchus_nerka")

#see development for calculations in mol/J
eta_OA= -(-H_d_p$n_CX*H_d_p$eta_XA+H_d_p$n_CE/H_d_p$mu_E+H_d_p$n_CP*H_d_p$eta_PA+
            0.25*(-H_d_p$n_HX*H_d_p$eta_XA+H_d_p$n_HE/H_d_p$mu_E+H_d_p$n_HP*H_d_p$eta_PA)-
            0.5*(-H_d_p$n_OX*H_d_p$eta_XA+H_d_p$n_OE/H_d_p$mu_E+H_d_p$n_OP*H_d_p$eta_PA)-
            ((4*H_d_p$n_CN+H_d_p$n_HN-2*H_d_p$n_ON)/4*H_d_p$n_NN)*(-H_d_p$n_NX*H_d_p$eta_XA+
                                                                     H_d_p$n_NE/H_d_p$mu_E+H_d_p$n_NP*H_d_p$eta_PA))

eta_OD= -( -H_d_p$n_CE/H_d_p$mu_E-0.25*H_d_p$n_HE/H_d_p$mu_E+0.5*H_d_p$n_OE/H_d_p$mu_E+
             ((4*H_d_p$n_CN+H_d_p$n_HN-2*H_d_p$n_ON)/4*H_d_p$n_NN)*H_d_p$n_NE/H_d_p$mu_E)

eta_OG=-(H_d_p$n_CV*H_d_p$eta_VG-H_d_p$n_CE/H_d_p$mu_E+
           0.25*(H_d_p$n_HV*H_d_p$eta_VG-H_d_p$n_HE/H_d_p$mu_E)-
           0.5*(H_d_p$n_OV*H_d_p$eta_VG-H_d_p$n_OE/H_d_p$mu_E)-
           ((4*H_d_p$n_CN+H_d_p$n_HN-2*H_d_p$n_ON)/4*H_d_p$n_NN)*(H_d_p$n_NV*H_d_p$eta_VG-
                                                                    H_d_p$n_NE/H_d_p$mu_E))

f=0.8
TC=1
size_Oncor = seq(H_d_p$L_b,H_d_p$L_i*f,length=100)
Ww_Oncor = Ww(L=size_Oncor,f=0.8,p_am=H_d_p$p_Am
              ,v=H_d_p$v,n_HE=H_d_p$n_HE,n_OE=H_d_p$n_OE,n_NE=H_d_p$n_NE,mu_E=H_d_p$mu_E)


##O2 consumption estimation per day
O2_conso=data.frame(row.names = NULL)
Sm_val=NULL
for (i in 1: length(size_Oncor)){
  # i=1
  L=size_Oncor[i]
  if (L<H_d_p$L_j){
    Sm=L/H_d_p$L_b
    E_H=E_H_func(La=H_d_p$L_b,Lc=L,Lb=H_d_p$L_j,E_Ha=H_d_p$E_Hb,E_Hb=H_d_p$E_Hj)}
  if (H_d_p$L_j<=L& L<H_d_p$L_p){
    Sm=H_d_p$s_M
    E_H=E_H_func(La=H_d_p$L_j,Lc=L,Lb=H_d_p$L_p,E_Ha=H_d_p$E_Hj,E_Hb=H_d_p$E_Hp)
  }
  if (L>=H_d_p$L_p){
    Sm=H_d_p$s_M
    E_H=H_d_p$E_Hp
  }
  conso=resp_func(f=f,p_am=H_d_p$p_Am,L=L
                  ,v=H_d_p$v,Sm=Sm,Eg=H_d_p$E_G,k=H_d_p$kap
                  ,p_m=H_d_p$p_M,TC=TC,k_j=H_d_p$k_J,E_H=E_H
                  ,k_r=H_d_p$kap_R,eta_OA=eta_OA,eta_OD=eta_OD,eta_OG=eta_OG)
  O2_conso=rbind(O2_conso,data.frame(L=L,O2=conso))
  Sm_val=c(Sm_val,Sm)
}

O2_conso$Ww<- Ww_Oncor[[2]]
##O2 consumption estimation per hour
plot(O2_conso$Ww,O2_conso$O2/24,xlab="Ww (g)",ylab="O2 conso µmol/h")
##O2 consumption estimation per day
plot(O2_conso$Ww,O2_conso$O2,xlab="Ww (g)",ylab="O2 conso µmol/d")

#### 
molar_weight = 16*H_d_p$n_OE+16*H_d_p$n_OV

O2_conso$mat_quant <- O2_conso$Ww/molar_weight

plot(O2_conso$Ww,O2_conso$O2/O2_conso$mat_quant,xlab="Ww (g)",ylab="O2 conso d-1")


#lambda estimation
lambda_H_d <- lambda(L=size_Oncor,f=f,Sm=Sm_val
                     ,TC=TC,p_am=H_d_p$p_Am,v=H_d_p$v
                     ,k=H_d_p$kap,p_m=H_d_p$p_M
                     ,Eg=H_d_p$E_G,p_T=H_d_p$p_T)
lambda_H_d$Ww<-Ww_Oncor[[2]]
points(lambda_H_d$Ww,lambda_H_d$lam,xlab="Ww (g)",ylab="lambda (d-1)")

##FIG9A####
#comparing the ratio per day of the two 
ratio<- (O2_conso$O2/O2_conso$mat_quant)/lambda_H_d$lam

tiff("ratio_ONCOR_O2_lambda.tiff",  width =4 , height = 5, units = "in", res = 300)
par(mar=c(4,4,3,5))
plot(lambda_H_d$Ww,ratio,xlab="Wet Weight (g)",ylab="ratio conso O2/λ ",
     type="l",lwd=2)
dev.off()


##FIG 9B#### 
tiff("O2_lambda_Oncor_n_3.tiff",  width =4 , height = 5, units = "in", res = 300)
par(mar=c(4,4,3,5))  
plot(O2_conso$Ww,log10(O2_conso$O2/O2_conso$mat_quant),xlab="",ylab="",axes=F, type="l",col="black",lwd=2)


# legend axix lambda on left
axis(2, ylim=c(0,1),col="black")
mtext("O consumption ",side=2,line=2.5)
box() 
par(new=T) 
plot(O2_conso$Ww, log10(lambda_H_d$lam), pch=15,  xlab="", ylab="", axes=F, type="l", col="grey50",lwd=2)

# Legend axis Y on right
mtext("λ ",side=4,col="grey50",line=2.5)
axis(4, ylim=c(0,3),col="grey50",col.axis="grey50")

# AXE X
axis(1,pretty(range(lambda_H_d$Ww),10)) 
mtext("Wet weight (g)",side=1,col="black",line=2.5)
dev.off()
#----FIGURE 4----
H_d_p=subset(data,data$Species=="Oncorhynchus_nerka")

#max wet weight at f=1
Max_Ww_at_f1=Ww(L=Size_Hd,f=1,p_am=H_d_p$p_Am
                ,v=H_d_p$v,n_HE=H_d_p$n_HE
                ,n_OE=H_d_p$n_OE,n_NE=H_d_p$n_NE
                ,mu_E=H_d_p$mu_E)[[2]]

#running the graph
f=0.5
temp=10
TC=temp_correction(T_A=H_d_p$T_A,T_ref=H_d_p$T_ref,t=(temp+273.15))
Size_Hd= seq(H_d_p$L_b,H_d_p$L_i*f,length=100)


Max_Ww_at_f1=Ww(L=H_d_p$L_i,f=1,p_am=H_d_p$p_Am
                ,v=H_d_p$v,n_HE=H_d_p$n_HE
                ,n_OE=H_d_p$n_OE,n_NE=H_d_p$n_NE
                ,mu_E=H_d_p$mu_E)[[2]]



Sm_val=NULL
for (i in 1: length(Size_Hd)){
  # i=1
  L=Size_Hd[i]
  if (L<H_d_p$L_j){
    Sm=L/H_d_p$L_b}
  else{Sm=H_d_p$s_M}
  Sm_val=c(Sm_val,Sm)
}
# Sm_val=H_d_p$s_M

lambda_H_d <- lambda(L=Size_Hd,f=f,Sm=Sm_val
                     ,TC=TC,p_am=H_d_p$p_Am,v=H_d_p$v
                     ,k=H_d_p$kap,p_m=H_d_p$p_M
                     ,Eg=H_d_p$E_G,p_T=H_d_p$p_T)
lambda_H_d$Ww<- Ww(L=Size_Hd,f=f,p_am=H_d_p$p_Am
                   ,v=H_d_p$v,n_HE=H_d_p$n_HE
                   ,n_OE=H_d_p$n_OE,n_NE=H_d_p$n_NE
                   ,mu_E=H_d_p$mu_E)[[2]]


title=paste0("lambda_relation_On_f",f,"_TC_",temp,".tiff")
tiff(title,  width =4 , height = 5, units = "in", res = 300)
plot(lambda_H_d$Ww,lambda_H_d$lam
     ,xlim=c(0,Max_Ww_at_f1),ylim=c(0,0.3),xlab="Wet Weight (g)",ylab="λ ",type="l", lwd=2)
points(lambda_H_d$Ww,lambda_H_d$lamE,col="blue",type="l",lwd=2)
points(lambda_H_d$Ww,lambda_H_d$lamV,col="green",type="l",lwd=2)
dev.off()

#legend
tiff("legend_fig4.tiff",  width =4 , height = 5, units = "in", res = 300)
par(xpd=TRUE, mar=c(8,4,4,3))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c("λ ","λ ","λ"),cex=1,col=c("black","blue","green")
       ,lty=1,box.lty=0)
dev.off()

#----dataset 2 estimating wet weight-------


#Lymnaea_stagnalis Wet weight

param=subset(data,data$Species=="Lymnaea_stagnalis")
Wet_weight=Ww(L=0.1*0.341,f=1,p_am=param$p_Am
              ,v=param$v,n_HE=param$n_HE,n_OE=param$n_OE
              ,n_NE=param$n_NE,mu_E=param$mu_E)
#Engraulis_japonicus wet weigth
param=subset(data,data$Species=="Engraulis_japonicus")
Lw=c(9.05,9.5,9.7)*0.43546
Wet_weight=Ww(L=Lw,f=1,p_am=param$p_Am
              ,v=param$v,n_HE=param$n_HE,n_OE=param$n_OE
              ,n_NE=param$n_NE,mu_E=param$mu_E)
