
#----FUNCTIONS----

#estimating lambda
lambda <- function(L,f,Sm,TC,p_am,v,k,p_m,Eg,p_T){
  lambdaE=(v*Sm*TC)/L
  lambdaV=k*f*p_am*TC/(Eg*v)*(((v*Sm*Eg+p_T)*L^-1+p_m)/(k*f*p_am/v+Eg))-(p_m+p_T*L^-1)*TC/Eg
  lambda=lambdaE+lambdaV
  return(data.frame(L=c(L),lamE=c(lambdaE),lamV=c(lambdaV),lam=c(lambda)))
}

#temperature in K
temp_correction<-function(T_A,T_ref,t){
  return(exp((T_A/T_ref)-(T_A/t)))
}

# #wet weight
Ww<- function(L,f,p_am,v,n_HE,n_OE,n_NE,mu_E){
  We=12+n_HE+16*n_OE+14*n_NE
  Ww=(L)^3+(f*(p_am/v)*(L)^3)*(We/mu_E)
  return(data.frame(L=L,Ww=Ww))
}

##ww-> L structural
St_length <- function(Ww,f,p_am,v,n_HE,n_OE,n_NE,mu_E){
  We=12+n_HE+16*n_OE+14*n_NE
  denom=1+(f*p_am*We)/(v*mu_E)
  return((Ww/denom)^(1/3))
}



##Thomas and Crowther log10(t1/2) Ww relations  
#Carbon
TaC_t_1_2_C <- function(Ww,t){
  return(1.668+0.1935*log10(Ww)-0.0153*t)
}
#Nitrogen
TaC_t_1_2_N <- function(Ww,t){
  return(1.664+0.1933*log10(Ww)-0.0149*t)
}


#vander Zanden t1/2 Ww relations 
#val is the value 3.28 for invertebrates and 3.65 for vertebrates
VZ_ecto <- function(Ww,val ) {return(0.22*log10(Ww)+val/log(10))}
VZ_endo <-function(Ww){return(0.13*log10(Ww)+2.01/log(10))}

##HERE IS THE LOOP TO GET LAMBDA WET WEIGHT RELATIONSHIP FOR EACH SPECIES WE HAVE on obs
#num_points is the number of size tested between Lb and Li
#f_val a seq of tested f values between 0 and 1
#species name a c() of tested name species
#data the DEB database

relashionship <- function(num_points,f_val,species_names,data){
  rel_lam_ww <- data.frame(row.names = NULL)
  
  for (i in 1:length(species_names)){ #for every species
    # i=19
    name=species_names[i]
    DEB_param=subset(data,data$Species==name)
    TC=1
    for (k in 1:length(f_val)){ #for every f
      # k=6
      f=f_val[k]
      
      #TESTED ST SIZE. IF  std, stf, abj,asj sizes beteween Lb and Li; if stx, size between Lx and Li
      if (DEB_param$Model=="stx") {size=seq(DEB_param$L_x,DEB_param$L_i*f,length=num_points)
      }
      else {size=seq(DEB_param$L_b,DEB_param$L_i*f,length=num_points)
      }
      
      
      
      if (DEB_param$Model=="std"|DEB_param$Model=="stf"|DEB_param$Model=="stx"){ #for std, stf, stx
        lambda_value<- lambda(L=size
                              ,f=f,Sm=DEB_param$s_M,TC=TC,p_am=DEB_param$p_Am,
                              v=DEB_param$v,k=DEB_param$kap,p_m=DEB_param$p_M,Eg=DEB_param$E_G,p_T=DEB_param$p_T)
      }
      if (DEB_param$Model=="abj") { #for abj 
        
        lambda_value=data.frame(row.names = NULL)
        Sm= data.frame(row.names = NULL)
        for (j in 1:length(size)){
          if (size[j]<DEB_param$L_j){
            Sm=rbind(Sm,data.frame(Sm=size[j]/DEB_param$L_b))
            lambda_value=rbind(lambda_value,lambda(L=size[j]
                                                   ,f=f,Sm=size[j]/DEB_param$L_b,TC=TC,p_am=DEB_param$p_Am,
                                                   v=DEB_param$v,k=DEB_param$kap,p_m=DEB_param$p_M,Eg=DEB_param$E_G,p_T=DEB_param$p_T))}
          else {Sm=rbind(Sm,data.frame(Sm=DEB_param$s_M)) 
          lambda_value=rbind(lambda_value,lambda(L=size[j]
                                                 ,f=f,Sm=DEB_param$s_M,TC=TC,p_am=DEB_param$p_Am,
                                                 v=DEB_param$v,k=DEB_param$kap,p_m=DEB_param$p_M,Eg=DEB_param$E_G,p_T=DEB_param$p_T))}
        }
      }
      
      if (DEB_param$Model=="asj") { #for asj 
        
        lambda_value=data.frame(row.names = NULL)
        for (j in 1:length(size)){
          # j=2
          if (size[j]<DEB_param$L_s){
            Sm=1
            lambda_value=rbind(lambda_value,lambda(L=size[j]
                                                   ,f=f,Sm=Sm,TC=TC,p_am=DEB_param$p_Am,
                                                   v=DEB_param$v,k=DEB_param$kap,p_m=DEB_param$p_M,Eg=DEB_param$E_G,p_T=DEB_param$p_T))
            
          }
          if (size[j]>DEB_param$L_s & size[j]<DEB_param$L_j){
            Sm=size[j]/DEB_param$L_s
            lambda_value=rbind(lambda_value,lambda(L=size[j]
                                                   ,f=f,Sm=size[j]/DEB_param$L_b,TC=TC,p_am=DEB_param$p_Am,
                                                   v=DEB_param$v,k=DEB_param$kap,p_m=DEB_param$p_M,Eg=DEB_param$E_G,p_T=DEB_param$p_T))
            
          }
          if (size[j]>=DEB_param$L_j) {
            lambda_value=rbind(lambda_value,lambda(L=size[j]
                                                   ,f=f,Sm=DEB_param$s_M,TC=TC,p_am=DEB_param$p_Am,
                                                   v=DEB_param$v,k=DEB_param$kap,p_m=DEB_param$p_M,Eg=DEB_param$E_G,p_T=DEB_param$p_T))
          }
        }
      }
      
      lambda_value$t1_2<- log(2)/lambda_value$lam
      Wet_w=Ww(L=size,f=f,p_am=DEB_param$p_Am,v=DEB_param$v
               ,n_HE=DEB_param$n_HE,n_OE=DEB_param$n_OE,n_NE=DEB_param$n_NE,mu_E=DEB_param$mu_E)
      lambda_value$Ww<- Wet_w$Ww
      
      rel_lam_ww= rbind(rel_lam_ww ,data.frame(species=rep(name,num_points),f=rep(f,num_points),L=lambda_value$L,
                                               Ww=lambda_value$Ww,lambda=lambda_value$lam,
                                               t1_2=lambda_value$t1_2))
    }
  }
  return(rel_lam_ww)
}


#data the DEB database
#f is the value tested
#obs is the hlaf-life observation database (either oth or review_new)
size_est_corr <- function(f,obs,data){
  L_obs=data.frame(row.names = NULL)
  compting=0
  for (i in 1:nrow(obs)){
    # i=1
    name=obs$Species[i]
    DEB_param =subset(data,data$Species==name)
    obs_line = obs[i,]
    L_curr=St_length(Ww=obs_line$Mass,f=f,p_am=DEB_param$p_Am
                     ,v=DEB_param$v,n_HE=DEB_param$n_HE
                     ,n_OE=DEB_param$n_OE,n_NE=DEB_param$n_NE,mu_E=DEB_param$mu_E)
    
    if(L_curr>DEB_param$L_i){
      L_corr=DEB_param$L_i
      # compting=compting+1
    }
    else {L_corr=L_curr}
    
    L_obs=rbind(L_obs,data.frame(Lobs=L_curr,Li=DEB_param$L_i,Lx=DEB_param$L_x,Lb=DEB_param$L_b,L_corr=L_corr))
  }
  return (L_obs)
}

#Half life predictions using lambda-DEB
######/!\ in obs there must be a TC column
#/!\ in obs there must be a L_corr column with the size
pred_lam_half_life <- function(f,obs,data){
  t_1_2_pred=NULL
  for (i in 1:nrow(obs)){
    # i=1
    name=obs$Species[i]
    DEB_param =subset(data,data$Species==name)
    obs_line = obs[i,]
    #SM value
    if (DEB_param$Model=="std"|DEB_param$Model=="stf"|DEB_param$Model=="stx"){Sm=1}
    if(DEB_param$Model=="abj") {
      if(obs_line$L_corr<DEB_param$L_j){Sm=obs_line$L_corr/DEB_param$L_b}
      else {Sm=DEB_param$s_M}
    }
    if (DEB_param$Model=="asj"){
      if (obs_line$L_corr<DEB_param$L_s){Sm=1}
      if(obs_line$L_corr>DEB_param$L_s|obs_line$L_corr<DEB_param$L_j)
      {Sm=obs_line$L_corr/DEB_param$L_s}
      if (obs_line$L_corr>DEB_param$L_j){Sm=DEB_param$s_M}
      
    }
    
    pred_lam=lambda(L=obs_line$L_corr,f=f,Sm=Sm,TC=obs_line$TC
                    ,p_am=DEB_param$p_Am,v=DEB_param$v,k=DEB_param$kap
                    ,p_m=DEB_param$p_M,Eg=DEB_param$E_G,p_T=DEB_param$p_T)
    
    t_1_2_pred=c(t_1_2_pred,log(2)/pred_lam$lam)
    
  }
  return(t_1_2_pred)
}

###Respirations functions###

#function to interpolate E_H according to L. (linear interpolation)
#La <Lc <Lb
#we knwo E_Ha and E_Hb and looking fo E_Hc
E_H_func <- function(La,Lc,Lb,E_Ha,E_Hb){
  E_Hc= (1-abs(Lc-La)/abs(Lb-La))*E_Ha+(1-abs(Lb-Lc)/abs(Lb-La))*E_Hb
  return(E_Hc)
}

#fluxes functions
p_A_func <- function(p_am,Sm,f,TC,L){
  return (p_am*Sm*f*L^2*TC)
}

p_S_func <-function(p_m,L,TC){return(p_m*L^3*TC)}
p_J_func<- function(k_j,E_H){return(k_j*E_H)}

p_C_func <- function(f,p_am,L,v,Sm,Eg,k,p_m,TC) {
  E=f*p_am*L^3/v
  p_C=E*((v*TC*Sm*L^2*Eg+p_S_func(p_m,L,TC))/(E+Eg*L^3))
  return(p_C)
}

p_G_func<- function(f,p_am,L,v,Sm,Eg,k,p_m,TC){
  p_G=k*p_C_func(f,p_am,L,v,Sm,Eg,k,p_m,TC)-p_S_func(p_m,L,TC)
  return(p_G)
}

p_R_func<- function(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H){
  p_R=(1-k)*p_C_func(f,p_am,L,v,Sm,Eg,k,p_m,TC)-p_J_func(k_j,E_H)
  return(p_R)
}
p_D_func <- function(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H,k_r){
  p_D=p_S_func(p_m,L,TC)+p_J_func(k_j,E_H)+(1-k_r)*p_R_func(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H)
  return(p_D)
}


###respiration function in mol/J/d
resp_func<- function(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H,k_r,eta_OA,eta_OD,eta_OG){
  resp =eta_OA*p_A_func(p_am,Sm,f,TC,L)+eta_OD*p_D_func(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H,k_r)+
    eta_OG*p_G_func(f,p_am,L,v,Sm,Eg,k,p_m,TC)
  return(resp)
}

###showing the fluxes
pACSJGRD<-function(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H,k_r){
  fluxes= data.frame(pA=p_A_func(p_am,Sm,f,TC,L),
                     pC=p_C_func(f,p_am,L,v,Sm,Eg,k,p_m,TC),
                     pS=p_S_func(p_m,L,TC),pJ=p_J_func(k_j,E_H),
                     pG=p_G_func(f,p_am,L,v,Sm,Eg,k,p_m,TC),
                     pR=p_R_func(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H),
                     pD=p_D_func(f,p_am,L,v,Sm,Eg,k,p_m,TC,k_j,E_H,k_r))
  return(fluxes)
}