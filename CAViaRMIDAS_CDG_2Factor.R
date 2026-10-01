rm(list=ls())


source("CAViaR_MIDAS_Joint_VaR_ES_Rcpp.R")
load("DailyData.Rdata")

 

dates <- dimnames(Data)[[1]]
nrnames <- dimnames(Data)[[2]]
TT <- dim(Data)[1]
nr <- dim(Data)[2]


Models <- c("SAV","AS")
LongRun <- c("exponential" , "linear")
#LModels <- c(paste(Models[1],LongRun),paste(Models[2],LongRun))
nm <- length(Models) 
nl <- length(LongRun)  
#ml <- length(LModels)
 
APA <- c(0.01,0.025,0.05) ## VaR Probabilities
na <- length(APA)

QFL <- c(0.01,0.05,0.1,0.5,0.9,0.95,0.99)  
nq <- length(QFL)


ddts <- format(as.Date(dates),"%Y-%m")
dind <- matrix(0,length(dates),1)
dind[1] <- 1

for (i in 2:length(dates)) {
       if (ddts[i]==ddts[i-1]) {
             dind[i] <- dind[i-1]
           } else {
              dind[i] <- dind[i-1]+1
           }
}


Datta <- cbind(dind,Data)  ## The last month is 01/2026 for daily returns
dimnames(Datta)[[2]][1] <- "ind"
### matrix for backward macrodata
R <- 120              ## start at month 120

mm <- which(Datta[,1]==R)
loc <- mm[length(mm)]
sz <- dim(Datta)[1]
fp <- sz-loc

fpnames <- dimnames(Datta)[[1]][(loc+1):sz]

Datf <- Data[(loc+1):sz,]   #### forecast periods returns


##################### Load Monthly Macro Factors
load("CDG_QFactors2_Monthly.Rdata")  

Factors <- CDG_QFactors

nf <- length(Factors)
fdates <- dimnames(Factors[[nf]])[[1]]
mdts <- format(as.Date(fdates),"%Y-%m")
loc0 <- which(mdts=="1994-12")



K <- 2 ## number of factors

S <- 24 ## number of lags for weighting functions
trim <- S-1

TF <- FALSE ## TRUE: w1=1, w2 free; FALSE: w1 and w2 free



pnam <- param_names(model="SAV", K=K, fix_w1 =TF)
pnam_SAV <-c(pnam,paste(pnam,"_se"),"q_init","rq_loss") 
pnam <- param_names(model="AS", K=K, fix_w1 =TF)
pnam_AS <-c(pnam,paste(pnam,"_se"),"q_init","rq_loss") 
 
np_SAV <- length(pnam_SAV) #3+1
np_AS <- length(pnam_AS) #4+1

#############################################
############################################
OOS_CAViaR_MIDAS <- list(
       Y_OOS =Datf,
       MIDAS_Weights=array(0,c(fp,S,nq,na,nl,nm,nr),dimnames=list(fpnames,paste("Weights",1:S),QFL,APA,LongRun,Models,nrnames)),
       QES_Forecast = array(0,c(fp,na,2,nq,nl,nm,nr),dimnames=list(fpnames,APA,c("VaRF","ESF"),QFL,LongRun,Models,nrnames))
)


OOS_CAViaR_MIDAS_PARSSTAT <- list(
    SAV =  array(0,c(fp,np_SAV,nq,na,nl,nr),
      dimnames=list(fpnames,pnam_SAV,QFL,APA,LongRun,nrnames)), 
    AS =  array(0,c(fp,np_AS,nq,na,nl,nr),
      dimnames=list(fpnames,pnam_AS,QFL,APA,LongRun,nrnames))
)

######################### 
 for (r in seq_along(nrnames)) {  ## the assets :1:nr
   dat <- Datta[,c(1,r+1)]
 for (a in seq_along(APA)) { # <- 1    # 1:na     VaR levels
   tau <- APA[a]
   for (m in seq_along(Models)){
        model <- Models[m]
     for (n in seq_along(LongRun)){
          lr_type <- LongRun[n]
       for (p in seq_along(QFL)) {
                      
 for (s in 1: fp) {
   print(paste(nrnames[r],APA[a],QFL[p],s))
   ddt <- dat[s:(loc+s-1),]
   y <- ddt[,2]
   ddf_date <- dimnames(dat)[[1]][loc+s]
   ddf <- matrix(dat[loc+s,],1,,dimnames=list(ddf_date,c("month_id",nrnames[r]) ))
   
   uq_ddts <- unique(ddts)
   ind <- ddt[,1]
   did0 <- ind[1]
   did1 <- ind[length(ind)]
   did1_prev <- ind[length(ind)-1]
   didf <- ddf[,1]
   start_dmonth<- uq_ddts[did0]
   last_dmonth <- uq_ddts[did1]
   prev_last_dmonth <- uq_ddts[did1_prev]
  
   month_id <- ind-did0+1 
   q_init <- quantile(y,prob=tau)
 
   if (s==1) {
    mloc0 <- which(mdts==start_dmonth)
    mloc1 <- which(mdts==last_dmonth)  
    ms <- mloc1-loc0-R+1
    factors <- Factors[[ms]][(mloc0-S):(mloc1-1),1:K,p]
    
    if (didf==did1) {
    factorsF <- factors
    } else {
    factorsF <- Factors[[ms]][(mloc0-S):mloc1,1:K,p]
    }
   }
   
   if (s >1) {
        if (did1!=did1_prev){
             mloc0 <- which(mdts==start_dmonth)
             mloc1 <- which(mdts==last_dmonth)
             ms <- mloc1-loc0-R+1-1
             factors <- Factors[[ms]][(mloc0-S):(mloc1-1),1:K,p]  
             factorsF <- factors
        }
        if (did1==did1_prev) {
             mloc0 <- which(mdts==start_dmonth)
             mloc1 <- which(mdts==last_dmonth)
             ms <- mloc1-loc0-R+1-1
             factors <- Factors[[ms]][(mloc0-S):(mloc1-1),1:K,p] 
            if (didf==did1) {
             factorsF <- factors
            } else {
             factorsF <- Factors[[ms+1]][(mloc0-S):mloc1,1:K,p]
            }
        }
   }
   

 

############ Model estimation 
 temp <- estimate_caviar_midas_es_cpp(y=y, factors=factors,
                                  month_id=month_id,                                                tau = tau,model =model,
                                  lr_type =lr_type,
                                  trim=trim,
                                  fix_w1 = TF,
                                  S = S, n0 = 300,
                                  n_init = 5000, m_best = 10,
                                  n_burn = 0, maxit = 5000,
                                  seed = 42,compute_se=FALSE)
                               
 y_last <- y[length(y)]
 vef <- forecast_caviar_midas_es(fit=temp, 
               y_last=y_last, factors_new=factorsF,trim=trim) 
               
OOS_CAViaR_MIDAS[["QES_Forecast"]][s,a,,p,n,m,r] <- c(vef$VaR_forecast,vef$ES_forecast)


OOS_CAViaR_MIDAS[["MIDAS_Weights"]][s,,p,a,n,m,r] <- temp$midas_weights
rq_loss <- quantile_loss(y=y, Q=temp$Q, tau=tau) 
OOS_CAViaR_MIDAS_PARSSTAT[[m]][s,,p,a,n,r]<- c(temp$theta,temp$se,temp$q_init,rq_loss$mean)

}
}
}
}
}
}

save(OOS_CAViaR_MIDAS,file=paste("OOS_CAViaR_MIDAS_CDG_",K,".Rdata"))
save(OOS_CAViaR_MIDAS_PARSSTAT,file=paste("OOS_CAViaR_MIDAS_PARSSTAT_CDG_",K,".Rdata"))



 
