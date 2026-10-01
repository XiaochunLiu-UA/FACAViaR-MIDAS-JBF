 rm(list = ls())


################################################
######### Monthly Quantile Factors
FredMD <- read.csv("FRED_MD_Clean.csv")

dates <- FredMD[,1]
dat<- FredMD[,-c(1,2)]

vnam <- names(dat)
#####################################################################################

TAU = c(0.01,0.05,0.1,0.5,0.9,0.95,0.99) # compute q quantile positions

# Matlab variables
st <- which(dates=="2004-12-01") 
totalObs <-  dim(dat)[1] 

outSampleSize <- totalObs-st+1


rmax = 3
r = rmax          # number of factors in simulation
 
 
 CDG_QFactors <- list()
 CDG_Loadings <- list()
  
 for (i in 1:outSampleSize) {
 datt <- dat[1:(i+st-1),]
 
 CDG_QFactors[[i]] <- array(0,c(dim(datt)[1],r,length(TAU)),
        dimnames=list(dates[1:(i+st-1)],c("CQF(1)" ),
 CDG_Loadings[[i]] <- array(0,c(r,dim(datt)[2],length(TAU)),
        dimnames=list(c("CQF(1)" ),vnam,   TAU))
 
  for (j in 1:length(TAU)) {
  print(paste(i,"-",TAU[j]))
         tau <- TAU[j] 
         fname = paste0("lpfactor_",tau,"_",i,".csv")
         res <- read.csv(fname,header=FALSE)
        CDG_QFactors[[i]][,,j] <- as.matrix(res)
        
        fname = paste0("lpfactor_Loadings_",tau,"_",i,".csv")
         res <- read.csv(fname,header=FALSE)
        CDG_Loadings[[i]][,,j] <- as.matrix(res)
       }
 
}
 
 save(CDG_QFactors,file="CDG_QFactors3_Monthly.Rdata")
 save(CDG_Loadings,file="CDG_Loadings3_Monthly.Rdata")
 