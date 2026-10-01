################

library(fbi)
library(writexl)




###### monthly

dat <-  fredmd("2026-02.csv")

N <- dim(dat)[2]


dat1 <- tp_apc(X=dat[,-1],kmax=8)

dat2 <- dat1$data

dat3 <- dat
dat3[,2:N] <- dat2


x <- dat3[,"INDPRO"]
mn <- min(x)
mx <- max(x)
md <- median(x)
iqr <- IQR(x)
Lout <- abs((mn-md)/iqr)
Uout <- (mx-md)/iqr
if(Lout>=10) {
 loc <- which(x==mn)
 print(dat3[loc,1])
 print(dat3[loc,"INDPRO"])
}
### McCraken and Ng (2016) approach remove this outlier.
### However, it is inappropriate, because it is 4/2020 during COVID-19 important for risk modeling.
ro <- abs((x-md)/iqr)
sum(ro>=10)



################ after clean the data: 1/1959 to 1/2026

FRED_MD_Clean <- dat3


################################ arrange major macro variables
### in rolling window format
dates <- FRED_MD_Clean[,1] 
st <- which(dates=="2004-12-01") 
totalObs <-  dim(FRED_MD_Clean)[1]
outSampleSize <- totalObs-st+1
Macro_MD_RealTime <- list()

for (i in st:totalObs) {
  xx<-FRED_MD_Clean[1:i,c("INDPRO","CPIAUCSL","TB3MS")]
  Macro_MD_RealTime[[i-st+1]] <- xx
  dimnames(Macro_MD_RealTime[[i-st+1]])[[1]] <-FRED_MD_Clean[1:i,"date"] 
  
 }


dates <- as.Date(dat3[,1])
year_numeric <- as.numeric(format(dates, "%Y"))
month_numeric <- as.numeric(format(dates, "%m"))
ym <- year_numeric+(month_numeric-1)/12

FRED_MD_Clean <- cbind(dat3[,1],ym,dat3[,-1])
dimnames(FRED_MD_Clean)[[2]][1]<- "Dates"


save(FRED_MD_Clean,file="FRED_MD_Clean.Rdata")
save(Macro_MD_RealTime,file="Macro_MD_RealTime.Rdata")

write.csv(FRED_MD_Clean,file="FRED_MD_Clean.csv",row.names=F)
write_xlsx(FRED_MD_Clean, "FRED_MD_Clean.xlsx")


###############################################
################################

dat4 <- rm_outliers.fredmd(dat)
dat5 <- tp_apc(X=dat4[,-1],kmax=8)
dat6 <- dat5$data
dat7 <- dat
dat7[,2:N] <- dat6

x <- dat7[,"INDPRO"]
mn <- min(x)
mx <- max(x)
md <- median(x)
iqr <- IQR(x)
Lout <- abs((mn-md)/iqr)
Uout <- (mx-md)/iqr
if(Lout>=10) {
 loc <- which(x==mn)
 print(dat7[loc,1])
 print(dat7[loc,"INDPRO"])
}

### 4/2020 data is removed if using rm_outliers.fredmd.
## Apparently, this is not correct.



