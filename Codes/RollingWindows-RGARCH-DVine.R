##### Rolling Windows: VaR 7 cryptocurrencias
##### Carlos Truc??os and Aviral K. Tiwari
##### GAS-Skew-T

#### Load packages:
library(GAS)
library(VineCopula)
library(CDVine)
library(rugarch)
library(Rsolnp)
library(metRology)
library(RobGARCHBoot)

#### Load data: BTC, DASH, DGB, DOGE, ITC, MAID, VTC
#BTC = read.csv("/Users/ctruciosm/Dropbox/Academico/VineCopulas/all/btc.csv")
BTC = read.csv("Data/btc.csv")
#which(BTC$date == "2015-01-01")
BTC = BTC$price.USD.[2190:3815]
DASH = read.csv("Data/dash.csv")
#which(DASH$date == "2015-01-01")
DASH = DASH$price.USD.[348:1973]
DGB = read.csv("Data/dgb.csv")
#which(DGB$date == "2015-01-01")
DGB = DGB$price.USD.[357:1982]
DOGE = read.csv("Data/doge.csv")
#which(DOGE$date == "2015-01-01")
DOGE = DOGE$price.USD.[390:2015]
ITC = read.csv("Data/ltc.csv")
#which(ITC$date == "2015-01-01")
LTC = ITC$price.USD.[1183:2808]
MAID = read.csv("Data/maid.csv")
#which(MAID$date == "2015-01-01")
MAID = MAID$price.USD.[255:1880]
VTC = read.csv("Data/vtc.csv")
#which(VTC$date == "2015-01-01")
VTC= VTC$price.USD.[357:1982]

crypto = data.frame(BTC, DASH, DGB, DOGE, LTC, MAID, VTC) 
retornos = (log(crypto[2:1626,])-log(crypto[1:1625,]))*100



#### General Setting
Sigma = c()
N = length(retornos[,1])
L = 750
WR = N-L
M = 10000
residuals = precopula = precopulapobs = matrix(0,ncol=7,nrow=L) 
sigmafit = matrix(0,ncol=7,nrow=L+1)
paraminnov = matrix(0,nrow=7,ncol=3)
newresiduals = newreturns = matrix(0,ncol=7,nrow=M)
EWport = matrix(0,ncol=WR,nrow=M)
VaR = VaRT = ES = EST = matrix(0,ncol=3,nrow=WR)
mu = matrix(0,ncol=7,nrow = WR)
obport = sigmaport = c()
#### Windows
for (j in 1:WR){
  print(j)
  for (i in 1:7){
    mu[j,i] = mean(retornos[j:(j+L-1),i])
    y = scale(retornos[j:(j+L-1),i],center = TRUE,scale = FALSE)
    coeff = ROBUSTGARCH(y) 
    sigmafit[,i] = fitted_Vol(coeff, y)
    residuals[,i] = y/sigmafit[1:L,i]
    # Other distributions were also used in the paper
    paraminnov[i,] = rugarch::fitdist(distribution = "std", residuals[,i])$pars
    precopula[,i] = pdist(distribution = "std", residuals[,i], mu = paraminnov[i,1], sigma = paraminnov[i,2] ,shape = paraminnov[i,3])
  }
  RMVC = CDVineCopSelect(precopula, type = "DVine") 
  simRVC = CDVineSim(10000, family =  RMVC$family, par = RMVC$par, par2 = RMVC$par2, type = "DVine") 
  for (k in 1:7){
    newresiduals[,k] = qdist(distribution = "std",  simRVC[,k], mu = paraminnov[k,1], sigma = paraminnov[k,2], shape = paraminnov[k,3])
    newreturns[,k] = sigmafit[L+1,k]*newresiduals[,k] + mu[j,k]
  }
  R = cor(newreturns)
  H = diag(apply(newreturns,2,sd))%*%R%*%diag(apply(newreturns,2,sd))
  ones = rep(1,7)
  sigmaport[j] = sqrt(t(ones)%*%H%*%ones)
  EWport[,j] = apply(newreturns,1,mean)
  VaR[j,] = quantile(EWport[,j],c(0.01,0.025,0.05))  
  for (l in 1:3){
    ES[j,l] = mean(EWport[EWport[,j]<VaR[j,l],j])
  }
  obport[j] = mean(as.numeric(retornos[j+L,]))
}


write.table(data.frame(VaR),"VaRRGARCH2DVine2_std.txt",sep=",",dec=".",col.names=FALSE,row.names = FALSE)
write.table(data.frame(ES),"ESRGARCH2DVine2_std.txt",sep=",",dec=".",col.names=FALSE,row.names = FALSE)
