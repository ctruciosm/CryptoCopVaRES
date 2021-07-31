##### Rolling Windows: VaR 7 cryptocurrencias
##### Carlos Trucios et al. (2020)
##### GAS

#### Load packages:
library(GAS)
library(VineCopula)
library(rugarch)
library(CDVine)


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
residuals = precopula = matrix(0,ncol=7,nrow=L) 
sigmafit = matrix(0,ncol=7,nrow=L+1)
paramsstd = paramES = matrix(0,nrow=7,ncol=3)
newresiduals = newreturns = matrix(0,ncol=7,nrow=M)
EWport = matrix(0,ncol=WR,nrow=M)
VaR = VaRT =  ES = EST = matrix(0,ncol=3,nrow=WR)
GASSpec_ST <- UniGASSpec(Dist = "std", GASPar = list(scale = TRUE)) # other option in the paper is ast1
mu = matrix(0,ncol=7,nrow = WR)
obport = sigmaport = c()
#### Windows
for (j in 1:WR){
  print(j)
  for (i in 1:7){
    mu[j,i] = mean(retornos[j:(j+L-1),i])
    y = scale(retornos[j:(j+L-1),i],center = TRUE,scale = FALSE)
    Modelo = UniGASFit(GASSpec_ST, y)
    sigmafit[,i] = Modelo@GASDyn$mTheta[2,]
    residuals[,i] = y/Modelo@GASDyn$mTheta[2,1:L]
    paramsstd[i,] = rugarch::fitdist(distribution = "std", residuals[,i])$pars
    precopula[,i] = pdist(distribution = "std", residuals[,i], mu = paramsstd[i,1], sigma = paramsstd[i,2], shape = paramsstd[i,3])
  }
  RMVC = CDVineCopSelect(precopula, type = "DVine") 
  simRVC = CDVineSim(10000, family =  RMVC$family, par = RMVC$par, par2 = RMVC$par2, type = "DVine") 
  for (k in 1:7){
    newresiduals[,k] = qdist(distribution = "std",  simRVC[,k], mu = paramsstd[k,1], sigma = paramsstd[k,2], shape = paramsstd[k,3])
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


write.table(data.frame(VaR),"VaRGASDVine_std.txt",sep=",",dec=".",col.names=FALSE,row.names = FALSE)
write.table(data.frame(ES),"ESGASDVine_std.txt",sep=",",dec=".",col.names=FALSE,row.names = FALSE)

