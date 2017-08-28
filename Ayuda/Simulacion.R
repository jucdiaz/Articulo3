if (!require('devtools')) install.packages('devtools')
if (!require('rmutil')) install.packages('rmutil')
if (!require('boot')) install.packages('boot')
if (!require('numDeriv')) install.packages('numDeriv')
if (!require('GHQp')) install.packages('GHQp')
devtools::install_github('jucdiaz/ZOIP', force=TRUE)
library(ZOIP)

family<-'R-S'
N<-10
n<-c(5,10,20,30,50,75)
n.points<-c(3,5,10,15,25)
pruning<-c(TRUE,FALSE)

escenarios<-expand.grid(N,n,n.points,pruning)
colnames(escenarios)<-c('N','ni','n.points','pruning')
comb<-c(19,53,27)
escenarios<-escenarios[comb,]
dim(escenarios)
nrep<-1000


mu.b0<-NULL
mu.b1<-NULL
sigma.b0<-NULL
sigma.b1<-NULL
p0.b0<-NULL
p1.b0<-NULL
random.mu<-NULL
random.sigma<-NULL
logvero<-NULL
time<-NULL
num.iter<-NULL


formula.mu=Y~Total_mora
formula.sigma=~Total_mora
formula.p0=~1
formula.p1=~1
formula.random= ~ 1 | Ciudad

link=c('logit','logit','identity','identity')

optimizer<-'nlminb'

#-------------------------------------------------------------------
i<-1
while(i<=dim(escenarios)[1]){
  j<-1
  while(j<=nrep){
    try({ 

Ciudad <- rep(1:escenarios[i,1], each=escenarios[i,2])
Total_mora<-rexp(escenarios[i,2]*escenarios[i,1],rate=1.075)

b0i <- rep(rnorm(n=escenarios[i,1],sd=0.51), each=escenarios[i,2])
b1i <- rep(rnorm(n=escenarios[i,1],sd=0.4), each=escenarios[i,2])


neta <- (-1.13+b0i)+0.33*Total_mora
neta2<-(0.33+b1i)+0.14*Total_mora

mu <- 1 / (1 + exp(-neta))
sigma <- 1 / (1 + exp(-neta2))

p0 <- 0.23
p1 <- 0.07

mu[mu==1] <- 0.999 
mu[mu==0] <- 0.001 

sigma[sigma==1] <- 0.999
sigma[sigma==0] <- 0.001



Y <- rZOIP(n=length(mu), mu = mu, sigma = sigma ,p0=p0,p1=p1,family=family) 

data_sim<-data.frame(Y,Total_mora,Ciudad)

n.points <-escenarios[i,3]
pruning <- escenarios[i,4]

mod<-RMM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,
              formula.p1=formula.p1,data=data_sim,formula.random=formula.random,link=link,
              family=family,optimizer=optimizer,n.points=n.points,pruning=pruning)

mu.b0[j]<-mod$Fixed_Parameters.mu[1]
mu.b1[j]<-mod$Fixed_Parameters.mu[2]

sigma.b0[j]<-mod$Fixed_Parameters.sigma[1]
sigma.b1[j]<-mod$Fixed_Parameters.sigma[2]


p0.b0[j]<-mod$Fixed_Parameters.p0[1]

p1.b0[j]<-mod$Fixed_Parameters.p1[1]

random.mu[j]<-mod$Parameters.randoms[1,1]
random.sigma[j]<-mod$Parameters.randoms[2,1]

logvero[j]<-mod$logverosimilitud

time[j]<-mod$Time

num.iter[j]<-mod$num.iter

j=j+1
    },silent=TRUE)
  }
  
data_i<-as.data.frame(cbind(I=rep(i,nrep),N=rep(escenarios$N[i],nrep),ni=rep(escenarios$ni[i],nrep),
                            n.points=rep(escenarios$n.points[i],nrep),
                            pruning=rep(escenarios$pruning[i],nrep),
                            mu.b0,mu.b1,sigma.b0,sigma.b1,p0.b0,p1.b0,random.mu,random.sigma,
                            logvero,time,num.iter))
write.table(data_i,file='D:\\ZOIP_mix_JUAN.csv',append=TRUE,sep = ",",col.names = FALSE, row.names = FALSE)
i=i+1
}
i
j
escenarios[i,]
