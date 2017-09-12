library(ZOIP)
library(GHQp)
N<-21

Times <- c(2, 10, 20, 40) # cantidad de dias
Concentration.encoded <- rep(rep(-1:1, each=length(Times)), each=N/3) #3 por la cantidad de niveles en xij
subject <- rep(1:N, each=length(Times)) # numero de sujetos en la muestra repetidos tantas veces haya dias
Days <- rep(Times, times=N)
b0i <- rep(rnorm(n=N,sd=1), each=length(Times))
b1i <- rep(rnorm(n=N,sd=0.5), each=length(Times))

neta <- (1.6+b0i)-1.3*log(Days)#-0.85*(log(Days))^2+0.96*Concentration.encoded # Original 5
neta2<-(0.1+b1i)-0.8*log(Days) # 0.1 en 'R-S' y 'simplex' 8 en 'F-C' y 'Original' 
#neta3<-0.1-1.4*log(Days)
#neta4<-0.2-2.2*log(Days)
# hay 60 datos correspondientes a 15 datos por cada periodo de tiempo
mu <- 1 / (1 + exp(-neta)) # aplicamos funcion link
#mu<-exp(neta)
sigma <- 1 / (1 + exp(-neta2))
#sigma <-exp(neta2)
p0 <- 0.1#1 / (1 + exp(-neta3))
p1 <- 0.1#1 / (1 + exp(-neta4))

mu[mu==1] <- 0.999 # Transforma datos inflados en unos
mu[mu==0] <- 0.001 # transforma datos inflados en ceros

sigma[sigma==1] <- 0.999
sigma[sigma==0] <- 0.001
# 
# p0[p0==1] <- 0.999
# p0[p0==0] <- 0.001
# 
# p1[p1==1] <- 0.999
# p1[p1==0] <- 0.001


family<-'R-S'

Y <- rZOIP(n=length(mu), mu = mu, sigma = sigma ,p0=p0,p1=p1,family=family) # simulacion de datos con dispersion fija, gran detalle que la dispersion este fija (como hacer si no?)

base<-data.frame(Y,Days,Concentration.encoded,subject)


formula.mu <- Y ~ log(Days) #+ I(log(Days)^2) + Concentration.encoded
formula.sigma<-~log(Days)
formula.p0<-~1#log(Days)
formula.p1<-~1#log(Days)

formula.random=~1 | subject
#formula.random.sigma=~1 | subject
#it's very important that fulfilled formula.random.mu equal to formula.random.mu


#link<-c('logit','logit','logit','logit')
link<-c('logit','logit','identity','identity')


optimizer<-'nlminb'
n.points <-11
pruning <- TRUE

mod<-RMM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,formula.p1=formula.p1,data=base,
         formula.random=formula.random,link=link,family=family,optimizer=optimizer,
         n.points=n.points,pruning=pruning)

#-------------------------------------

RMM.ZOIP<-function(formula.mu,formula.sigma=~1,formula.p0=~1,formula.p1=~1,data,formula.random
                   ,link=c('identity','identity','identity','identity'),family='R-S',optimizer='nlminb',
                   n.points=11,pruning=TRUE){
  
  if (any(family != 'R-S') && any(family != 'F-C') && any(family != 'Original') && any(family != 'Simplex'))
    stop(paste("family must be in R-S, F-C, Original, Simplex", "\n", ""))
  if(any(length(as.character(attr(terms(formula.mu),'variable')))==1))
    stop(paste("formula.mu must have response variable","\n",""))
  if(any(attr(terms(formula.sigma),'response')>=1))
    stop(paste("sigma can't have response variable in formula.sigma","\n",""))
  if(any(attr(terms(formula.p0),'response')>=1))
    stop(paste("p0 can't have response variable in formula.p0","\n",""))
  if(any(attr(terms(formula.p1),'response')>=1))
    stop(paste("p1 can't have response variable in formula.p1","\n",""))
  if(any(dim(data)[2])<1)
    stop(paste("Object 'data' must have data"),"\n","")
  if(any(is.data.frame(data)==F))
    stop(paste("Object 'data' must be a data frame class","\n",""))
  if(any(length(link)<4))
    stop(paste("link must have four enlace functions", "\n",""))
  if(any(link[1]!='identity' && link[1]!='logit' && link[1]!='log'))
    stop(paste("link for mu must be in identity, logit, log", "\n",""))
  if(any(link[2]!='identity' && link[2]!='logit' && link[2]!='log'))
    stop(paste("link for sigma must be in identity, logit, log", "\n",""))
  if(any(link[3]!='identity' && link[3]!='logit'))
    stop(paste("link for p0 must be in identity, logit", "\n",""))
  if(any(link[4]!='identity' && link[4]!='logit'))
    stop(paste("link for p1 must be in identity, logit", "\n",""))
  if(any(length(as.character(attr(terms(formula.mu),'variable')))==2 && link[1]!='identity'))
    stop(paste("mu don't have covariables then link must be identity", "\n",""))
  if(any(family!='Original' && length(as.character(attr(terms(formula.mu),'variable')))>2 && link[1]!='logit'))
    stop(paste("If family is diferent a Original and mu have covariables then link must be logit", "\n",""))
  if(any(family=='Original'&& length(as.character(attr(terms(formula.mu),'variable')))>2 && link[1]!='log'))
    stop(paste("If family is Original and mu have covariables then link must be exp", "\n",""))
  if(any(length(as.character(attr(terms(formula.sigma),'variable')))==1 && link[2]!='identity'))
    stop(paste("sigma don't have covariables then link must be identity", "\n",""))
  if(any(family!='R-S' && length(as.character(attr(terms(formula.sigma),'variable')))>1 && link[2]!='log'))
    stop(paste("If family is diferent a R-S and sigma have covariables then link must be exp", "\n",""))
  if(any(family=='R-S'&& length(as.character(attr(terms(formula.sigma),'variable')))>1 && link[2]!='logit'))
    stop(paste("If family is R-S and sigma have covariables then link must be logit", "\n",""))
  if(any(length(as.character(attr(terms(formula.p0),'variable')))==1 && link[3]!='identity'))
    stop(paste("p0 don't have covariables then link must be identity", "\n",""))
  if(any(length(as.character(attr(terms(formula.p0),'variable')))>1 && link[3]!='logit'))
    stop(paste("p0 have covariables then link must be logit", "\n",""))
  if(any(length(as.character(attr(terms(formula.p1),'variable')))==1 && link[4]!='identity'))
    stop(paste("p1 don't have covariables then link must be identity", "\n",""))
  if(any(length(as.character(attr(terms(formula.p1),'variable')))>1 && link[4]!='logit'))
    stop(paste("p1 have covariables then link must be logit", "\n",""))
  if(any(optimizer!='nlminb' && optimizer!='optim'))
    stop(paste("optimizer should be 'nlminb' or 'optim'", "\n",""))
  
  if(any(length(as.character(attr(terms(formula.random),'variable')))!=2))
    stop(paste("formula.random don't have ~1 and/or group variable", "\n",""))
  if(any(grep("|",attr(terms(formula.random),'term.labels'))!=1))
    stop(paste("formula.random don't have symbol '|' for indicate group variable", "\n",""))
  if(any(class(pruning)!="logical"))
    stop(paste("pruning should be a logical object", "\n",""))
  if(any(n.points==0))
    stop(paste("n.points must have Greater than zero", "\n",""))
  
quad <- GHQp ::GHQ(n=n.points, ndim=2, pruning=pruning)

matri<-model.matrix.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma
                  ,formula.p0=formula.p0,formula.p1=formula.p1, data=data,
                  formula.random=formula.random)

opt<-fit.ZOIP2(matri=matri,link=link,family=family,optimizer=optimizer)

opt$par
theta0 <- c( opt$par, -1, -1) 
names(theta0) <- c(colnames(matri$mat.mu),colnames(matri$mat.sigma),colnames(matri$mat.p0),colnames(matri$mat.p1),"log(t1)","log(t2)")

nparm.mu <- ncol(matri$mat.mu)
nparm.sigma <- ncol(matri$mat.sigma)
nparm.p0 <- ncol(matri$mat.p0)
nparm.p1 <- ncol(matri$mat.p1)

lower.val=c(rep(ifelse((link[1]=='logit'|| link[1]=='log'),-Inf,1e-16),nparm.mu),rep(ifelse((link[2]=='logit' || link[2]=='log'),-Inf,1e-16),nparm.sigma),
            rep(ifelse(link[3]=='logit',-Inf,1e-16),nparm.p0),rep(ifelse(link[4]=='logit',-Inf,1e-16),nparm.p1))

upper.mu<-if(link[1]=='logit' || ((link[1]=='log' || link[1]=='identity') && family=='Original')){
  Inf
}else 0.999999999

upper.sigma<-if(link[2]=='logit' || link[2]=='log' || (link[2]=='identity' && family!='R-S')){
  Inf
}else 0.999999999

upper.val=c(rep(upper.mu,nparm.mu),rep(upper.sigma,nparm.sigma),
            rep(ifelse(link[3]=='logit',Inf,0.999999999),nparm.p0),rep(ifelse(link[4]=='logit',Inf,0.999999999),nparm.p1))

lower.val<-c(lower.val,-Inf,-Inf)
upper.val<-c(upper.val,Inf,Inf)



tf<-system.time( fit <- nlminb(theta0, llM, Y=matri$y, mat.mu=matri$mat.mu, mat.sigma=matri$mat.sigma,
              mat.p0=matri$mat.p0, mat.p1=matri$mat.p1,inter.ran=matri$inter.ran,
              quad=quad, link=link, family=family,
              control=list(eval.max=10000,iter.max=10000,trace=0,rel.tol=1e-5,x.tol=1e-3,xf.tol=1e-7),
              lower=lower.val,upper=upper.val)
)

elem.mu<-fit$par[1:nparm.mu]
elem.sigma<-fit$par[(nparm.mu+1):(nparm.mu+nparm.sigma)]
elem.p0<-fit$par[(nparm.mu+nparm.sigma+1):(nparm.mu+nparm.sigma+nparm.p0)]
elem.p1<-fit$par[(nparm.mu+nparm.sigma+nparm.p0+1):(nparm.mu+nparm.sigma+nparm.p0+nparm.p1)]

log_elem.ran.mu<-fit$par[(nparm.mu+nparm.sigma+nparm.p0+nparm.p1)+1]
log_elem.ran.sigma<-fit$par[(nparm.mu+nparm.sigma+nparm.p0+nparm.p1)+2]
elem.ran.mu<-exp(log_elem.ran.mu)
elem.ran.sigma<-exp(log_elem.ran.sigma)

elem.random<-cbind(rbind(elem.ran.mu,elem.ran.sigma),rbind(log_elem.ran.mu,log_elem.ran.sigma))
rownames(elem.random)<-c("Random Intercept mu","Random Intercept sigma")
colnames(elem.random)<-c("","log(.)")

elem.time<-tf[1]
names(elem.time)<-NULL

result<-list(Fixed_Parameters.mu=NULL,Fixed_Parameters.sigma=NULL,Fixed_Parameters.p0=NULL,Fixed_Parameters.p1=NULL
             ,Parameters.randoms=NULL,logverosimilitud=NULL,message=NULL,Time=NULL,num.iter=NULL)

result$Fixed_Parameters.mu <- elem.mu
result$Fixed_Parameters.sigma <- elem.sigma
result$Fixed_Parameters.p0 <- elem.p0
result$Fixed_Parameters.p1 <- elem.p1
result$Parameters.randoms <- elem.random
result$logverosimilitud <- fit$objective
result$Time <- elem.time
result$message <- fit$message
result$num.iter <- fit$iterations

return(result)
}
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
##Funciones en orden de ejecucion-------------------------------------------------------
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------

model.matrix.ZOIP <- function(formula.mu,formula.sigma,formula.p0,formula.p1, data=NULL
                              ,formula.random) {
  stopifnot (class(formula.mu) == 'formula')
  stopifnot (class(formula.sigma) == 'formula')
  stopifnot (class(formula.p0) == 'formula')
  stopifnot (class(formula.p1) == 'formula')
  stopifnot (class(formula.random) == 'formula')
  #stopifnot (class(formula.random.sigma) == 'formula')
  response <- all.vars(formula.mu)[1]
  formula.sigma <- as.formula(paste(response, paste(as.character(formula.sigma),
                                                    collapse='')))
  formula.p0 <- as.formula(paste(response, paste(as.character(formula.p0),
                                                 collapse='')))
  formula.p1 <- as.formula(paste(response, paste(as.character(formula.p1),
                                                 collapse='')))
  mat.mu <- model.matrix(formula.mu, data)
  mat.sigma <- model.matrix(formula.sigma, data)
  mat.p0 <- model.matrix(formula.p0, data)
  mat.p1 <- model.matrix(formula.p1, data)
  
  inter.ran_aux <- paste0("data$",all.vars(formula.random)[1])
  inter.ran<-as.factor(eval(parse(text=inter.ran_aux)))

  #inter.ran.sigma_aux <- paste0("data$",all.vars(formula.random.sigma)[1])
  #inter.ran.sigma<-as.factor(eval(parse(text=inter.ran.sigma_aux)))  
  
  y <- model.frame(formula.mu, data=data)[, 1]
  matri<-list(y=y,mat.mu=mat.mu, mat.sigma=mat.sigma,mat.p0=mat.p0,mat.p1=mat.p1
              ,inter.ran=inter.ran)
  return(matri)
}


##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
## log-likelihood function FiXED EFFECTS -----------------------------------------------
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------

fit.ZOIP2<-function(matri,link,family,optimizer){
  nparm.mu <- ncol(matri$mat.mu)
  nparm.sigma <- ncol(matri$mat.sigma)
  nparm.p0 <- ncol(matri$mat.p0)
  nparm.p1 <- ncol(matri$mat.p1)
  
  X.mu <- matri$mat.mu
  X.sigma <- matri$mat.sigma
  X.p0 <- matri$mat.p0
  X.p1 <- matri$mat.p1
  y <- matri$y
  
  
  val.inic<-rep(0.1,nparm.mu+nparm.sigma+nparm.p0+nparm.p1)
  
  lower.val=c(rep(ifelse((link[1]=='logit'|| link[1]=='log'),-Inf,1e-16),nparm.mu),rep(ifelse((link[2]=='logit' || link[2]=='log'),-Inf,1e-16),nparm.sigma),
              rep(ifelse(link[3]=='logit',-Inf,1e-16),nparm.p0),rep(ifelse(link[4]=='logit',-Inf,1e-16),nparm.p1))
  
  upper.mu<-if(link[1]=='logit' || ((link[1]=='log' || link[1]=='identity') && family=='Original')){
    Inf
  }else 0.999999999
  
  upper.sigma<-if(link[2]=='logit' || link[2]=='log' || (link[2]=='identity' && family!='R-S')){
    Inf
  }else 0.999999999
  
  upper.val=c(rep(upper.mu,nparm.mu),rep(upper.sigma,nparm.sigma),
              rep(ifelse(link[3]=='logit',Inf,0.999999999),nparm.p0),rep(ifelse(link[4]=='logit',Inf,0.999999999),nparm.p1))
  
  if (optimizer == 'nlminb') {
    opt <- nlminb(start=val.inic, objective=ll.ZOIP2,
                  y=y, X.mu=X.mu, X.sigma=X.sigma,X.p0=X.p0,X.p1=X.p1,
                  link=link,family=family,lower=lower.val,upper=upper.val)
    opt$objective <- -opt$objective
  }
  
  if (optimizer == 'optim') {
    
    opt <- optim(par=val.inic, fn=ll.ZOIP2,
                 y=y, X.mu=X.mu, X.sigma=X.sigma,X.p0=X.p0,X.p1=X.p1,
                 link=link,family=family,lower=lower.val,upper=upper.val)
    opt$objective <- -opt$value
  }
  
  return(opt)
  
}



ll.ZOIP2<-function(theta,y,X.mu,X.sigma,X.p0,X.p1,link,family){
  betas.mu <- matrix(theta[1:ncol(X.mu)], ncol=1)
  betas.sigma <- matrix(theta[seq(ncol(X.mu)+1,ncol(X.mu)+ncol(X.sigma))], ncol=1)
  betas.p0 <- matrix(theta[seq(ncol(X.mu)+ncol(X.sigma)+1,ncol(X.mu)+ncol(X.sigma)+ncol(X.p0))], ncol=1)
  betas.p1 <- matrix(theta[seq(ncol(X.mu)+ncol(X.sigma)+ncol(X.p0)+1,ncol(X.mu)+ncol(X.sigma)+ncol(X.p0)+ncol(X.p1))], ncol=1)
  
  if(link[1]=='identity'){
    mu<-X.mu %*% betas.mu
  }else if(link[1]=='logit'){
    mu<-1 / (1 + exp(- X.mu %*% betas.mu))
  }else if(link[1]=='log'){
    mu<-exp(X.mu %*% betas.mu)
  }
  
  if(link[2]=='identity'){
    sigma<-X.sigma %*% betas.sigma
  }else if(link[2]=='logit'){
    sigma<-1 / (1 + exp(- X.sigma %*% betas.sigma))
  }else if(link[2]=='log'){
    sigma<-exp(X.sigma %*% betas.sigma)
  }
  
  if(link[3]=='identity'){
    p0<-X.p0 %*% betas.p0
  }else if(link[3]=='logit'){
    nu<- exp(X.p0 %*% betas.p0)
    #p0<-1 / (1 + exp(- X.p0 %*% betas.p0))
  }
  
  if(link[4]=='identity'){
    p1<-X.p1 %*% betas.p1
  }else if(link[4]=='logit'){
    tau<-exp(X.p1 %*% betas.p1)
    #p1<-1 / (1 + exp(- X.p1 %*% betas.p1))
  }
  
  if(link[3]=='logit' && link[4]=='logit'){
    p0<-nu/(1+nu+tau)
    p1<-tau/(1+nu+tau)
  }else if(link[3]=='logit' && link[4]=='identity'){
    p0<-(nu-(p1*nu))/(nu+1)
  }else if(link[3]=='identity' && link[4]=='logit'){
    p1<-(tau-(p0*tau))/(tau+1)
  }
  
  ll<-sum(dZOIP(x=y,mu=mu,sigma=sigma,p0=p0,p1=p1,family=family,log=TRUE))
  -ll
}


##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------
## log-likelihood function MIXED EFFECTS -----------------------------------------------
##--------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------

llM <- function(theta, Y, mat.mu, mat.sigma, mat.p0, mat.p1, inter.ran, quad,link,family) {

  nbeta1 <- ncol(mat.mu) 
  nbeta2 <- ncol(mat.sigma)
  nbeta3 <- ncol(mat.p0)
  nbeta4 <- ncol(mat.p1)
  beta1 <- theta[1:nbeta1] 
  beta2 <- theta[(nbeta1+1):(nbeta1+nbeta2)]
  beta3 <- theta[(nbeta1+nbeta2+1):(nbeta1+nbeta2+nbeta3)]
  beta4 <- theta[(nbeta1+nbeta2+nbeta3+1):(nbeta1+nbeta2+nbeta3+nbeta4)]
  
  t1 <- exp(theta[nbeta1+nbeta2+nbeta3+nbeta4+1])
  t2 <- exp(theta[nbeta1+nbeta2+nbeta3+nbeta4+2])
  
  i <- levels(inter.ran) # despues cuando programemos q puede haber intercepto en una y en la otra no esto se debe cambiar
  
  ll <- lapply(i, llind, Y=Y, mat.mu=mat.mu, mat.sigma=mat.sigma, mat.p0=mat.p0, mat.p1=mat.p1,
               beta1=beta1,beta2=beta2,beta3,beta4=beta4, inter.ran=inter.ran,t1=t1, t2=t2,
               quad=quad,link=link,family=family) 
  
  -sum(log(unlist(ll)))
}


llind <- function(i, Y, mat.mu, mat.sigma, mat.p0, mat.p1,
                  beta1, beta2, beta3, beta4, inter.ran, t1,t2, quad, link, family) {
  
  nparm.mu <- ncol(mat.mu)
  nparm.sigma <- ncol(mat.sigma)
  nparm.p0 <- ncol(mat.p0)
  nparm.p1 <- ncol(mat.p1)
  
  y  <-  Y[inter.ran==i]
  X.mu <- mat.mu[inter.ran==i,]
  X.sigma <- mat.sigma[inter.ran==i,]
  X.p0 <- mat.p0[inter.ran==i,]
  X.p1 <- mat.p1[inter.ran==i,]

  
  lower.val=c(rep(ifelse((link[1]=='logit'|| link[1]=='log'),-Inf,1e-16),nparm.mu),rep(ifelse((link[2]=='logit' || link[2]=='log'),-Inf,1e-16),nparm.sigma),
              rep(ifelse(link[3]=='logit',-Inf,1e-16),nparm.p0),rep(ifelse(link[4]=='logit',-Inf,1e-16),nparm.p1))
  
  upper.mu<-if(link[1]=='logit' || ((link[1]=='log' || link[1]=='identity') && family=='Original')){
    Inf
  }else 0.999999999
  
  upper.sigma<-if(link[2]=='logit' || link[2]=='log' || (link[2]=='identity' && family!='R-S')){
    Inf
  }else 0.999999999
  
  upper.val=c(rep(upper.mu,nparm.mu),rep(upper.sigma,nparm.sigma),
              rep(ifelse(link[3]=='logit',Inf,0.999999999),nparm.p0),rep(ifelse(link[4]=='logit',Inf,0.999999999),nparm.p1))
  
  lower.val<-c(lower.val,-Inf,-Inf)
  upper.val<-c(upper.val,Inf,Inf)
    
  
  opt <- optim(par=c(0,0), fn=integrando,
               y=y, X.mu=X.mu, X.sigma=X.sigma, X.p0=X.p0, X.p1=X.p1,
               beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
               t1=t1, t2=t2, link=link, family=family,
               log=TRUE,
               hessian=TRUE, method="L-BFGS-B",
               control=list(fnscale=-1),lower=lower.val,upper=upper.val)
  
  x.hat <- opt$par
  Q   <- solve(-opt$hessian)
  Q12 <- chol(Q)
  Z   <- x.hat + sqrt(2) * t(Q12%*%t(quad$nodes))
  norma <- exp(-rowSums(quad$nodes^2))
  temp <- integrando(Z, y=y, X.mu=X.mu, X.sigma=X.sigma, X.p0=X.p0, X.p1=X.p1,
                     beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
                     t1=t1, t2=t2, link=link, family=family, log=FALSE)
  integral <- 2 * det(Q) * sum(quad$product * temp / norma)
  return(integral)
}


integrando <- function(u, y, X.mu,X.sigma, X.p0, X.p1, 
                       beta1, beta2, beta3, beta4, t1,t2, link, family,log=TRUE) {
  
  if(class(dim(u)) == "NULL"){u <- matrix(u,nrow=1,ncol=2)}
  
  ll <- apply(u,1,function(ui){
    
    if(link[1]=='identity'){
      mu<-X.mu %*% beta1 + ui[1]
    }else if(link[1]=='logit'){
      mu<-1 / (1 + exp(- X.mu %*% beta1 + ui[1]))
    }else if(link[1]=='log'){
      mu<-exp(X.mu %*% beta1 + ui[1])
    }
    
    if(link[2]=='identity'){
      sigma<-X.sigma %*% beta2 + ui[2]
    }else if(link[2]=='logit'){
      sigma<-1 / (1 + exp(- X.sigma %*% beta2 + ui[2]))
    }else if(link[2]=='log'){
      sigma<-exp(X.sigma %*% beta2 + ui[2])
    }
    
    if(link[3]=='identity'){
      p0<-as.matrix(X.p0) %*% beta3
    }else if(link[3]=='logit'){
      nu<- exp(X.p0 %*% beta3)
    }
    
    if(link[4]=='identity'){
      p1<-as.matrix(X.p1) %*% beta4
    }else if(link[4]=='logit'){
      tau<-exp(X.p1 %*% beta4)
    }
    
    if(link[3]=='logit' && link[4]=='logit'){
      p0<-nu/(1+nu+tau)
      p1<-tau/(1+nu+tau)
    }else if(link[3]=='logit' && link[4]=='identity'){
      p0<-(nu-(p1*nu))/(nu+1)
    }else if(link[3]=='identity' && link[4]=='logit'){
      p1<-(tau-(p0*tau))/(tau+1)
    }   
    
    
    # mu    <- 1 / (1 + exp(- c(x1%*%beta1) + ui)) #1 / (1 + exp(- c(x1%*%beta1) + rowSums(x1*ui) )) # 1 / (1 + exp(- c(x1%*%beta1) + ui))
    # 
    # sigma     <- 1 / (1 + exp(- x2 %*% beta2 ))
    # 
    # nu<- exp(x3 %*% beta3)
    # tau<-exp(x4 %*% beta4)
    # 
    # p0<-nu/(1+nu+tau)
    # p1<-tau/(1+nu+tau)
    
    temp1 <- sum( dZOIP(x=y, mu=mu, sigma=sigma,p0=p0,p1=p1, family=family, log=TRUE) )
    temp2<-dnorm(ui[1],mean=0,sd=t1,log=TRUE)
    temp3<-dnorm(ui[2],mean=0,sd=t2,log=TRUE)
    temp1 + temp2 + temp3 }) 
  if(log == FALSE) ll <- exp(ll)
  return(ll)
}




