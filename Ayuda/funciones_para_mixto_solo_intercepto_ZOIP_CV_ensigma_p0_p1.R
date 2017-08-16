library(ZOIP)
N<-60
n.points=15
pruning=T


Times <- c(2, 10, 20, 40) # cantidad de dias
Concentration.encoded <- rep(rep(-1:1, each=length(Times)), each=N/3) #3 por la cantidad de niveles en xij
subject <- rep(1:N, each=length(Times)) # numero de sujetos en la muestra repetidos tantas veces haya dias
Days <- rep(Times, times=N)
b0i <- rep(rnorm(n=N,sd=0.19), each=length(Times))
neta <- 1.61+b0i+0.91*log(Days)-0.85*(log(Days))^2+0.96*Concentration.encoded
neta2<-0.5-0.8*log(Days)
neta3<-0.1-1.4*log(Days)
neta4<-0.2-2.2*log(Days)
# hay 60 datos correspondientes a 15 datos por cada periodo de tiempo
mu <- 1 / (1 + exp(-neta)) # aplicamos funcion link
sigma <- 1 / (1 + exp(-neta2))
p0 <- 1 / (1 + exp(-neta3))
p1 <- 1 / (1 + exp(-neta4))

mu[mu==1] <- 0.999 # Transforma datos inflados en unos
mu[mu==0] <- 0.001 # transforma datos inflados en ceros

sigma[sigma==1] <- 0.999
sigma[sigma==0] <- 0.001

p0[p0==1] <- 0.999
p0[p0==0] <- 0.001

p1[p1==1] <- 0.999
p1[p1==0] <- 0.001


Y <- rZOIP(n=length(mu), mu = mu, sigma = sigma ,p0=p0,p1=p1,family='R-S') # simulacion de datos con dispersion fija, gran detalle que la dispersion este fija (como hacer si no?)

base<-data.frame(Y,Days,Concentration.encoded,subject)

###ORGANIZACION DE DATOS # esta debe ser la primera funcion
formula <- Y ~ log(Days) + I(log(Days)^2) + Concentration.encoded # preguntar por que debe estar en clase AsIs
formula.sigma<-~log(Days)
formula.p0<-~log(Days)
formula.p1<-~log(Days)

X1 <- model.matrix(formula)
X2<-model.matrix(formula.sigma)
X3<-model.matrix(formula.p0)
X4<-model.matrix(formula.p1)

datos<-list(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,subject=as.factor(subject)) # contruccion de la forma de la base simulada


quad <- GHQ(n=n.points, ndim=1, pruning=pruning)

nparm.mu <- ncol(datos$X1)
nparm.sigma <- ncol(datos$X2)
nparm.p0 <- ncol(datos$X3)
nparm.p1 <- ncol(datos$X4)

lower.val.initial<-c(rep(-Inf,nparm.mu+nparm.sigma+nparm.p0+nparm.p1))
upper.val.initial<-c(rep(Inf,nparm.mu+nparm.sigma+nparm.p0+nparm.p1))

theta0 <- fit.initial(datos) 

lower.val<-c(lower.val.initial,-Inf)
upper.val<-c(upper.val.initial,Inf)

fit <- nlminb(theta0, llM, Y=datos$Y, X1=datos$X1, X2=datos$X2, X3=datos$X3, X4=datos$X4
              ,subject=datos$subject, quad=quad,control=list(eval.max=10000,
                           iter.max=10000,trace=0),
              lower=lower.val,upper=upper.val)





#Funciones en orden de ejecucion-------------------------------------------------------------------------
#link=c('logit','identity','identity','identity')
#RM.ZOIP(formula.mu=formula,formula.sigma=~1,formula.p0=~1,formula.p1=~1,data=base,link=link,family='R-S')

fit.initial <- function(datos) {

  aux <- nlminb(rep(0.1,nparm.mu+nparm.sigma+nparm.p0+nparm.p1), llF, Y=datos$Y, X1=datos$X1, X2=datos$X2, X3=datos$X3, X4=datos$X4,
                lower=lower.val.initial,upper=upper.val.initial,control=list(eval.max=10000,
                                                                             iter.max=10000,trace=0)) # ajustamos bajo el optimizador nlimnb la maxima verosimilitud de los paramestos, la funcion llf recive por aparte la y y las covariables, en nuestro caso previamente resivira una formula
  # debe tener 5 parametros ya que 4 estan asociados a los efectos fijos (como si no existieran los aleatorios) y uno a la estimacion de la dispersion, a pesar de que es fija.
  aux$par
  theta0 <- c( aux$par, -0.69) # define valor iniciales para los parametros alatorios Ya hay 5 +2=7 parametros incluyendo los alatorios
  names(theta0) <- c(colnames(datos$X1), colnames(datos$X2),colnames(datos$X3),colnames(datos$X4),"log(t1)") # coloca nombres 
  theta0 # resultado de la funcion 
}


llF <- function(theta, Y, X1,X2,X3,X4) {
  # numero de covariables
  nbeta1 <- ncol(X1) 
  nbeta2 <- ncol(X2)
  nbeta3 <- ncol(X3)
  nbeta4 <- ncol(X4)
  beta1  <- theta[1:nbeta1] # definicion de thetas a hallar asociados con el numero de covariables
  beta2  <- theta[(nbeta1+1):(nbeta1+nbeta2)] # esta con exp por el dominio o el enlace que debe tener la dispersion en la simplex
  beta3<-theta[(nbeta1+nbeta2+1):(nbeta1+nbeta2+nbeta3)]
  beta4<-theta[(nbeta1+nbeta2+nbeta3+1):(nbeta1+nbeta2+nbeta3+nbeta4)]
  
  #meter funciones de enalace
  mu     <- 1 / (1 + exp(- X1 %*% beta1 )) # transaformacion de mu para que quede con thetas a hallar
  sigma     <- 1 / (1 + exp(- X2 %*% beta2 ))
  
  nu<- exp(X3 %*% beta3)
  tau<-exp(X4 %*% beta4)

  p0<-nu/(1+nu+tau)
  p1<-tau/(1+nu+tau)
  
  -sum(dZOIP(x=Y, mu=mu, sigma=sigma,p0=p0,p1=p1, family='R-S', log=TRUE)) # la funcion de verosimilitud a maximizar
}




#-----------------
# log-likelihood function MIXED EFFECTS
llM <- function(theta, Y, X1,X2, X3, X4, subject, quad) {
  N <- nlevels(subject)
  nbeta1 <- ncol(X1) 
  nbeta2 <- ncol(X2)
  nbeta3 <- ncol(X3)
  nbeta4 <- ncol(X4)
  beta1  <- theta[1:nbeta1] #definicon de los thetas para efectos fijos
  beta2  <- theta[(nbeta1+1):(nbeta1+nbeta2)] # esta con exp por el dominio o el enlace que debe tener la dispersion en la simplex
  beta3<-theta[(nbeta1+nbeta2+1):(nbeta1+nbeta2+nbeta3)]
  beta4<-theta[(nbeta1+nbeta2+nbeta3+1):(nbeta1+nbeta2+nbeta3+nbeta4)]
  
  t1     <- exp(theta[nbeta1+nbeta2+nbeta3+nbeta4+1] ) # definicion de theta para el intercepto aleatorio, asociada a la funcion link log como el dominio del parametro
  #t2     <- exp( theta[nbeta1+3] ) # definicion de theta para la pendiente aleatoria, asociada a la funcion link log como normal
  #i <- 1:N # una secuencia de 1 hasta el numero de sujetos # Esto debe ser mas automatico
  i <- as.numeric(levels(subject))
  ll <- lapply(i, llind, Y=Y, X1=X1,X2=X2, X3=X3, X4=X4,
               beta1=beta1,beta2=beta2,beta3,beta4=beta4,subject=subject,t1=t1,quad=quad) # a cada i le aplicara la funcion llind, i varia de 1 a el numero de sujetos
  ## la funcion llind requiere como onjetos la Y, X1, sigma, subject, beta1, t1,t2, quad
  -sum(log(unlist(ll)))  # funcion a optimizar -sum del log (siempre es log ?) por que es normal? es la funcion log verosimilitud por eso log
}


llind <- function(i, Y, X1,X2,X3,X4,
                  beta1,beta2,beta3,beta4, subject, t1, quad) {
  y  <-  Y[subject==i] # las y del sujeto i
  x1 <- X1[subject==i,] # las covariables del sujeto i 
  x2 <- X2[subject==i,]
  x3 <- X3[subject==i,]
  x4 <- X4[subject==i,]
  
  
  opt <- optim(par=c(0), fn=integrando,
               y=y, x1=x1,x2=x2,x3=x3,x4=x4,
               beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
               t1=t1,
               log=TRUE,
               hessian=TRUE, method="BFGS", ## por que optim y no nlimnb?
               control=list(fnscale=-1),lower=lower.val,upper=upper.val) ## vamoa optimizar la funcion integrando, tien dos incognitas asociadas a los efectos aleatorios
  x.hat <- opt$par # entrega ls estimaciones donde ocurren los maximos de las x's
  Q   <- solve(-opt$hessian) # calcula la funcion inversa de la negativa de la matriz hessian 
  Q12 <- chol(Q) # calcula la descomposion de choleskey
  Z   <- x.hat + sqrt(2) * t(Q12%*%t(quad$nodes)) # calcula los z* centrados 
  norma <- exp(-rowSums(quad$nodes^2))
  temp <- integrando(Z, y=y, x1=x1,x2=x2,x3=x3,x4=x4, beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
                     t1=t1, log=FALSE)
  integral <- 2 * det(Q) * sum(quad$product * temp / norma) # creo que las 3 lineas anteriores es para calcular la aproximacion de hermite sobre la integral
}


integrando <- function(u,y,x1,x2, x3, x4, beta1,beta2,beta3,beta4, t1,log=TRUE) {
  if(class(dim(u)) == "NULL"){u <- matrix(u,nrow=1,ncol=1)}
  ll <- apply(u,1,function(ui){
    mu    <- 1 / (1 + exp(- c(x1%*%beta1) + rowSums(x1*ui) )) #Calcula la mu para cada fila de observaciones del sujeto i, con su respectiva funcion enlace
    sigma     <- 1 / (1 + exp(- x2 %*% beta2 ))

    nu<- exp(x3 %*% beta3)
    tau<-exp(x4 %*% beta4)
    
    p0<-nu/(1+nu+tau)
    p1<-tau/(1+nu+tau)
    
    temp1 <- sum( dZOIP(x=y, mu=mu, sigma=sigma,p0=p0,p1=p1, family='R-S', log=TRUE) ) # colocamos los mu y sigma donde sigma es theta y mu es una combinacion lineal de thetas por las covariables
    temp2<-dnorm(ui,mean=0,sd=t1^2,log=TRUE)
     temp1 + temp2  } ) # suma las dos densidades
  if(log == FALSE) ll <- exp(ll) #optimiza si el log==False lo hace con exp que viene siendo sin log
  return(ll)
}
