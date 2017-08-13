library(ZOIP)
N<-120
n.points=11
pruning=T


Times <- c(1, 10, 20, 40) # cantidad de dias
Concentration.encoded <- rep(rep(-1:1, each=length(Times)), each=N/3) #3 por la cantidad de niveles en xij
subject <- rep(1:N, each=length(Times)) # numero de sujetos en la muestra repetidos tantas veces haya dias
Days <- rep(Times, times=N)
b0i <- rep(rnorm(n=N,sd=0.19), each=length(Times))
neta <- 1.61+b0i+0.91*log(Days)-0.85*(log(Days))^2+0.96*Concentration.encoded
# hay 60 datos correspondientes a 15 datos por cada periodo de tiempo
mu <- 1 / (1 + exp(-neta)) # aplicamos funcion link
mu[mu==1] <- 0.999 # Transforma datos inflados en unos
mu[mu==0] <- 0.001 # transforma datos inflados en ceros
require(VGAM)
Y <- rZOIP(n=length(mu), mu = mu, sigma = 0.2,p0=0.05,p1=0.05,family='R-S') # simulacion de datos con dispersion fija, gran detalle que la dispersion este fija (como hacer si no?)

base<-data.frame(Y,Days,Concentration.encoded,subject)

###ORGANIZACION DE DATOS # esta debe ser la primera funcion
formula <- Y ~ log(Days) + I(log(Days)^2) + Concentration.encoded # preguntar por que debe estar en clase AsIs
X1 <- model.matrix(formula)

datos<-list(Y=Y,X1=X1,subject=as.factor(subject)) # contruccion de la forma de la base simulada


quad <- GHQ(n=n.points, ndim=1, pruning=pruning)
theta0 <- fit.initial(datos) 
lower.val<-c(-Inf,-Inf,-Inf,-Inf,0.0001,0.0001,0.0001,-Inf)
upper.val<-c(Inf,Inf,Inf,Inf,0.9999,0.9999,0.9999,Inf)
fit <- nlminb(theta0, llM, Y=datos$Y, X1=datos$X1,subject=datos$subject, quad=quad,control=list(eval.max=10000,
                           iter.max=10000,trace=0),
              lower=lower.val,upper=upper.val)






#Funciones en orden de ejecucion-------------------------------------------------------------------------
#link=c('logit','identity','identity','identity')
#RM.ZOIP(formula.mu=formula,formula.sigma=~1,formula.p0=~1,formula.p1=~1,data=base,link=link,family='R-S')

fit.initial <- function(datos) {
  require(VGAM)
  aux <- nlminb(c(0.1,0.1,0.1,0.1,0.1,0.1,0.1), llF, Y=datos$Y, X1=datos$X1,lower=c(-Inf,-Inf,-Inf,-Inf,0.0001,0.0001,0.0001),upper=c(Inf,Inf,Inf,Inf,0.9999,0.9999,0.9999)) # ajustamos bajo el optimizador nlimnb la maxima verosimilitud de los paramestos, la funcion llf recive por aparte la y y las covariables, en nuestro caso previamente resivira una formula
  # debe tener 5 parametros ya que 4 estan asociados a los efectos fijos (como si no existieran los aleatorios) y uno a la estimacion de la dispersion, a pesar de que es fija.
  aux$par
  theta0 <- c( aux$par, -0.69) # define valor iniciales para los parametros alatorios Ya hay 5 +2=7 parametros incluyendo los alatorios
  names(theta0) <- c(colnames(datos$X1), "sigma","p0","p1","log(t1)") # coloca nombres 
  theta0 # resultado de la funcion 
}


llF <- function(theta, Y, X1) {
  nbeta1 <- ncol(X1) # numero de covariables
  beta1  <- theta[1:nbeta1] # definicion de thetas a hallar asociados con el numero de covariables
  sigma  <- theta[nbeta1+1] # esta con exp por el dominio o el enlace que debe tener la dispersion en la simplex
  p0<-theta[nbeta1+2]
  p1<-theta[nbeta1+3]
  
  mu     <- 1 / (1 + exp(- X1 %*% beta1 )) # transaformacion de mu para que quede con thetas a hallar

  -sum(dZOIP(x=Y, mu=mu, sigma=sigma,p0=p0,p1=p1, family='R-S', log=TRUE)) # la funcion de verosimilitud a maximizar
}

#-----------------
# log-likelihood function MIXED EFFECTS
llM <- function(theta, Y, X1, subject, quad) {
  N <- nlevels(subject) # el numero de sujetos que hay en la base que es 15 esto es necesario?
  nbeta1 <- ncol(X1) # numero de parametros a estimar sin tener encuenta la catidad de efectos aleatorios
  beta1  <- theta[1:nbeta1] #definicon de los thetas para efectos fijos
  sigma  <- theta[nbeta1+1]# definicion de theta para la dispersion, asoociada a la funcion link log
  p0<-theta[nbeta1+2]
  p1<-theta[nbeta1+3]
  
  t1     <- exp( theta[nbeta1+4] ) # definicion de theta para el intercepto aleatorio, asociada a la funcion link log como el dominio del parametro
  #t2     <- exp( theta[nbeta1+3] ) # definicion de theta para la pendiente aleatoria, asociada a la funcion link log como normal
  i <- 1:N # una secuencia de 1 hasta el numero de sujetos # Esto debe ser mas automatico
  ll <- lapply(i, llind, Y, X1, sigma,p0,p1, subject, beta1, t1, quad) # a cada i le aplicara la funcion llind, i varia de 1 a el numero de sujetos
  ## la funcion llind requiere como onjetos la Y, X1, sigma, subject, beta1, t1,t2, quad
  -sum(log(unlist(ll)))  # funcion a optimizar -sum del log (siempre es log ?) por que es normal? es la funcion log verosimilitud por eso log
}


llind <- function(i, Y, X1, sigma, p0, p1, subject, beta1, t1, quad) {
  y  <-  Y[subject==i] # las y del sujeto i
  x1 <- X1[subject==i,] # las covariables del sujeto i 
  opt <- optim(par=c(0), fn=integrando,
               y=y, x1=x1,
               beta1=beta1, sigma=sigma, p0=p0, p1=p1,
               t1=t1,
               log=TRUE,
               hessian=TRUE, method="BFGS", ## por que optim y no nlimnb?
               control=list(fnscale=-1),lower=lower.val,upper=upper.val) ## vamoa optimizar la funcion integrando, tien dos incognitas asociadas a los efectos aleatorios
  x.hat <- opt$par # entrega ls estimaciones donde ocurren los maximos de las x's
  Q   <- solve(-opt$hessian) # calcula la funcion inversa de la negativa de la matriz hessian 
  Q12 <- chol(Q) # calcula la descomposion de choleskey
  Z   <- x.hat + sqrt(2) * t(Q12%*%t(quad$nodes)) # calcula los z* centrados 
  norma <- exp(-rowSums(quad$nodes^2))
  temp <- integrando(Z, y=y, x1=x1, beta1=beta1, sigma=sigma, p0=p0, p1=p1,
                     t1=t1, log=FALSE)
  integral <- 2 * det(Q) * sum(quad$product * temp / norma) # creo que las 3 lineas anteriores es para calcular la aproximacion de hermite sobre la integral
}


integrando <- function(u,y,x1,beta1,sigma, p0, p1, t1,log=TRUE) {
  if(class(dim(u)) == "NULL"){u <- matrix(u,nrow=1,ncol=1)}
  ll <- apply(u,1,function(ui){
    mu    <- 1 / (1 + exp(- c(x1%*%beta1) + rowSums(x1*ui) )) #Calcula la mu para cada fila de observaciones del sujeto i, con su respectiva funcion enlace
    temp1 <- sum( dZOIP(x=y, mu=mu, sigma=sigma,p0=p0,p1=p1, family='R-S', log=TRUE) ) # colocamos los mu y sigma donde sigma es theta y mu es una combinacion lineal de thetas por las covariables
    temp2<-dnorm(ui,mean=0,sd=t1^2,log=TRUE)
    #temp2 <- d2norm(ui, l1=t1^2, l2=t2^2, rho=0, log=TRUE) # realizara d2norm por que hay variables aletorias,rho=0 por que los interceptos no estas correlacionados
    temp1 + temp2  } ) # suma las dos densidades
  if(log == FALSE) ll <- exp(ll) #optimiza si el log==False lo hace con exp que viene siendo sin log
  return(ll)
}
