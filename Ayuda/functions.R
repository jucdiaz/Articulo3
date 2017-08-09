#===============================================================================
#================================ FUNCTIONS ====================================
#===============================================================================
#===============================================================================
# Simplex regression model
# Fixed effects and random intercept/slope uncorrelated
# Sigma parameter constant and unknown

# log-likelihood function FIXED EFFECTS
llF <- function(theta, Y, X1) {
  nbeta1 <- ncol(X1) # numero de covariables
  beta1  <- theta[1:nbeta1] # definicion de thetas a hallar asociados con el numero de covariables
  sigma  <- exp( theta[nbeta1+1] ) # esta con exp por el dominio o el enlace que debe tener la dispersion en la simplex
  mu     <- itrans.mu( X1 %*% beta1 ) # transaformacion de mu para que quede con thetas a hallar
  -sum( dsimplex(x=Y, mu=mu, dispersion=sigma, log=TRUE) ) # la funcion de verosimilitud a maximizar
}
#-------------------------------------------------------------------------------
# log-likelihood function MIXED EFFECTS
llM <- function(theta, Y, X1, subject, quad) {
  N <- nlevels(subject) # el numero de sujetos que hay en la base que es 15 esto es necesario?
  nbeta1 <- ncol(X1) # numero de parametros a estimar sin tener encuenta la catidad de efectos aleatorios
  beta1  <- theta[1:nbeta1] #definicon de los thetas para efectos fijos
  sigma  <- exp( theta[nbeta1+1] ) # definicion de theta para la dispersion, asoociada a la funcion link log
  t1     <- exp( theta[nbeta1+2] ) # definicion de theta para el intercepto aleatorio, asociada a la funcion link log como normal
  t2     <- exp( theta[nbeta1+3] ) # definicion de theta para la pendiente aleatoria, asociada a la funcion link log como normal
  i <- 1:N # una secuencia de 1 hasta el numero de sujetos
  ll <- lapply(i, llind, Y, X1, sigma, subject, beta1, t1, t2, quad) # a cada i le aplicara la funcion llind, i varia de 1 a el numero de sujetos
  ## la funcion llind requiere como onjetos la Y, X1, sigma, subject, beta1, t1,t2, quad
  -sum(log(unlist(ll)))  # funcion a optimizar -sum del log (siempre es log ?) por que es normal?
}
#-------------------------------------------------------------------------------
# llfunc function
llind <- function(i, Y, X1, sigma, subject, beta1, t1, t2, quad) {
  y  <-  Y[subject==i] # las y del sujeto i
  x1 <- X1[subject==i,] # las covariables del sujeto i 
  opt <- optim(par=c(0,0), fn=integrando,
               y=y, x1=x1,
               beta1=beta1, sigma=sigma,
               t1=t1, t2=t2,
               log=TRUE,
               hessian=TRUE, method="BFGS", ## por que optim y no nlimnb?
               control=list(fnscale=-1)) ## vamoa optimizar la funcion integrando, tien dos incognitas asociadas a los efectos aleatorios
  x.hat <- opt$par # entrega ls estimaciones donde ocurren los maximos de las x's
  Q   <- mysolve(-opt$hessian) # calcula la funcion inversa de la negativa de la matriz hessian 
  Q12 <- chol(Q) # calcula la descomposion de choleskey
  Z   <- x.hat + sqrt(2) * t(Q12%*%t(quad$nodes)) # calcula los z* centrados 
  norma <- exp(-rowSums(quad$nodes^2))
  temp <- integrando(Z, y=y, x1=x1, beta1=beta1, sigma=sigma,
                     t1=t1, t2=t2, log=FALSE)
  integral <- 2 * det(Q) * sum(quad$product * temp / norma) # creo que las 3 lineas anteriores es para calcular la aproximacion de hermite sobre la integral
}
#-------------------------------------------------------------------------------
# Integrand function
integrando <- function(u,y,x1,beta1,sigma,t1,t2,log=TRUE) {
  if(class(dim(u)) == "NULL"){u <- matrix(u,nrow=1,ncol=2)}
  ll <- apply(u,1,function(ui){
    mu    <- itrans.mu( c(x1%*%beta1) + rowSums(x1*ui) ) #Calcula la mu para cada fila de observaciones del sujeto i, con su respectiva funcion enlace
    temp1 <- sum( dsimplex(x=y, mu=mu, dispersion=sigma, log=TRUE) ) # colocamos los mu y sigma donde sigma es theta y mu es una combinacion lineal de thetas por las covariables
    temp2 <- d2norm(ui, l1=t1^2, l2=t2^2, rho=0, log=TRUE) # realizara d2norm por que hay variables aletorias,rho=0 por que los interceptos no estas correlacionados
    temp1 + temp2  } ) # suma las dos densidades
  if(log == FALSE) ll <- exp(ll) #optimiza si el log==False lo hace con exp que viene siendo sin log
  return(ll)
}# esta es como la funcion g del texto
#-------------------------------------------------------------------------------
# Our dmvnorm function
d2norm <- function(x,l1,l2,rho,log=FALSE){ # l1 y l2 are variances
  if(class(dim(x)) == "NULL"){x <- matrix(x,1,2)}
  dens <- apply(x,1,function(xi){
    t1 <- 1/(2*pi*sqrt(l1*l2*(1-rho^2)))
    t2 <- 1/(2*(1-rho^2))
    t3 <- xi[1]/sqrt(l1)
    t4 <- xi[2]/sqrt(l2)
    t1*exp(-t2*(t3^2 + t4^2 -2*rho*t3*t4)) } )
  if (log) dens <- log(dens)
  dens
}# Creoq ue esta funcion es la normal mutivariada para dos, preguntar por que la personalizo y no utilizo la 
#funcion dmvnorm, creo que esta funcion cambiaria a una normal normal i.e un dnorm
#-------------------------------------------------------------------------------
# Inverse of 2x2 matrix
mysolve <- function(M) {
  determ <- M[1]*M[4] - M[2]^2 # calcula el determinante
  temp.M1 <- M[1]
  M[1] <- -M[4]
  M[4] <- -temp.M1
  -M / determ
}# este calculo de inversa si es necesario para cuando solo hay intercepto, por que no utilizo la funcion solve?
#-------------------------------------------------------------------------------
# Functions to transform rho
transf.rho     <- function(rho) log((rho+1)/(1-rho))
inv.transf.rho <- function(x)   -1+2*(exp(x)/(exp(x)+1))
curve(transf.rho(x), -1, 1)
curve(inv.transf.rho(x), -10, 10) # lo de las transaformaciones rho no le encuentro sentido dentro del documento
#===============================================================================
#======================== END OF FUNCTIONS =====================================
#===============================================================================
#--------------------------------------------------------------------------
itrans.mu <- function(x) atan(x)/pi + 0.5 # esto no lo entiendo para que es ? creo que es el link de cauchy
#--------------------------------------------------------------------------
# Function to create the dataset
#--------------------------------------------------------------------------
simul_data <- function(N) {
  Times <- c(1, 10, 20, 40) # cantidad de dias
  Concentration.encoded <- rep(rep(-1:1, each=length(Times)), each=N/3) #3 por la cantidad de niveles en xij
  subject <- rep(1:N, each=length(Times)) # numero de sujetos en la muestra repetidos tantas veces haya dias
  Days <- rep(Times, times=N)
  b0i <- rep(rnorm(n=N,sd=0.19), each=length(Times))
  b1i <- rep(rnorm(n=N,sd=0.18), each=length(Times))
  neta <- 1.61+b0i+0.91*log(Days)-0.85*(log(Days))^2+(0.96+b1i)*Concentration.encoded # hay 60 datos correspondientes a 15 datos por cada periodo de tiempo
  mu <- itrans.mu(neta) # aplicamos funcion link
  mu[mu==1] <- 0.999 # Transforma datos inflados en unos
  mu[mu==0] <- 0.001 # transforma datos inflados en ceros
  require(VGAM)
  Y <- rsimplex(n=length(mu), mu = mu, dispersion = 3.78) # simulacion de datos con dispersion fija, gran detalle que la dispersion este fija (como hacer si no?)
  formula <- Y ~ log(Days) + I(log(Days)^2) + Concentration.encoded
  X1 <- model.matrix(formula)
  zero.obs <- Y %in% c(0,1) # verifica que no haya cero o unos en y
  X1 <- X1[!zero.obs,] # coge todas las covariables que no son unos ni ceros
  Y  <-  Y[!zero.obs] #coge todas la variables respuestas que no son unos ni ceros
  subject <- subject[!zero.obs] # coge todos las observaciones que no sean unos ni ceros
  
  list(Y=Y,X1=X1,subject=as.factor(subject)) # contruccion de la forma de la base simulada
}
#--------------------------------------------------------------------------
# Fitting with fixed effects only two obtain start values
#--------------------------------------------------------------------------
fit.initial <- function(datos) {
  require(VGAM)
  aux <- nlminb(c(3,0,-0.3,0.5,1), llF, Y=datos$Y, X1=datos$X1,
                control=list(eval.max=1000, iter.max=1000)) # ajustamos bajo el optimizador nlimnb la maxima verosimilitud de los paramestos, la funcion llf recive por aparte la y y las covariables, en nuestro caso previamente resivira una formula
  # debe tener 5 parametros ya que 4 estan asociados a los efectos fijos (como si no existieran los aleatorios) y uno a la estimacion de la dispersion, a pesar de que es fija.
  aux$par
  theta0 <- c( aux$par, -0.69, -0.69) # define valor iniciales para los parametros alatorios Ya hay 5 +2=7 parametros incluyendo los alatorios
  names(theta0) <- c(colnames(datos$X1), "log(sigma)", "log(t1)", "log(t2)") # coloca nombres 
  theta0 # resultado de la funcion 
}
#--------------------------------------------------------------------------
# Fitting with random effects
#--------------------------------------------------------------------------
fit.random <- function(N, n.points, pruning, trace=1) {
  quad <- GHQ(n=n.points, ndim=2, pruning=pruning) # cacula los puntos de gauss hermite con o sin pruning depende del usuario ndim2 por la catidad de efectos alatorios que hay
  datos <- simul_data(N) # la base de entrada debe tener la estructura de lista con un objeto que tenga la Y un vector, una matriz de covariables dado por model.matrix, la identificacion de los sujetos o observaciones (para nuestro caso no creo que sea necesaria) 
  theta0 <- fit.initial(datos) # ajuste de los efectos fijos aca utilizariamos la funcion RMZOIP, se tomaran como valores inciales
  t <- system.time( fit <- nlminb(theta0, llM, Y=datos$Y, X1=datos$X1,
                                  subject=datos$subject, quad=quad,
                                  control=list(eval.max=10000,
                                               iter.max=10000,trace=trace) ) )# los valores inciales son theta cero que incluye por aparte y al final la estimacion de los efectos alatorios
  ##Utilizo de nuevo la optimizacion de funcion para maximizar la funcion LLM
  ## con el vector de datos de Y y la matriz de datos x1 dado por model.matrix, necesita los sujetos
  ##quad necesita todos los objetos de la lista de la cuadratura de gauss hermite
  ##poner cuidado en el trace que siempre esta en cero y el numero de iteraciones y evaluacion limite son 10 mil
  
  c(time=t[1], param=fit$par)
}
#--------------------------------------------------------------------------
require(GHQp)
#--------------------------------------------------------------------------
simulation <- function(Nrep, N, n.points, pruning) {
  result <- NULL
  for (i in 1:Nrep) { result <- rbind(result, fit.random(N, n.points, pruning, trace=0)) }
  write(t(result), ncolumns=8,
        file=paste("N=",N,"_n.points=",n.points,"_pruning=",pruning,".txt", sep=""))
}
#--------------------------------------------------------------------------
simulation2 <- function(Nrep, N, n.points, pruning) {
  name=paste("N=",N,"_n.points=",n.points,"_pruning=",pruning,".txt", sep="")
  result <- NULL
  for (i in 1:Nrep) { result <- rbind(result, fit.random(N, n.points, pruning, trace=0)) }
  write( t(result ), ncolumns=8, file=name, append=TRUE)
}
#--------------------------------------------------------------------------
mrd <- function(dataset) {
  theta <- c(1.61, 0.91, 0.50, 0.96, log(3.78), log(0.19), log(0.18))
  mrd.i <- colMeans(t(abs(t(dataset[,-1])-theta)/abs(theta)))
  norm.theta <- sqrt(sum(theta^2))
  mrd.theta <- mean(apply(t(t(dataset[,-1])-theta), 1,
                          function(x) sqrt(sum(x^2))) / norm.theta)
  c(mrd.i, theta=mrd.theta)
}
#--------------------------------------------------------------------------
# Estas dos funciones son iguales, hacen lo mismo.
# La primera recibe tres elementos y la segunda solo un vector de 3 componentes
final.result1 <- function (N,n.points,pruning) {
  file <- paste("N=",N,"_n.points=",n.points,"_pruning=",pruning,".txt", sep="")
  dataset <- read.table(file)
  names(dataset) <- c('time','b0','blogdays','blog2days','bxij',
                      'logsigma','logsigmaint','logsigmapend')
  mean.time <- mean(dataset[,1])
  round(c(mrd(dataset), mean.time=mean.time), digits=2)
}
#------
# Esta función recibe el N, n.points y pruning y calcula mrd y mean para los parámetros
final.result2 <- function (x) {
  N <- x[1] ; n.points <- x[2]; pruning <- as.logical(x[3])
  file <- paste("N=",N,"_n.points=",n.points,"_pruning=",pruning,".txt", sep="")
  dataset <- read.table(file)
  names(dataset) <- c('time','b0','blogdays','blog2days','bxij',
                      'logsigma','logsigmaint','logsigmapend')
  mean.time <- mean(dataset[,1])
  round(c(mrd(dataset), mean.time=mean.time), digits=2)
}
#--------------------------------------------------------------------------
verificador_n <- function(x) {
  N <- x[1]; n.points <- x[2]; pruning <- as.logical(x[3])
  name=paste("N=",N,"_n.points=",n.points,"_pruning=",pruning,".txt", sep="")
  datos <- read.table(name)
  dimension <- dim(datos)
  names(dimension) <- c('NREP','nvar')
  dimension
}
#--------------------------------------------------------------------------
tabla.latex <- function(N) {
  casos <- data.frame(N=rep(c(15,21,30,60,90,120), each=6), 
                      n.points=rep(c(3,5,7,11,15,21), times=6),
                      pruning=rep(c(FALSE,TRUE), each=3, times=3))
  x <- apply(casos[casos$N==N,], 1, final.result2)
  names <- c("$\\beta_0$ &","$\\beta_1$ &","$\\beta_2$ &","$\\beta_3$ &",
             "$\\log(\\sigma)$ &","$\\log(\\tau_1)$ &","$\\log(\\tau_2)$ &",
             "Mean relative distance &", "Mean time (sec) &")
  y <- array(NA, dim=c(9,1))
  for (i in 1:9)
    y[i, ] <- paste(names[i], paste(x[i,], collapse=" & "), "\\\\", ifelse(i>=7," \\hline", ""))
  write(t(y), ncolumns=1, file=paste("tabla_latex_N=",N,".txt",sep=""))
}

# pruebita de tabla final
tabla.latex <- function(N) {
  NN=c(15,21,30,60,90,120)
  n.points=c(3,11,21)
  pruning=c(FALSE,TRUE)
  casos <- expand.grid(N=NN, n.points=n.points, pruning=pruning)
  
  x <- apply(casos[casos$N==N,], 1, final.result2)
  names <- c("$\\beta_0$ &","$\\beta_1$ &","$\\beta_2$ &","$\\beta_3$ &",
             "$\\log(\\sigma)$ &","$\\log(\\tau_1)$ &","$\\log(\\tau_2)$ &",
             "Mean relative distance &", "Mean time (sec) &")
  y <- array(NA, dim=c(9,1))
  for (i in 1:9)
    y[i, ] <- paste(names[i], paste(x[i,], collapse=" & "), "\\\\", ifelse(i>=7," \\hline", ""))
  write(t(y), ncolumns=1, file=paste("tabla_latex_N=",N,".txt",sep=""))
}
#--------------------------------------------------------------------------


# Funciones para crear la figura de tiempo vs n para cada N ---------------

tiempos <- function(x) {
  N <- x[1] ; n.points <- x[2]; pruning <- as.logical(x[3])
  file <- paste("N=",N,"_n.points=",n.points,"_pruning=",pruning,".txt", sep="")
  dataset <- read.table(file)
  mean.time <- mean(dataset[,1])
  data.frame(mean.time=mean.time, n.points, pruning, N)
}

tiempos.all <- function() {
  Times <- NULL
  NN <- c(15, 21, 30, 60, 90, 120)
  nn <- c(3, 5, 7, 11, 15, 21)
  pruning <- 0:1
  for (i in NN)
    for (j in nn)
      for (k in pruning)
        Times <- rbind(Times, tiempos(c(i, j, k)))
  Times
}







