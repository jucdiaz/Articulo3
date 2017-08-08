#===============================================================================
#================================ FUNCTIONS ====================================
#===============================================================================
#===============================================================================
# Simplex regression model
# Fixed effects and random intercept/slope uncorrelated
# Sigma parameter constant and unknown

# log-likelihood function FIXED EFFECTS
llF <- function(theta, Y, X1) {
  nbeta1 <- ncol(X1)
  beta1  <- theta[1:nbeta1]
  sigma  <- exp( theta[nbeta1+1] )
  mu     <- itrans.mu( X1 %*% beta1 )
  -sum( dsimplex(x=Y, mu=mu, dispersion=sigma, log=TRUE) )
}
#-------------------------------------------------------------------------------
# log-likelihood function MIXED EFFECTS
llM <- function(theta, Y, X1, subject, quad) {
  N <- nlevels(subject)
  nbeta1 <- ncol(X1)
  beta1  <- theta[1:nbeta1]
  sigma  <- exp( theta[nbeta1+1] )
  t1     <- exp( theta[nbeta1+2] )
  t2     <- exp( theta[nbeta1+3] )
  i <- 1:N
  ll <- lapply(i, llind, Y, X1, sigma, subject, beta1, t1, t2, quad)
  -sum(log(unlist(ll)))
}
#-------------------------------------------------------------------------------
# llfunc function
llind <- function(i, Y, X1, sigma, subject, beta1, t1, t2, quad) {
  y  <-  Y[subject==i]
  x1 <- X1[subject==i,]
  opt <- optim(par=c(0,0), fn=integrando,
               y=y, x1=x1,
               beta1=beta1, sigma=sigma,
               t1=t1, t2=t2,
               log=TRUE,
               hessian=TRUE, method="BFGS",
               control=list(fnscale=-1))
  x.hat <- opt$par
  Q   <- mysolve(-opt$hessian)
  Q12 <- chol(Q)
  Z   <- x.hat + sqrt(2) * t(Q12%*%t(quad$nodes))
  norma <- exp(-rowSums(quad$nodes^2))
  temp <- integrando(Z, y=y, x1=x1, beta1=beta1, sigma=sigma,
                     t1=t1, t2=t2, log=FALSE)
  integral <- 2 * det(Q) * sum(quad$product * temp / norma)
}
#-------------------------------------------------------------------------------
# Integrand function
integrando <- function(u,y,x1,beta1,sigma,t1,t2,log=TRUE) {
  if(class(dim(u)) == "NULL"){u <- matrix(u,nrow=1,ncol=2)}
  ll <- apply(u,1,function(ui){
    mu    <- itrans.mu( c(x1%*%beta1) + rowSums(x1*ui) )
    temp1 <- sum( dsimplex(x=y, mu=mu, dispersion=sigma, log=TRUE) )
    temp2 <- d2norm(ui, l1=t1^2, l2=t2^2, rho=0, log=TRUE)
    temp1 + temp2  } )
  if(log == FALSE) ll <- exp(ll)
  return(ll)
}
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
}
#-------------------------------------------------------------------------------
# Inverse of 2x2 matrix
mysolve <- function(M) {
  determ <- M[1]*M[4] - M[2]^2
  temp.M1 <- M[1]
  M[1] <- -M[4]
  M[4] <- -temp.M1
  -M / determ
}
#-------------------------------------------------------------------------------
# Functions to transform rho
transf.rho     <- function(rho) log((rho+1)/(1-rho))
inv.transf.rho <- function(x)   -1+2*(exp(x)/(exp(x)+1))
curve(transf.rho(x), -1, 1)
curve(inv.transf.rho(x), -10, 10)
#===============================================================================
#======================== END OF FUNCTIONS =====================================
#===============================================================================
#--------------------------------------------------------------------------
itrans.mu <- function(x) atan(x)/pi + 0.5
#--------------------------------------------------------------------------
# Function to create the dataset
#--------------------------------------------------------------------------
simul_data <- function(N) {
  Times <- c(1, 10, 20, 40)
  Concentration.encoded <- rep(rep(-1:1, each=length(Times)), each=N/3)
  subject <- rep(1:N, each=length(Times))
  Days <- rep(Times, times=N)
  b0i <- rep(rnorm(n=N,sd=0.19), each=length(Times))
  b1i <- rep(rnorm(n=N,sd=0.18), each=length(Times))
  neta <- 1.61+b0i+0.91*log(Days)-0.85*(log(Days))^2+(0.96+b1i)*Concentration.encoded
  mu <- itrans.mu(neta)
  mu[mu==1] <- 0.999
  mu[mu==0] <- 0.001
  require(VGAM)
  Y <- rsimplex(n=length(mu), mu = mu, dispersion = 3.78)
  formula <- Y ~ log(Days) + I(log(Days)^2) + Concentration.encoded
  X1 <- model.matrix(formula)
  zero.obs <- Y %in% c(0,1)
  X1 <- X1[!zero.obs,]
  Y  <-  Y[!zero.obs]
  subject <- subject[!zero.obs]
  
  list(Y=Y,X1=X1,subject=as.factor(subject))
}
#--------------------------------------------------------------------------
# Fitting with fixed effects only two obtain start values
#--------------------------------------------------------------------------
fit.initial <- function(datos) {
  require(VGAM)
  aux <- nlminb(c(3,0,-0.3,0.5,1), llF, Y=datos$Y, X1=datos$X1,
                control=list(eval.max=1000, iter.max=1000))
  aux$par
  theta0 <- c( aux$par, -0.69, -0.69)
  names(theta0) <- c(colnames(datos$X1), "log(sigma)", "log(t1)", "log(t2)")
  theta0
}
#--------------------------------------------------------------------------
# Fitting with random effects
#--------------------------------------------------------------------------
fit.random <- function(N, n.points, pruning, trace=1) {
  quad <- GHQ(n=n.points, ndim=2, pruning=pruning)
  datos <- simul_data(N)
  theta0 <- fit.initial(datos)
  t <- system.time( fit <- nlminb(theta0, llM, Y=datos$Y, X1=datos$X1,
                                  subject=datos$subject, quad=quad,
                                  control=list(eval.max=10000,
                                               iter.max=10000,trace=trace) ) )
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







