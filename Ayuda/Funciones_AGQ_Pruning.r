# Installing the appropriate package
list.of.packages <- c("mvtnorm", "statmod", "scatterplot3d", "GHQp")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# With this instruction we can obtain the quadrature points and weights with n=5
require(statmod)
quad <- gauss.quad(n=5, kind="hermite")
quad
#-------------------------------------------------------------------------------
# Example 1
#-------------------------------------------------------------------------------
# Defining g1(x)
g1     <- function(x) exp(-(x-1)^2)

# To obtain figure (1)
curve(g1, -4, 6, ylim=c(0,1), ylab=expression(g[1](x)), las=1)

# Adding the quadrature points
points(x=quad$nodes, y=rep(0,5), pch=19, cex=1.2)
legend('topright', bty='n', legend='Quadrature points', pch=19, pt.cex=1.2)

# Approximating the integral of g1(x) by expression (1)
sum(quad$weights * g1(quad$nodes) * exp(quad$nodes^2))
# True value with integrate function
integrate(f=g1, lower=-Inf, upper=Inf) 

#-------------------------------------------------------------------------------
# Example 2
#-------------------------------------------------------------------------------
# Defining g2(x) and log(g2(x))
g2     <- function(x) exp(-5*(x-3)^2)
log.g2 <- function(x) -5*(x-3)^2

# To obtain figure (2)
par(mfrow=c(1,2))
curve(g2, -3, 5, ylim=c(0,1), ylab=expression(g[2](x)), las=1)

# Adding the quadrature points
points(x=quad$nodes, y=rep(0,5), pch=19, cex=1.2)
legend('topright', bty='n', legend=expression(p[i]), pch=19, pt.cex=1.2)

# Approximating the integral of g2(x) by expression (1)
sum(quad$weights * g2(quad$nodes) * exp(quad$nodes^2))
# True value with integrate function
integrate(f=g2, lower=-Inf, upper=Inf)

# Using the Adaptive version
      opt <- optim(fn=log.g2,par=0,control=list(fnscale=-1),hessian=T,method='CG')
    x.hat <- opt$par
sigma.hat <- as.numeric(sqrt(-1/opt$hessian))
       xi <- sqrt(2)*sigma.hat*quad$nodes + x.hat
       wi <- quad$weights

# Plotting g2(x) with new quadrature points       
curve(g2, -3, 5, ylim=c(0,1), ylab=expression(g[2](x)), las=1)
points(x=xi, y=rep(0,5), pch='*', cex=2)
legend('topright', bty='n', legend=expression(p[i]^symbol("*")), pch='*', pt.cex=2)

# Approximating the integral of g2(x) by expression (2)
sqrt(2) * sigma.hat * sum( wi * g2(xi) * exp(quad$nodes^2) )
# True value with integrate function
integrate(f=g2, lower=-Inf, upper=Inf) # valor de la integral con integrate

#-------------------------------------------------------------------------------
# Example 3
#-------------------------------------------------------------------------------
require(mvtnorm)
# Defining the mean vector an covariance matrix
mu    <- c(-2,-2)
sigma <- matrix(c(1,-0.5,-0.5,2),2,2)

# Defining g(x1,x2) and log(g(x1,x2))
    g <- function(x,mu,sigma) dmvnorm(x=c(x[1],x[2]), mean=mu, sigma=sigma)
log.g <- function(x,mu,sigma) dmvnorm(x=c(x[1],x[2]), mean=mu, sigma=sigma, log=TRUE)

# To plot the figure (3)
x1 <- seq(mu[1]-4,mu[1]+4,length.out=50)
x2 <- seq(mu[2]-4,mu[2]+4,length.out=50)
z  <- matrix(NA,ncol=length(x2),nrow=length(x1))
for(j in 1:length(x1)) for(k in 1:length(x2)) z[j,k] <- dmvnorm(x=c(x1[j],x2[k]), mean=mu, sigma=sigma)
par(mfrow=c(1,2))
plot(x = 0, y = 0,type = "n", xlim = c(-6.2,2.2), ylim = c(-6,2), xlab = expression(x[1]), ylab = expression(x[2]), las=1)
contour(x1, x2, z, col = 'gray', lty = "solid", add = TRUE, nlevels=15)
text(1, -5, expression(g(x[1],x[2])), col='gray')

# Obtaining the quadrature points
require(statmod)
  nq <- 5
quad <- gauss.quad(n=nq, kind="hermite")
# The next matrices contain the cartesian products with the nodes and weights
  CPnodes <- as.matrix(expand.grid(quad$nodes,quad$nodes))
CPweights <- as.matrix(expand.grid(quad$weights,quad$weights))

# Adding the quadrature points
points(CPnodes, pch=19, cex=1.2)
legend('bottomleft', bty='0', legend=expression(z[i]), pch=19, pt.cex=1.2)

# To obtain the mode of g(x1,x2) 
opt  <- optim(fn=log.g, par=c(0,0), mu=mu, sigma=sigma, control=list(fnscale=-1), hessian=T)
moda <- matrix(opt$par, byrow=T, nrow=nrow(CPnodes), ncol=2)
Q <- solve(-opt$hessian)
B <- chol(Q)
Z <- moda + sqrt(2) * t(B%*%t(CPnodes))

# Adding the new quadrature points that were transformed
plot(x = 0, y = 0,type = "n", xlim = c(-6.2,2.2), ylim = c(-6,2), xlab = expression(x[1]), ylab = expression(x[2]), las=1)
contour(x1, x2, z, col = 'gray', lty = "solid", add = TRUE, nlevels=15)
text(1,-5,expression(g(x[1],x[2])),col='gray')
points(Z, pch='*', cex=2)
legend('bottomleft', bty='o', legend=expression(z[i]^symbol("*")), pch='*', pt.cex=2)

# Obtaining the integral with expression (3), Multivariate Gaussian Quadrature
sum( apply(X=CPnodes,MARGIN=1,FUN=g,mu=mu,sigma=sigma) * exp(rowSums(CPnodes^2)) * apply(CPweights,1,prod) )

# Obtaining the integral with expression (4), Multivariate Adaptive Gaussian Quadrature
sqrt(det(Q)) * 2 *  sum( apply(X=Z,MARGIN=1,FUN=g,mu=mu,sigma=sigma) * exp(rowSums(CPnodes^2)) * apply(CPweights,1,prod) )

#-------------------------------------------------------------------------------
# To obtain figure (4)
#-------------------------------------------------------------------------------
# It is necessary to load the function GHQ from the another supplement
require(GHQp)
par(mfrow=c(2,2))
plot(GHQ(15,2,FALSE)$nodes,pch=20,xlab='',ylab='',main='Without pruning, n=15 and q=2',las=1)
plot(GHQ(15,2,TRUE)$nodes, pch=20,xlab='',ylab='',main='With pruning, n=15 and q=2',las=1)
require(scatterplot3d)
datos <- GHQ(15,3,FALSE)$nodes
scatterplot3d(datos, type="p", highlight.3d=TRUE,angle=55, scale.y=0.7, pch=16,
              main='Without pruning, n=15 and q=3', cex.symbols=0.4,xlab='',ylab='',zlab='')
datos <- GHQ(15,3,TRUE)$nodes
scatterplot3d(datos, type="p", highlight.3d=TRUE,angle=55, scale.y=0.7, pch=16,
              main='With pruning, n=15 and q=3', cex.symbols=0.4,xlab='',ylab='',zlab='')
#-------------------------------------------------------------------------------
# Application

#-----------------------------------------------------------------------
# Here start the application
#-----------------------------------------------------------------------

link <- "https://dl.dropboxusercontent.com/u/9602072/eye.txt"
datos <- read.table(file=link, sep="", header=T)
head(datos)

# Creando el gráfico de evolución
require(lattice)

attach(datos)
concentration <- as.factor( paste('Concentration ',datos$Concentration,'%',sep='') )

xyplot(Proportion ~ Days | concentration, data=datos, as.table=TRUE,
       xlab='Days', ylab=expression("Proportion of gas remained ( " * y[ij] * " )"),
       groups=Id, type='b', col=1, pch=19, cex=0.3,
       par.settings = list(strip.background=list(col="white")))

#===============================================================================
#=============================== FUNCTIONS =====================================
#===============================================================================
# Simplex regression model
# Fixed effects and random intercept/slope uncorrelated
# Sigma parameter constant and unknown

# log-likelihood function FIXED EFFECTS
# log-likelihood function FIXED EFFECTS
llF <- function(theta, Y, X1) {
  nbeta1 <- ncol(X1)
  beta1  <- theta[1:nbeta1]
  sigma  <- exp( theta[nbeta1+1] )
  mu     <- itrans.mu( X1 %*% beta1 )
  -sum( dsimplex(x=Y, mu=mu, dispersion=sigma, log=TRUE) )
}
# log-likelihood function MIXED EFFECTS
llM <- function(theta, Y, X1, subject, quad) {
  N <- nlevels(subject)
  nbeta1 <- ncol(X1)
  beta1  <- theta[1:nbeta1]
  sigma  <- exp( theta[nbeta1+1] )
  t1     <- exp( theta[nbeta1+2] )
  t2     <- exp( theta[nbeta1+3] )
  i <- 1:N
  ll <- lapply(i, llind, Y, X1, sigma, subject, beta1, t1, t2)
  -sum(log(unlist(ll)))
}
# llfunc function
llind <- function(i, Y, X1, sigma, subject, beta1, t1, t2) {
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
# Inverse of 2x2 matrix
mysolve <- function(M) {
  determ <- M[1]*M[4] - M[2]^2
  temp.M1 <- M[1]
  M[1] <- -M[4]
  M[4] <- -temp.M1
  -M / determ
}

# Functions to transform rho
transf.rho     <- function(rho) log((rho+1)/(1-rho))
inv.transf.rho <- function(x)   -1+2*(exp(x)/(exp(x)+1))
#===============================================================================
#======================== END OF FUNCTIONS =====================================
#===============================================================================

# EYES dataset from Qui et al (2008)
#-----------------------------------------------------------------
# Verifying the required packages
list.of.packages <- c("statmod","GHQp","VGAM")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#-----------------------------------------------------------------
# Required packages
require(statmod)
require(GHQp)
#-----------------------------------------------------------------
# Loading the dataset
#link <- "http://dl.dropboxusercontent.com/u/9602072/Eye.txt"
#datos <- read.table(file=link, header=T, sep="")

formula <- Proportion ~ log(Days) + I(log(Days)^2) + Concentration.encoded
X1 <- model.matrix(formula, data=datos)
Y  <- datos$Proportion
subject <- as.factor(datos$Id)
zero.obs <- datos$Proportion %in% c(0,1)
X1 <- X1[!zero.obs,]
Y  <-  Y[!zero.obs]
subject <- subject[!zero.obs]

itrans.mu <- function(x) atan(x)/pi + 0.5
#-----------------------------------------------------------------
# Fitting with fixed effects only two obtain start values
#-----------------------------------------------------------------
require(VGAM)
aux <- nlminb(c(3,0,-0.3,0.5,1), llF, Y=Y, X1=X1,
              control=list(eval.max=1000, iter.max=1000))
aux$par
theta0 <- c( aux$par, -0.69, -0.69)
names(theta0) <- c(colnames(X1), "log(sigma)", "log(t1)", "log(t2)")
#-----------------------------------------------------------------
# Fitting with random effects
#-----------------------------------------------------------------
n.points <- 3
quad <- GHQ(n=n.points, ndim=2, pruning=FALSE)
t3 <- system.time( fit3 <- nlminb(theta0, llM, Y=Y, X1=X1,
                                  subject=subject, quad=quad,
                                  control=list(eval.max=10000,
                                               iter.max=10000,trace=1) ) )
n.points <- 5
quad <- GHQ(n=n.points, ndim=2, pruning=FALSE)
t5 <- system.time( fit5 <- nlminb(theta0, llM, Y=Y, X1=X1,
                                  subject=subject, quad=quad,
                                  control=list(eval.max=10000,
                                               iter.max=10000,trace=1) ) )
n.points <- 7
quad <- GHQ(n=n.points, ndim=2, pruning=FALSE)
t7 <- system.time( fit7 <- nlminb(theta0, llM, Y=Y, X1=X1,
                                  subject=subject, quad=quad,
                                  control=list(eval.max=10000,
                                               iter.max=10000,trace=1) ) )
n.points <- 3
quad <- GHQ(n=n.points, ndim=2, pruning=TRUE)
t3P <- system.time( fit3P <- nlminb(theta0, llM, Y=Y, X1=X1,
                                    subject=subject, quad=quad,
                                    control=list(eval.max=10000,
                                                 iter.max=10000,trace=1) ) )
n.points <- 5
quad <- GHQ(n=n.points, ndim=2, pruning=TRUE)
t5P <- system.time( fit5P <- nlminb(theta0, llM, Y=Y, X1=X1,
                                    subject=subject, quad=quad,
                                    control=list(eval.max=10000,
                                                 iter.max=10000,trace=1) ) )
n.points <- 7
quad <- GHQ(n=n.points, ndim=2, pruning=TRUE)
t7P <- system.time( fit7P <- nlminb(theta0, llM, Y=Y, X1=X1,
                                    subject=subject, quad=quad,
                                    control=list(eval.max=10000,
                                                 iter.max=10000,trace=1) ) )
#-----------------------------------------------------------------------
# To replicate Table 2 of paper
#-----------------------------------------------------------------------
results <- rbind( cbind(fit3$par, fit5$par, fit7$par, fit3P$par, fit5P$par, fit7P$par),
                  Time=c(t3[1], t5[1], t7[1], t3P[1], t5P[1], t7P[1]) )
colnames(results) <- c('k=3','k=5','k=7','k=3P','k=5P','k=7P')
round(results, 4)
#-----------------------------------------------------------------------

