library(ZOIP)
N<-21 # Numeros de grupos o sujetos

Times <- c(2, 10, 20, 40) # cantidad de dias

subject <- rep(1:N, each=length(Times)) 
# numero de sujetos en la muestra repetidos tantas veces haya dias

Days <- rep(Times, times=N)
b0i <- rep(rnorm(n=N,sd=1), each=length(Times))
b1i <- rep(rnorm(n=N,sd=0.5), each=length(Times))

neta <- (1.6+b0i) -1.3*log(Days)
neta2 <-(0.1+b1i) -0.8*log(Days)


mu <- 1 / (1 + exp(-neta))
sigma <- 1 / (1 + exp(-neta2))

p0 <- 0.1
p1 <- 0.1

mu[mu==1] <- 0.999
mu[mu==0] <- 0.001

sigma[sigma==1] <- 0.999
sigma[sigma==0] <- 0.001

family<-'R-S'

Y <- rZOIP(n=length(mu), mu = mu, sigma = sigma , 
           p0=p0, p1=p1, family = family) 
base<-data.frame(Y, Days, subject)

formula.mu <- Y ~ log(Days)
formula.sigma <- ~ log(Days)
formula.p0 <- ~ 1
formula.p1 <- ~ 1

formula.random = ~1 | subject

link <- c('logit', 'logit', 'identity', 'identity')

optimizer <- 'nlminb'
n.points <-11
pruning <- TRUE

mod<-RMM.ZOIP(formula.mu = formula.mu, formula.sigma = formula.sigma, 
              formula.p0 = formula.p0, formula.p1 = formula.p1,
              data=base, formula.random = formula.random, link = link,
              family = family, optimizer = optimizer,
              n.points = n.points, pruning = pruning)
mod
summary(mod)
