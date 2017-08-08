#--------------------------------------------------------------
setwd("C:/Users/Freddy/Dropbox/Investigacion/IAGQ/Paper examples/Simulations")
source("C:/Users/Freddy/Dropbox/Investigacion/IAGQ/Paper examples/Simulations/functions.R")
#--------------------------------------------------------------

Nrep <- 1000

N <- 15
simulation2(Nrep=Nrep, N=N, n.points=11, pruning=T)
simulation2(Nrep=Nrep, N=N, n.points=15, pruning=T)
simulation2(Nrep=Nrep, N=N, n.points=21, pruning=T)
simulation2(Nrep=Nrep, N=N, n.points=11, pruning=F)
simulation2(Nrep=Nrep, N=N, n.points=15, pruning=F)
simulation2(Nrep=Nrep, N=N, n.points=21, pruning=F)

N <- 21
simulation2(Nrep=Nrep, N=N, n.points=11, pruning=T)
simulation2(Nrep=Nrep, N=N, n.points=15, pruning=T)
simulation2(Nrep=Nrep, N=N, n.points=21, pruning=T)
simulation2(Nrep=Nrep, N=N, n.points=11, pruning=F)
simulation2(Nrep=Nrep, N=N, n.points=15, pruning=F)
simulation2(Nrep=Nrep, N=N, n.points=21, pruning=F)

#--------------------------------------------------------------
#--------------------------------------------------------------
