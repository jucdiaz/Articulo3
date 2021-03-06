library(lattice)

datos_sim<-read.table(file.choose(),sep=",")
head(datos_sim)
colnames(datos_sim)<-c('N','ni','n_points','pruning','mu_b0','mu_b1','sigma_b0',
                       'sigma_b1','p0_b0','p1_b0','random_mu','random_sigma',
                       'logveromilitud','time','num_iter')
head(datos_sim)
datos_sim$pruning2[datos_sim$pruning==0]<-'Con Pruning'
datos_sim$pruning2[datos_sim$pruning==1]<-'Sin Pruning'
head(datos_sim)
attach(datos_sim)

b0i <- 0.51
b1i <- 0.4

mu.beta0<--1.13
mu.beta1<-0.33
  
sigma.beta0<-0.33
sigma.beta1<-0.14

p0 <- 0.23
p1 <- 0.07


Plot_mape_mix<-function(beta_sim,beta_real,main,xlab,ylab){
  
  n2<-c(5,20,50)
  n.points2<-c(3,10,20)
  pruning3<-c(FALSE,TRUE)
  
  dat<-expand.grid(n2,n.points2,pruning3)
  colnames(dat)<-c('ni','n_points','pruning')
  
  dat$pruning2[dat$pruning==0]<-'Con Pruning'
  dat$pruning2[dat$pruning==1]<-'Sin Pruning'
  
  base<-tapply(abs((beta_sim-beta_real)/beta_real),list(ni,n_points,pruning),median)
  
  i<-1
  r<-1
  while(i<=length(pruning3)){
    j=1
    while(j<=length(n.points2)){
      k=1
      while(k<=length(n2)){
        dat$valor[r]<-base[k,j,i]
        k=k+1
        r=r+1
      }
      j=j+1
    }
    i=i+1
  }

xyplot(valor ~ ni | as.factor(n_points), data = dat,group=pruning2,
         strip = strip.custom(factor.levels=c('Q=3','Q=10','Q=20'),par.strip.text = list(aspect = "iso",cex = 0.75)),
         type='b',layout = c(3,1),lty=c(1,2),col=c(1,2),lwd=c(1,1),pch=13,main=main,xlab=xlab,ylab=ylab,
         key = list(text=list(levels(as.factor(dat$pruning2))),
                    space = "right",
                    points=list(pch=13),
                    lines=list(lty=c(1,2),col=c(1,2),lwd=c(1,1),pch=13)
                    )
         )

}

Plot_mape_mix(beta_sim=mu_b0,beta_real=mu.beta0,main=expression(paste('Estimaci�n de ',beta,'0 de ',mu)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')
Plot_mape_mix(beta_sim=mu_b1,beta_real=mu.beta1,main=expression(paste('Estimaci�n de ',beta,'1 de ',mu)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')
Plot_mape_mix(beta_sim=sigma_b0,beta_real=sigma.beta0,main=expression(paste('Estimaci�n de ',beta,'0 de ',sigma)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')
Plot_mape_mix(beta_sim=sigma_b1,beta_real=sigma.beta1,main=expression(paste('Estimaci�n de ',beta,'1 de ',sigma)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')
Plot_mape_mix(beta_sim=p0_b0,beta_real=p0,main=expression(paste('Estimaci�n de ',beta,'0 de ',p0)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')
Plot_mape_mix(beta_sim=p1_b0,beta_real=p1,main=expression(paste('Estimaci�n de ',beta,'0 de ',p1)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')
Plot_mape_mix(beta_sim=random_mu,beta_real=b0i,main=expression(paste('Estimaci�n de ',lambda,'1 de ',mu)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')
Plot_mape_mix(beta_sim=random_sigma,beta_real=b1i,main=expression(paste('Estimaci�n de ',lambda,'2 de ',sigma)),xlab='Tama�o muestra por ciudad', ylab='Mediana del error relativo')


Plot_mix<-function(var,main,xlab,ylab){
  
  n2<-c(5,20,50)
  n.points2<-c(3,10,20)
  pruning3<-c(FALSE,TRUE)
  
  dat<-expand.grid(n2,n.points2,pruning3)
  colnames(dat)<-c('ni','n_points','pruning')
  
  dat$pruning2[dat$pruning==0]<-'Con Pruning'
  dat$pruning2[dat$pruning==1]<-'Sin Pruning'
  
  base<-tapply(var,list(ni,n_points,pruning),median)
  
  i<-1
  r<-1
  while(i<=length(pruning3)){
    j=1
    while(j<=length(n.points2)){
      k=1
      while(k<=length(n2)){
        dat$valor[r]<-base[k,j,i]
        k=k+1
        r=r+1
      }
      j=j+1
    }
    i=i+1
  }

  xyplot(valor ~ ni | as.factor(n_points), data = dat,group=pruning2,
	
         strip = strip.custom(factor.levels=c('Q=3','Q=10','Q=20'),par.strip.text = list(aspect = "iso",cex = 0.75)),
         type='b',layout = c(3,1),lty=c(1,2),col=c(1,2),lwd=c(1,1),pch=13,main=main,xlab=xlab,ylab=ylab,
         key = list(text=list(levels(as.factor(dat$pruning2))),
                    space = "right",
                    points=list(pch=13),
                    lines=list(lty=c(1,2),col=c(1,2),lwd=c(1,1),pch=13)
         )
  )
  
}

Plot_mix(logveromilitud,main='',xlab='Tama�o muestra por ciudad', ylab='Log-verosimilitud')
Plot_mix(time,main='',xlab='Tama�o muestra por ciudad', ylab='Mediana tiempo de ejecuci�n (Seg)')
Plot_mix(num_iter,main='',xlab='Tama�o muestra por ciudad', ylab='Mediana nro de iteraciones')

tapply(mu_b0,list(n_points,pruning),median)
tapply(mu_b1,list(n_points,pruning),median)
tapply(sigma_b0,list(n_points,pruning),median)
tapply(sigma_b1,list(n_points,pruning),median)
tapply(p0_b0,list(n_points,pruning),median)
tapply(p1_b0,list(n_points,pruning),median)
tapply(random_mu,list(n_points,pruning),median)
tapply(random_sigma,list(n_points,pruning),median)
tapply(time,list(n_points,pruning),median)
tapply(num_iter,list(n_points,pruning),median)

tapply(mu_b0,list(ni,pruning),median)
tapply(mu_b1,list(ni,pruning),median)
tapply(random_mu,list(ni,pruning),median)
tapply(sigma_b0,list(ni,pruning),median)
tapply(sigma_b1,list(ni,pruning),median)
tapply(random_sigma,list(ni,pruning),median)
tapply(p0_b0,list(ni,pruning),median)
tapply(p1_b0,list(ni,pruning),median)
tapply(time,list(ni,pruning),median)
tapply(num_iter,list(ni,pruning),median)

