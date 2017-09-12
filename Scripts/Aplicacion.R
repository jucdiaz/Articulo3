if (!require('devtools')) install.packages('devtools')
devtools::install_github('jucdiaz/ZOIP', force=TRUE)
library(ZOIP)



data_bank<-read.csv(file.choose(),header=T)

head(data_bank)
colnames(data_bank)

y_i<-data_bank$Porcen_USO
data_bank$CUPO_TDC_ENTIDAD
data_bank$TOTAL_CARTERA
data_bank$SCORE

# data<-as.data.frame(cbind(y_i,CUPO_TDC_ENTIDAD=data_bank$CUPO_TDC_ENTIDAD,
#                           TOTAL_CARTERA=data_bank$TOTAL_CARTERA,PROM_Cuotas=data_bank$PROM_Cuotas))
# data<-as.data.frame(cbind(y_i,SCORE=data_bank$SCORE/1000,PROM_Cuotas=data_bank$PROM_Cuotas
#                           ,CUPO_TDC_ENTIDAD=log(data_bank$CUPO_TDC_ENTIDAD+1)))

data<-as.data.frame(cbind(y_i,TOTAL_MORA=log(data_bank$TOTAL_MORA+1),cant_meses=log(data_bank$Meses_Emision_base+1),
                          ingresos=log(data_bank$Valor_Ingresos+1)))
#TOTAL_MORA= Total Mes en Mora: Número de meses que ha estado en mora la tdc (durante toda la vida de la tarjeta).
head(data)
dim(data)
formula.mu=y_i~TOTAL_MORA
formula.sigma=~TOTAL_MORA
formula.p0=~1
formula.p1=~1
link=c('logit','logit','identity','identity')
family='R-S'
system.time(mod<-RM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,formula.p1=formula.p1,data=data,link=link,family=family))
summary(mod)
mod


data_mix<-as.data.frame(cbind(y_i,TOTAL_MORA=log(data_bank$TOTAL_MORA+1),Ciudad=data_bank$CIUDAD))

head(data_mix)

table(data_mix$Ciudad)

par(mfrow=c(1,2))
plot(density(data_mix$TOTAL_MORA))

formula.mu=y_i~TOTAL_MORA
formula.sigma=~TOTAL_MORA
formula.p0=~1
formula.p1=~1
formula.random= ~ 1 | Ciudad

link=c('logit','logit','identity','identity')
family='R-S'

optimizer<-'nlminb'
n.points <-11
pruning <- TRUE

mod<-RMM.ZOIP(formula.mu=formula.mu,formula.sigma=formula.sigma,formula.p0=formula.p0,
              formula.p1=formula.p1,data=data_mix,formula.random=formula.random,link=link,
              family=family,optimizer=optimizer,n.points=n.points,pruning=pruning)

library(MASS)

(ajuste <- fitdistr(data_mix$TOTAL_MORA,"exponential"))
plot(density(rexp(1000,ajuste$estimate)))

par(mfrow=c(1,1))
plot(data_mix$Dias_mora,data_mix$y_i, col=data_mix$Ciudad)

plot(density(data_mix$y_i))
# ajuste <- fitdistr(data_mix$Dias_mora+1,"weibull")
# ajuste
# 
# plot(density(rweibull(1000,1.89569063,2.17982587)-1))
# 
# ajuste <- fitdistr(data_mix$Dias_mora+1,"log-normal")
# ajuste
# 
# plot(density(rweibull(1000,0.50232353,0.53323966)-1))
