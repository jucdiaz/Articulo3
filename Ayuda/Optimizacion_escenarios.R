

data<-read.table(file.choose(),sep=',')
colnames(data)<-c('N','ni','n.points','Pruning','intercept.mu','b1.mu',
                  'intercept.sigma','b1.sigma','p0','p1','inter.ran.mu',
                  'inter.ran.sigma','logvero','time','iter')
head(data)
I<-rep(1:60,each=10)
datos<-data.frame(I,time=data$time,Esce=paste(data$ni,data$n.points,data$Pruning))
head(datos)
base<-tapply(datos$time,datos$Esce,mean)

base<-data.frame(tiempo=base,esce=1:60)
rownames(base)<-1:60

N<-10
n<-c(5,10,20,30,50,75)
n.points<-c(3,5,10,15,25)
pruning<-c(TRUE,FALSE)

escenarios<-expand.grid(N,n,n.points,pruning)
colnames(escenarios)<-c('N','ni','n.points','pruning')

(base_ord <- base[ do.call(order, base), ])

base_ord2<-data.frame(base_ord,time_estim=(((base_ord$tiempo*1000)/60)/60)/24)

i=1
suma=base_ord2$time_estim[1]
j=1
while(i<=60){
  if(suma<=22){
  base_ord2$grupo[i]<-j
  aux<-base_ord2$time_estim[i+1]
  suma<-suma+aux
  i=i+1
  }else{
    suma=base_ord2$time_estim[i]
    j=j+1
  }
}
