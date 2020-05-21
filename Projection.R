
###################################################################################################
##---------------K constant----------------##
##################################################################################################

rm(list = ls(all=T))
Don1 <- read.csv(file="Ndeaths.csv",header = T)
Don <- Don1[1:45,]
SS<-getInitial(TotDeaths~SSlogis(time,alpha,xmid,scale),data=Don)
#y <- dailyoD$cumdO[1:21]
y1 <- Don$TotDeaths
#y <- Don$inHosp
#times <- timoD[1:21]
times <- Don$time
#times <- 1:length(y1)

#we used a different parametrization
K_start<-SS["alpha"]
#K_start<-y[1]
R_start<-1/SS["scale"]
N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
#the formula for the model
log_formula<-as.formula(TotDeaths~K*N0*exp(R*times)/(K+N0*(exp(R*times)-1)))

library(nls2)
#fit the model
m<-nls2(log_formula,data=Don,start=list(K=K_start,R=R_start,N0=N0_start))
Y.predi <- round(predict(m),1)
summary(m)
plot(times,y1,xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario")
lines(times,predict(m),col="blue",lty=2,lwd=3)
legend("topleft", legend=c("DeathObs", "DeathsPredi"),
       lty=c(1,1), pch=c(1,5), col=c("black", "blue"))


#new.data <- data.frame(time = c(35,36,37,38,39,40,41))
new.data <- data.frame(time = c(46:53))
p1 <- predict(m, new.data)
yp1 <- c(Y.predi,p1)
predict(m, new.data,interval = "confidence")
library(xtable)

library(nlstools)

cc <- confint2(m)
xtable(cc)

##-----------Plotting with predicted----------------------------#
alpha <- coef(m)  #extracting coefficients
K <- alpha["K"]
R <- alpha["R"]
N0 <- alpha["N0"]

##--lower limits
Kl <-cc[1]
Rl <- cc[2]
N0l <- cc[3]
##----Upper limits
Ku <- cc[4]
Ru <- cc[5]
N0u <- cc[6]
plot(TotDeaths ~ time, data = Don, main = "Fitted logistic growth model with CCC", 
     xlab = "Days", ylab = "Deaths", xlim = c(0, 53), ylim = c(0, 2000))  # Census data
curve(K*N0*exp(R*x)/(K+N0*(exp(R*x)-1)) , add = T, col = "red")  # Fitted model
curve(Kl*N0l*exp(Rl*x)/(Kl+N0l*(exp(Rl*x)-1)) , add = T, col = "green",lty=2)  # Fitted model
curve(Ku*N0u*exp(Ru*x)/(Ku+N0u*(exp(Ru*x)-1)) , add = T, col = "green",lty=2)  # Fitted model


points(35, 1123.241 , pch = "*", cex = 1.3,col="blue")
points(36, 1157.225, pch = "*", cex = 1.3,col="blue")
points(37, 1188.671 , pch = "*", cex = 1.3,col="blue")
points(38, 1217.603, pch = "*", cex = 1.3,col="blue")
points(39, 1244.080 , pch = "*", cex = 1.3,col="blue")

library(MASS)
confint(m)

library(minpack.lm)
mm <- nlsLM(log_formula,data=Don,start=list(K=K_start,R=R_start,N0=N0_start))
confint2(m)
predict(m,new.data)
predict(mm,new.data)

library(nls2)
nls2(log_formula,data=Don,start=list(K=K_start,R=R_start,N0=N0_start))


confint2(m)

O2K.boot1 <- nlsBoot(m, niter = 100)
plot(m)
plot(O2K.boot1, type = "boxplot", ask = FALSE)
summary(O2K.boot1)

preview(log_formula, Don,list(K=K_start,R=R_start,N0=N0_start))
O2K.nls1 <- nls(log_formula, start=list(K=K_start,R=R_start,N0=N0_start), data = Don)
overview(O2K.nls1)
plotfit(O2K.nls1, smooth = TRUE)


predict(m, new.data,interval = "confidence",level = 0.95)
plotfit(m,smooth = TRUE)

#############################################################################################
##------------With K estimated--------------------------##
############################################################################################

Y <- function(x,K,N,a,r,G){
aa= (N*G*K/(K+(G-K)*exp(-a*x)))/(N+((G*K)/(K+(G-K)*exp(-a*t))-N)*(exp(-r*x)))
return(aa)
}

#log_formula3<-formula(y1~(N*2500*K/(K+(2500-K)*exp(-a*times)))/(N+((2500*K)/(K+(2500-K)*exp(-a*times))-N)*(exp(-r*times))))
#fit the model
#m3<-nls(TotalDeaths~(N*2500*K/(K+(2500-K)*exp(-a*time)))/(N+((2500*K)/(K+(2500-K)*exp(-a*time))-N)*(exp(-r*time))),data=Don,start=init)

#require(minpack.lm)
library(nls2)
G=16608
#G=12245
#init <- list(N=30, K=900,a=0.19,r=0.16)
init <- list(N=30, K=900,a=0.14,r=0.12)
init1 <- list(N0=30, G0=900,c=0.14,r=0.12)
xx <- TotDeaths~(N*G*K/(K+(G-K)*exp(-a*time)))/(N+((G*K)/(K+(G-K)*exp(-a*time))-N)*(exp(-r*time)))
#xx1= TotDeaths~(G/(1+(((G/G0)-1)*(r/(r-c)*exp(-c*time)))+(((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*time)))
xx1 <- TotDeaths ~ (G/(1+((((G/G0)-1)*(r/(r-c)))*exp(-c*time))+
                         ((((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*time))))
fit3.1 <- nls2(xx1,data=Don,start=init1)
#fit3 <- nlsLM(xx,data=Don,start=init)
summary(fit3.1)
#y3 <- round(predict(fit3),1)
y3.1 <- predict(fit3.1)
#cf <-confint(fit3.1)
cf2 <- confint2(fit3.1)
cc1<- round(confint2(fit3.1),1)
#new.data3 <- data.frame(time = c(35,36,37,38,39,40,41))
new.data3 <- data.frame(time = c(46:53))
#p3 <- predict(fit3, new.data3)
p3.1 <- predict(fit3.1, new.data3)
#yp3 <- c(y3,p3)
yp3.1 <- c(y3.1,p3.1)

##-----------Plotting with predicted----------------------------#
beta <- coef(fit3.1)  #extracting coefficients
N0 <- beta["N0"]
G0 <- beta["G0"]
r <- beta["r"]
c <- beta["c"]
b <- c(N0,G0,r,c)
#a <- c(32.70516525,448.42255262, 0.17403172,0.03832116)
##--lower limits
tm <- Don$time
N0l <-cf2[1]
G0l <- cf2[2]
cl <- cf2[3]
rl <- cf2[4]
ll <- c(N0l,G0l,cl,rl)
##----Upper limits
N0u <- cf2[5]
G0u <- cf2[6]
cu <- cf2[7]
ru <- cf2[8]
u <- c(N0u,G0u,cu,ru)
conf.band <- function(param,x){
  N0 <- param[1]
  G0 <- param[2]
  c <- param[3]
  r <- param[4]
  y<- G/(1+((((G/G0)-1)*(r/(r-c)))*exp(-c*x))+((((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*x))) 
  return(y)
}
Uper <- conf.band(param = u,x=tm)
Lw <- conf.band(param = ll,x=tm)
bb <- conf.band(param = a,x=tm)

plot(TotDeaths ~ time, data = Don, main = "Fitted logistic growth model with vcc, G=16608", 
     xlab = "Days", ylab = "Deaths", xlim = c(0, 53),ylim=c(0,2500))  # Census data
curve((G/(1+((((G/G0)-1)*(r/(r-c)))*exp(-c*x))+((((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*x)))) ,
      add = T, col = "red")  # Fitted model
aa=curve((G/(1+((((G/G0l)-1)*(rl/(rl-cl)))*exp(-cl*x))+((((G/N0l)-1)-((G/G0l)-1)*(rl/(rl-cl)))*exp(-rl*x)))) ,
      add = T, col = "green",lty=2)  # Fitted model
curve( (G/(1+((((G/G0u)-1)*(ru/(ru-cu)))*exp(-cu*x))+((((G/N0u)-1)-((G/G0u)-1)*(ru/(ru-cu)))*exp(-ru*x)))) ,
       add = T, col = "green",lty=2)  # Fitted model

##------------------------------------------------------------------------------------------------

G1=5000
G=G1
init <- list(N=30, K=1452,a=0.14,r=0.14)
#init <- list(N=30, K=1452,a=0.16,r=0.14)
init1 <- list(N0=30, G0=900,c=0.14,r=0.12)
#xx2 <- TotDeaths~(N*G1*K/(K+(G1-K)*exp(-a*time)))/(N+((G1*K)/(K+(G1-K)*exp(-a*time))-N)*(exp(-r*time)))
#fit31 <- nlsLM(xx2,data=Don,start=init,control = nls.control(maxiter = 500))
fit31.1 <- nls2(xx1,data=Don,start=init1,control = nls.control(maxiter = 500))
summary(fit31.1)
y31 <- predict(fit31)
 p31<- predict(fit31, new.data3)
 y31.1 <- predict(fit31.1)
 p31.1<- predict(fit31.1, new.data3) 
yp31 <- c(y31,p31)
yp31.1 <- c(y31.1,p31.1)
cf21 <- confint2(fit31.1)
beta1 <- coef(fit31.1)  #extracting coefficients
N0 <- beta1["N0"]
G0 <- beta1["G0"]
r <- beta1["r"]
c <- beta1["c"]

##--lower limits
tm <- Don$time
N0l <-cf21[1]
G0l <- cf21[2]
cl <- cf21[3]
rl <- cf21[4]
##----Upper limits
N0u <- cf21[5]
G0u <- cf21[6]
cu <- cf21[7]
ru <- cf21[8]

plot(TotDeaths ~ time, data = Don, main = "Fitted logistic growth Model with vcc, G=5000", 
     xlab = "Days", ylab = "Deaths", xlim = c(0, 53),ylim=c(0,2500))  # Census data
curve((G/(1+((((G/G0)-1)*(r/(r-c)))*exp(-c*x))+((((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*x)))) ,
      add = T, col = "red")  # Fitted model
curve((G/(1+((((G/G0l)-1)*(rl/(rl-cl)))*exp(-cl*x))+((((G/N0l)-1)-((G/G0l)-1)*(rl/(rl-cl)))*exp(-rl*x)))) ,
         add = T, col = "green",lty=2)  # Fitted model
curve( (G/(1+((((G/G0u)-1)*(ru/(ru-cu)))*exp(-cu*x))+((((G/N0u)-1)-((G/G0u)-1)*(ru/(ru-cu)))*exp(-ru*x)))) ,
       add = T, col = "green",lty=2)  # Fitted model

##----------------------------------------------------------------------------------------------

#G2=2342
G2=2500
G=G2
init <- list(N=30, K=500,a=0.08,r=0.1)
#init1 <- list(N=28, K=900,a=0.12,r=0.14)
xx3 <- TotDeaths~(N*G2*K/(K+(G2-K)*exp(-a*time)))/(N+((G2*K)/(K+(G2-K)*exp(-a*time))-N)*(exp(-r*time)))
fit32 <- nlsLM(xx3,data=Don,start=init,control = nls.control(maxiter = 500))
init1 <- list(N0=30, G0=500,c=0.09,r=0.15)
fit32.1 <- nls2(xx1,data=Don,start=init1,control = nls.control(maxiter = 500))

summary(fit32.1)
y32<- predict(fit32)
p32<- predict(fit32, new.data3)
yp32 <- c(y32,p32)
cf3 <- confint2(fit32.1)

y32.1<- predict(fit32.1)
p32.1<- predict(fit32.1, new.data3)
yp32.1 <- c(y32.1,p32.1)
cf22 <- confint2(fit32.1)

beta2 <- coef(fit32.1)  #extracting coefficients
N0 <- beta2["N0"]
G0 <- beta2["G0"]
r <- beta2["r"]
c <- beta2["c"]

##--lower limits
tm <- Don$time
N0l <-cf3[1]
G0l <- cf3[2]
cl <- cf3[3]
rl <- cf3[4]
##----Upper limits
N0u <- cf3[5]
G0u <- cf3[6]
cu <- cf3[7]
ru <- cf3[8]

plot(TotDeaths ~ time, data = Don, main = "Fitted logistic Growth Model with vcc, G=2500", 
     xlab = "Days", ylab = "Deaths", xlim = c(0, 53),ylim=c(0,2500))  # Census data
curve((G/(1+((((G/G0)-1)*(r/(r-c)))*exp(-c*x))+((((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*x)))) ,
      add = T, col = "red")  # Fitted model
aa=curve((G/(1+((((G/G0l)-1)*(rl/(rl-cl)))*exp(-cl*x))+((((G/N0l)-1)-((G/G0l)-1)*(rl/(rl-cl)))*exp(-rl*x)))) ,
         add = T, col = "green",lty=2)  # Fitted model
curve( (G/(1+((((G/G0u)-1)*(ru/(ru-cu)))*exp(-cu*x))+((((G/N0u)-1)-((G/G0u)-1)*(ru/(ru-cu)))*exp(-ru*x)))) ,
       add = T, col = "green",lty=2)  # Fitted model


#hh <-c(Don$TotDeaths,kk)

plot(times,y1,xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario")
lines(times,predict(fit3),col="green",lty=1,lwd=3)
legend("topleft", legend=c("DeathObs", "DeathsPredi"),
       lty=c(1,1), pch=c(1,5), col=c("black", "green"))

plot(times,y1,xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario")
lines(times,predict(m2),col="red",lty=2,lwd=3)
lines(times,predict(m),col="blue",lty=2,lwd=3)
lines(times,predict(fit3),col="green",lty=2,lwd=3)
legend("topleft", legend=c("Kconstant", "Kvarying","K_logistic"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red","green"))

#-------------------------------------------------------------------------------------------------
##-----------------All parameters estimated-----------------------------------#

library(nls2)
#st2 <- data.frame(N=c(20,30), K=c(600,900),a=c(0.19,0.16),r=c(0.14,0.12),G=c(3000,20000))
#st1 <- expand.grid(N = seq(18, 30, len = 4), r=seq(0.18, 0.4, len = 4),
#                   K = seq(500, 900, len = 4), a = seq(0.1, 0.2, len = 4),G = seq(1000, 2000, len = 4))
#y <- Don$TotalDeaths
#tm <- Don$time
#aa= TotDeaths~(N*G*K/(K+(G-K)*exp(-a*time)))/(N+((G*K)/(K+(G-K)*exp(-a*time))-N)*(exp(-r*time)))
#m3<-nls2(aa,data = Don,start = st2, algorithm = "brute-force")
#st <- list(N=30,K=900,a=0.19,r=0.14,G=5000)
#aa= TotDeaths~(N*G*K/(K+(G-K)*exp(-a*time)))/(N+((G*K)/(K+(G-K)*exp(-a*time))-N)*(exp(-r*time)))
#m6<- nlsLM(aa,data=Don,start =st,control = nls.control(maxiter = 1000))

st1 <- list(N0=30,G0=700,c=0.08,r=0.16,G=3000)

aa1= TotDeaths~(G/(1+((G/G0)-1)*(r/(r-c)*exp(-c*time))+(((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*time)))
m6b<- nls2(aa1,data=Don,start =st1,control = nls.control(maxiter = 1000))
summary(m6b)


y4 <-predict(m6b)
new.data3 <- data.frame(time = 46:60)
p4 <- predict(m6b, new.data3)
yp4 <- c(y4,p4)
cf4 <- confint2(m6b)


beta3 <- coef(m6b)  #extracting coefficients
N0 <- beta3["N0"]
G0 <- beta3["G0"]
r <- beta3["r"]
c <- beta3["c"]
G  <- beta3["G"]

##--lower limits
tm <- Don$time
N0l <-cf4[1]
G0l <- cf4[2]
cl <- cf4[3]
rl <- cf4[4]
Gl <-cf4[5]
##----Upper limits
N0u <- cf4[6]
G0u <- cf4[7]
cu <- cf4[8]
ru <- cf4[9]
Gu <- cf4[10]

plot(TotDeaths ~ time, data = Don1, main = "Fitted logistic Growth Model with vcc", 
     xlab = "Days", ylab = "Deaths", xlim = c(0, 53),ylim=c(0,2500),cex = 1.5)  # Census data
curve((G/(1+((((G/G0)-1)*(r/(r-c)))*exp(-c*x))+((((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*x)))) ,
      add = T, col = "red",lwd=2)  # Fitted model
curve((Gl/(1+((((Gl/G0l)-1)*(rl/(rl-cl)))*exp(-cl*x))+((((Gl/N0l)-1)-((Gl/G0l)-1)*(rl/(rl-cl)))*exp(-rl*x)))) ,
         add = T, col = "green",lty=2,lwd=2)  # Fitted model
curve( (Gu/(1+((((Gu/G0u)-1)*(ru/(ru-cu)))*exp(-cu*x))+((((Gu/N0u)-1)-((Gu/G0u)-1)*(ru/(ru-cu)))*exp(-ru*x)))) ,
       add = T, col = "green",lty=2,lwd=2)  # Fitted model

points(46, 1790.6 , pch = "*", cex = 1.5,col="blue")
points(47, 1840.0, pch = "*", cex = 1.5,col="blue")
points(48, 1887.8 , pch = "*", cex = 1.5,col="blue")
points(49, 1933.9, pch = "*", cex = 1.5,col="blue")
points(50, 1978.2 , pch = "*", cex = 1.5,col="blue")
points(51, 2020.6 , pch = "*", cex = 1.5,col="blue")
points(52, 2061.2 , pch = "*", cex = 1.5,col="blue")
points(53, 2099.9 , pch = "*", cex = 1.5,col="blue")

plot(times,y1,xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario")
lines(times,predict(m6),col="yellow",lty=2,lwd=3)
legend("topleft", legend=c("DeathObs", "DeathsPredi"),
       lty=c(1,1), pch=c(1,5), col=c("black", "yellow"))

plot(times,y1,xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario")
lines(times,predict(m2),col="red",lty=2,lwd=3)
lines(times,predict(m),col="blue",lty=2,lwd=3)
lines(times,predict(fit3),col="green",lty=2,lwd=3)
lines(times,predict(m6),col="yellow",lty=2,lwd=3)
legend("topleft", legend=c("Kconstant", "Kvarying","Kestimated","KGestimated"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red","green","yellow"))

tim<- Don1$time
y <- Don1$TotDeaths

plot(tim,y,type="l",xlab = "days",ylim = c(15,2000),ylab="Total deaths", main = "Total deaths  in Ontario",lwd=3)
points(y)
lines(tim,yp1,col="red",lty=2,lwd=3)
lines(tim,yp3,col="blue",lty=2,lwd=3)
lines(tim,yp31,col="green",lty=2,lwd=3)
lines(tim,yp32,col="yellowgreen",lty=2,lwd=3)
lines(tim,yp4,col="orange3",lty=2,lwd=3)
abline(v=34,col="yellow3",lty=2)
#text(x=36,y=500,labels = "30/04/20",col=4,font = 4)
mtext(text = "30-04-20",side=1,adj=0.78,col="yellow3")
legend("topleft", legend=c("SimpleLog", "MaxG","5kG","2KG","AllEst"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red","green","yellowgreen","orange3"))



##-----------Plotting with predicted----------------------------#
al6 <- coef(m6)  #extracting coefficients
N <- al6["N"]
K <- al6["K"]
a <- al6["a"]
r <- al6["r"]
G.est <- al6["G"]
plot(TotDeaths ~ time, data = Don, main = "Logistic Growth Model of Deaths in Ontario", 
     xlab = "Days", ylab = "Deaths", xlim = c(0, 45), ylim = c(0, 1500))  # Census data
curve((N*G.est*K/(K+(G.est-K)*exp(-a*x)))/(N+((G.est*K)/(K+(G.est-K)*exp(-a*x))-N)*(exp(-r*x))), 
      add = T, col = "green")  # Fitted model
points(35, 1183.575 , pch = "*", cex = 1.3,col="red")
points(36, 1247.775, pch = "*", cex = 1.3,col="red")
points(37, 1314.967 , pch = "*", cex = 1.3,col="red")
points(38, 1385.347, pch = "*", cex = 1.3,col="red")
points(39, 1459.116 , pch = "*", cex = 1.3,col="red")




Pred.D <- data.frame(yObse=y1,y.constK=round(as.numeric(m$m$fitted()),2),
                     #y.varyK=round(as.numeric(m2$m$fitted()),2),
                     y.estK=round(as.numeric(fit3$m$fitted()),2),
                     y.estK1=round(as.numeric(fit31$m$fitted()),2),
                     y.estK2=round(as.numeric(fit32$m$fitted()),2)
                     )

pred.pre <- data.frame(yObse=c(951,996,1082,1121,1176,1216,1300,1361,NA,NA),y.constK=round(predict(m, new.data),2),
                       y.estK=round(predict(fit3, new.data3),2),
                       y.estK1=round(predict(fit31, new.data3),2),
                       y.estK2=round(predict(fit32, new.data3),2))

as.matrix(AIC(m, fit3,fit31,fit32,m6))
as.matrix(BIC(m, fit3,fit31,fit32,m6))
#-------------------------------------------------------------------------------------------
y <- Don$TotDeaths
x <- Don$time
g <-  y~G1 +( G2/(1+exp(-a*(x-b))))
ff <- nls(g,start= c(G1=1452,a=1,G2=5000,b=10))

g1 <-  y~ G2/(1+exp(-a*(x-b)))
ff1 <- nls(g1,start= c(a=1,G2=5000,b=10))
vcov(ff1)
vcov(ff)

as.matrix(AIC(ff, ff1))
as.matrix(BIC(ff, ff1))
plot(x,y)
lines(predict(ff),col="red",lwd=2)
lines(predict(ff1),col="green",lwd=2)
#################################################################################################
#################################################################################################



Don <- read.csv(file="death.csv",header = T)
library(astsa)

newdeaths <- ts(Don$TotalDeaths)
newHosp <- ts(Don$inHosp)

lag2.plot (newHosp, newdeaths, 10)

alldata=ts.intersect(newdeaths, Hosplag1 = lag(newHosp,-1),Hosplag2 = lag(newHosp,-2),
                     Hosplag3 = lag(newHosp,-3), Hosplag4 = lag(newHosp,-4), 
                     Hosplag5 = lag(newHosp,-5),Hosplag6 = lag(newHosp,-6),
                     Hosplag7 = lag(newHosp,-7))

library(useful)
#mydata7<-shift.column(data=Don, columns="inHosp", len=7, up=TRUE)
mydata71<-shift.column(data=Don, columns="TotalDeaths", len=7, up=TRUE)

library(nls2)
st3<- data.frame(N=c(30,20),r=c(0.14,0.2))
lo <- TotalDeaths.Shifted ~ (N*inHosp)/(N+(inHosp-N)*exp(-r*time))
mm<-nls2(lo,data = mydata71,start = st3, algorithm = "brute-force")

mm1<-nls(lo,data=mydata71,start=list(r=0.2,N=30))
y.pre <- predict(mm1)
new.data1 <- data.frame(inHosp = 1349, time=23)
predict(mm1, new.data1)

SS<-getInitial(inHosp.Shifted~SSlogis(time,alpha,xmid,scale),data=mydata7)
y1 <- mydata7$TotalDeaths
y <- mydata7$inHosp.Shifted



SS<-getInitial(TotalDeaths~SSlogis(time,alpha,xmid,scale),data=Don)
#y <- dailyoD$cumdO[1:21]
y <- Don$TotalDeaths
z<- Don$inHosp
#times <- timoD[1:21]
times <- Don$time
#we used a different parametrization
K_start<-SS["alpha"]
R_start<-1/SS["scale"]
N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
N0_start<- 55
#the formula for the model
log_formula<-formula(y~K*N0*exp(R*times)/(K+N0*(exp(R*times)-1)),data=Dt)
#fit the model
m<-nls(log_formula,data=Dt,start=list(K=K_start,R=R_start,N0=N0_start))
#estimated parameters
summary(m)
plot(times,y,xlab = "days",ylab="Total deaths", main = "Total deaths in Ontario")
lines(times,predict(m),col="red",lty=2,lwd=3)

new.data <- data.frame(times = c(30,37))
predict(m, new.data)

library(ggplot2)

Don[-24,] %>%
  gather(variable,value, TotalDeaths, inHosp) %>%
  ggplot(aes(x=time, y=value, colour=variable)) +
  geom_line(size=2)+
  ggtitle("Cumulative Daily Aggregated values")

dH <-diff(newHosp[-24],lag = 1)
dD <-diff(newdeaths[-24],lag = 1)
time <- 1:22
DD <- data.frame(time,dH,dD)
DD1 <- data.frame(time,a=log(dH),b=log(dD+0.0001))
DD %>%
  gather(variable,value, dD, dH) %>%
  ggplot(aes(x=time, y=value, colour=variable)) +
  geom_line(size=2)+
  ggtitle(" Daily  Reports or Difference of order 1")

DD1 %>%
  gather(variable,value, b, a) %>%
  ggplot(aes(x=time, y=value, colour=variable)) +
  geom_line()

##-------------------Richards growth equation--------------
library(growthmodels)
growth <- richard(0:10, 10, 0.5, 0.3, 0.5)

##--------------------Lag 7--------------------------
Don <- read.csv(file="death.csv",header = T)

##-------------------------------------------------------------

library(astsa)

newdeaths <- ts(Don$TotalDeaths)
newHosp <- ts(Don$inHosp)

lag2.plot (newHosp, newdeaths, 10)

alldata=ts.intersect(newdeaths, Hosplag1 = lag(newHosp,-1),Hosplag2 = lag(newHosp,-2),
                     Hosplag3 = lag(newHosp,-3), Hosplag4 = lag(newHosp,-4), 
                     Hosplag5 = lag(newHosp,-5),Hosplag6 = lag(newHosp,-6),
                     Hosplag7 = lag(newHosp,-7))

yDeaths <- as.numeric(alldata[,1])
yHosp <- as.numeric(alldata[,8])
timpo <- 1:length(yHosp)
Dt<- data.frame(y=yDeaths,z=yHosp,tim=timpo)

SS<-getInitial(y~SSlogis(tim,alpha,xmid,scale),data=Dt)
#y1 <- mydata7$TotalDeaths
y <- Dt$y
tim <-Dt$tim
#times <- timoD[1:21]
#times <- Don$time[-24]
#we used a different parametrization
#K_start<-SS["alpha"]
K_start<-SS["alpha"]
R_start<-1/SS["scale"]
N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
#the formula for the model
log_formula<-formula(y~K*N0*exp(R*tim)/(K+N0*(exp(R*tim)-1)))
#fit the model
m<-nls(log_formula,start=list(K=K_start,R=R_start,N0=N0_start))

##-------------------------------------------------

ff <- function(a,b,x,t){
  N = a*x/(a+(x-a)*exp(-b*t))
  return(N)
}
y <- Dt$y
z<- Dt$z
f<-ff(a=94,b=0.14,x=Dt$z,t=tim)

##########################################################
modeldf <- data.frame(rev=c(17906.4, 5303.72, 2700.58 ,1696.77 ,947.53 ,362.03), weeks=c(1,2,3,4,5,6))

require(minpack.lm)
fit <- nlsLM(rev ~ a*weeks^b, data=modeldf, start = list(a=1,b=1))
fit1 <- nlsLM(y~z*N0*exp(R*tim)/(z+N0*(exp(R*tim)-1)),data=Dt,start=list(R=R_start,N0=N0_start))



m1<-nls(y~z*N0*exp(R*tim)/(z+N0*(exp(R*tim)-1)),data=Dt,start=list(R=R_start,N0=N0_start))
m2<-nls(TotalDeaths~inHosp*N0*exp(R*time)/(inHosp+N0*(exp(R*time)-1)),data=Don,start=list(R=R_start,N0=N0_start))
new.data <- data.frame(z = 1200, tim=26)
predict(m1, new.data)
require(broom)
fit_data <- augment(fit)
fit_data1 <- augment(fit1)
plot(.fitted~y, data=fit_data1)

predict(fit1,newdata = c(1200,26))

y.fit <- as.numeric(m2$m$fitted())
y.obs <- Don$TotalDeaths

plot(y.obs,main = "fitted model")
lines(y.fit,col="red")

######################################################################################
##-----------------------------------------------------

##----------------------------------------------------------------

library(nls.multstart)
lgst <- function(K,N,a,r,G,t){
  aa= (N*G*K/(K+(G-K)*exp(-a*t)))/(N+((G*K)/(K+(G-K)*exp(-a*t))-N)*(exp(-r*t)))
  return(aa)
}

fitt <- nls(TotDeaths~lgst(K,N,a,r,G,t=time),data=Don,
                      start = list(K=250, N=2, a=0.09, r=0.12,G=2500)
                      )



##----------------Richard's growth model-------------------

Rich <- function(input,itian,cary,gr,scp){
  rgm <- cary*(1-exp(-scp*gr*input)*(1-(itian/cary)^(-scp)))^(1/scp)
  return(rgm)
}
rr <- formula(TotDeaths~ (A+(K-A)/(1+N0*exp(-r*time))^(1/b)))

st <- list(N0=30,A=500,b=0.1,r=0.08,K=2000)

m7<- nlsLM(rr,data=Don,start =st,control = nls.control(maxiter = 1000))
m7

library(nls2)
m8 <-  nls(rr,data=Don,start =st)

########################################################
##----Bootstrap Negative Binomial---------------------##
########################################################

yd <- diff(y4)
yd <- c(y4[1],yd)
pd <- diff(p4)
pdd <- c(50.914647,pd)

dd<- diff(Don1$TotDeaths[1:45])
dd<- c(Don1$TotDeaths[1],dd)
hist(dd)
hist(dd,probability = T,main = "Daily reported deaths")
barplot(dd)

# load library
library(fitdistrplus)
require(npsurv)
require(lsei)
library(MASS)
# fit the negative binomial distribution
fit <- fitdist(dd, "nbinom",method = "mle")
summary(fit)
fit1 <- fitdist(dd, "nbinom",method = "mme")
summary(fit1)
plot(fit)
plot(fit1)
p <- 2.860865/(2.860865+32.968695)

p1 <- 3.875167/(3.875167+36.925000)

ff <- rnbinom(34,size =2.860865,prob = p )

FF<- matrix(NA,nrow=45,ncol=1000)
for (i in 1:45) {
  for (j in 1:1000) {
    
    FF[i,j] <- rnbinom(1,mu=fit$estimate[2],size =fit$estimate[1] )
  }
}


fm <- colMeans(FF)



ff1 <- rnbinom(34,size =3.875167,prob = p1 )

#--Only simulate form predicted values
Y<- matrix(NA,nrow=1000,ncol=15)
for (j in 1:15) {
  for (i in 1:1000) {
    
    Y[i,j] <- rnbinom(1,mu=pdd[j],size =fit$estimate[1] )
  }
}

yy <- colMeans(Y)

Z<- matrix(NA,nrow=1000,ncol=60)

for (i in 1:1000) {
  
  Z[i,] <- c(yd,Y[i,])
}


ZZ<- matrix(NA,nrow=1000,ncol=60)

  for (i in 1:1000) {
    
    ZZ[i,] <- cumsum(Z[i,])
}

zz <- colMeans(ZZ)
plot(Don1$TotDeaths, type="n", ylab="Predicted/observed", xlab="Days",ylim = c(0,2600),xlim=c(0,60),
     main = "Cumulative number of Deaths in Ontario")
for(i in 1:1000){
  y <- ZZ[i,] # mme
  #y1 <- rnbinom(45,mu=38.33,size=3.135)
  lines(y, col="red")
  #lines(cumsum(y1), col="blue" )
}
points(Don1$TotDeaths, type="b")
lines(yp4,col="blue",lwd=3)
legend("topleft", legend=c("Fitted", "simu"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red"))


tim <- 1:34

library(ggplot2)
library(reshape2)

##--------Daily data
dat <- data.frame(dd,ff,tim)
df2 <- melt(dat, id.vars='tim')

ggplot(df2, aes(x=tim, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge')

ggplot(df2, aes(x=tim, y=value, fill=variable)) +
  geom_bar(stat='identity')+
 facet_wrap(~ variable) + 
scale_x_continuous(breaks=seq(1,34,2))

##----Cumulative data
obd <- Don1$TotDeaths[1:34]
bestd <-cumsum(ff)
ylogist <- Y.predi
dt1 <- data.frame(obd,ylogist,bestd,tim)
dt3 <- melt(dt1, id.vars='tim')
ggplot(dt3, aes(x=tim, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge')+
  ggtitle(" Cumulative number of reported and estimated deaths in Ontario")

ggplot(dt3, aes(x=tim, y=value, fill=variable)) +
  geom_bar(stat='identity')+
  facet_wrap(~ variable) + 
  scale_x_continuous(breaks=seq(1,34,2))+
  ggtitle(" Cumulative number of deaths in Ontario")


plot(tim,y1,type="l",xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario",lwd=3)
points(y1)
lines(tim,Y.predi,col="red",lty=1,lwd=3)
lines(tim,y3,col="blue",lty=1,lwd=3)
legend("topleft", legend=c("reported","ccc","vcc"),
       lty=c(1,1), pch=c(2,5), col=c( "black","red","blue"),
       adj = c(0,1))




fitb <- bootdist(fit1) 

summary(fitb)
plot(fitb)

fbL <- rnbinom(40,size =2.407482,prob = 0.07152888 )
fbU <- rnbinom(40,size =7.416094,prob = 0.1468792 )
plotdist(dd,discrete = TRUE)

fit2 <- fitdist(dd, "pois",method = "mle")

# get the fitted densities. mu and size from fit.
fitD <- dnbinom(0:100, size=fit$estimate[1], mu=fit$estimate[2])
fitD3 <- dnbinom(0:100, size=4.33, mu=38.33)
fitD1 <- dpois(0:100, lambda =fit1$estimate)
# add fitted line (blue) to histogram
lines(fitD, lwd="3", col="blue")
lines(fitD1, lwd="3", col="red")
lines(fitD3, lwd="3", col="green")
legend("topleft", legend=c("Obs", "NBmle","NBmmnt","Pois"),
       lty=c(1,1), pch=c(1,5), col=c("black", "blue","green","red"))

y <- rnbinom(40,mu=36.920581,size=2.870386)
y1 <- rpois(40,36.925)

plot(Don1$TotDeaths[1:40])
lines(cumsum(y),col="red")
lines(cumsum(y1),col="green")

legend("topleft", legend=c("Kconstant", "Kvarying","Kestimated","KGestimated"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red","green","yellow"))

plot(dd, type="n", ylab="Predicted/observed", xlab="Days")
for(i in 1:100){
  y <- rnbinom(40,size =3.875167,prob = 0.09497919 ) # mme
  #y1 <- rnbinom(45,mu=38.33,size=3.135)
  lines(y, col="red")
  #lines(cumsum(y1), col="green" )
}
points(dd, type="b")



plot(Don1$TotDeaths[1:40], type="n", ylab="Predicted/observed", xlab="Days",ylim = c(0,1800))
for(i in 1:100){
  #y <- Y[i,] # mme
  y1 <- rnbinom(60,size =3.875167,prob = 0.09497919 )
  #lines(y, col="blue")
  lines(cumsum(y1), col="green" )
}
points(Don1$TotDeaths[1:40], type="b", col="red")



plot(dd,type="l")
lines(Y[1,],col="red")
lines(Y[90,],col="blue")

plot(dd, type="n", ylab="Predicted/observed", xlab="Days")
for(i in 1:100){
  y <- rpois(40,36.920581)
  #y1 <- rnbinom(40,mu=36.920581,size=2.870386)
  #points(y, col=grey(.5))
  lines(y, col="green" )
  #lines(y1, col="green" )
}
points(dd, type="b", col="red")


#############################################################################################
##----------------Growth and Diffusion phenomena----------------
############################################################################################
parameters <- data.frame(c =c(0.1,0.1,0,0.025),G= c(300,180,100,20))


gth <- function(param,t){
  G <- param[1]
  G0 <- param[2]
  c <- param[3]
  r <- param[4]
  N0 <- param[5]
  m1 <- (((G/N0)-1)-((G/G0)-1)*(r/(r-c)))
  m2 <- (((G/G0)-1)*(r/(r-c)))
N <- G/(1+(m2*exp(-c*t))+(m1*exp(-r*t)))
return(N)
}

time <- 1:40
pmt1<- c(300,100,0.1,0.25,10)
gro1 <- gth(param = pmt1,t=time)

pmt2<- c(180,100,0.1,0.25,10)
gro2 <- gth(param = pmt2,t=time)

pmt3<- c(100,100,0,0.25,10)
gro3 <- gth(param = pmt3,t=time)

pmt4<- c(20,100,0.025,0.25,10)
gro4 <- gth(param = pmt4,t=time)

plot(time,gro1,type = "l",ylim = c(0,300))
lines(gro2,lty=1,col="red")
lines(gro3,lty=1,col="blue")
lines(gro4,lty=1,col="green")
##------------Logistically variable carrying capacity---------------------

Kc <- function(prm,t){
  G= prm[1]
  G0=prm[2]
  c=prm[3]
  K <- G/(1+(((G/G0)-1)*exp(-c*t)))
  return(K)
}
pm1<- c(300,100,0.1)
lkc1 <- Kc(prm=pm1,t=time)

pm2<- c(180,100,0.1)
lkc2 <- Kc(prm=pm2,t=time)

pm3<- c(100,100,0.1)
lkc3 <- Kc(prm=pm3,t=time)

pm4<- c(20,100,0.1)
lkc4 <- Kc(prm=pm4,t=time)

lines(lkc1,lty=2)
lines(lkc2,col="red",lty=2)
lines(lkc3,col="blue",lty=2)
lines(lkc4,col="green",lty=2)
##----------------------------------------------------------------------------

st <- list(N=30,K=900,a=0.19,r=0.14,G=5000)
st1 <- list(N0=30,G0=900,c=0.19,r=0.14,G=5000)

aa= TotDeaths~(N*G*K/(K+(G-K)*exp(-a*time)))/(N+((G*K)/(K+(G-K)*exp(-a*time))-N)*(exp(-r*time)))
aa1= TotDeaths~(G/(1+((G/G0)-1)*(r/(r-c)*exp(-c*time))+(((G/N0)-1)-((G/G0)-1)*(r/(r-c)))*exp(-r*time)))
m6<- nlsLM(aa,data=Don,start =st,control = nls.control(maxiter = 1000))
m6b<- nlsLM(xx1,data=Don,start =st1,control = nls.control(maxiter = 1000))

pm6 <- predict(m6)
p6 <- predict(m6,new.data3)
yp6 <- c(pm6,p6)
pm6b <- predict(m6b)
p6b <- predict(m6b,new.data3)
yp6b <- c(pm6b,p6b)
plot(pm6)
lines(pm6b,col="red")

as.matrix(AIC(m, fit3,fit31,fit32,m6))


tim<- Don1$time
tim1<- 1:50
y <- Don1$TotDeaths
d1 <- c(1765, 1798, 1825,1858,1881)
y1 <- c(y,d1)
plot(tim1,y1,type="l",xlab = "days",ylim = c(15,2000),ylab="Total deaths", main = "Total deaths  in Ontario",lwd=3)
points(y1)
lines(tim1,yp1,col="red",lty=2,lwd=3)
lines(tim1,yp3,col="blue",lty=2,lwd=3)
lines(tim1,yp31,col="green",lty=2,lwd=3)
lines(tim1,yp32,col="yellowgreen",lty=2,lwd=3)
lines(tim1,yp6,col="orange3",lty=2,lwd=3)
abline(v=34,col="yellow3",lty=2)
#text(x=36,y=500,labels = "30/04/20",col=4,font = 4)
mtext(text = "30-04-20",side=1,adj=0.78,col="yellow3")
legend("topleft", legend=c("SimpleLog", "MaxG","5kG","2KG","AllEst"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red","green","yellowgreen","orange3"))



plot(tim1,y1,type="l",xlab = "days",ylim = c(15,2000),ylab="Total deaths", main = "Total deaths  in Ontario",lwd=3)
points(y1)
#lines(tim,yp1,col="red",lty=2,lwd=3)
lines(tim1,yp3.1,col="blue",lty=2,lwd=3)
lines(tim1,yp31.1,col="green",lty=2,lwd=3)
lines(tim1,yp32.1,col="yellowgreen",lty=2,lwd=3)
lines(tim1,yp6b,col="orange3",lty=2,lwd=3)
abline(v=34,col="yellow3",lty=2)
#text(x=36,y=500,labels = "30/04/20",col=4,font = 4)
mtext(text = "30-04-20",side=1,adj=0.78,col="yellow3")
legend("topleft", legend=c("SimpleLog", "MaxG","5kG","2KG","AllEst"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red","green","yellowgreen","orange3"))

##---------------------------------------------------------------------------------------
library(ggplot2)
curves <- data.frame(time,gro1,gro2,gro3,gro4,lkc1,lkc2,lkc3,lkc4)

p1 <- ggplot() +
  geom_line(data=curves, aes(x=time, y=gro1)) +
  geom_line(data=curves, aes(x=time, y=lkc1),colour="red")+
  ggtitle("Growth curve with G=300 and c=0.1")
p2 <- ggplot() +
  geom_line(data=curves, aes(x=time, y=gro2)) +
  geom_line(data=curves, aes(x=time, y=lkc2),colour="red")+
  ggtitle("Growth curve with G=180 and c=0.1")

p3 <- ggplot() +
  geom_line(data=curves, aes(x=time, y=gro3)) +
  geom_line(data=curves, aes(x=time, y=lkc3),colour="red")+
  ggtitle("Growth curve with G=100 and c=0")

p4 <- ggplot() +
  geom_line(data=curves, aes(x=time, y=gro4)) +
  geom_line(data=curves, aes(x=time, y=lkc4),colour="red")+
  ggtitle("Growth curve with G=20 and c=0.025")

library(tidyverse)
library(lubridate)
y1 <- Don1$TotDeaths
pd <- data.frame(y1,yp1,yp3,yp31,yp32,yp6)
Don1 <- Don1 %>%
  mutate(date =mdy(Date))
tim1 <- Don1$date
pd$tim1 <- tim1
pp1 <- ggplot() +
  geom_line(data=pd, aes(x=tim1, y=y1)) +
  geom_point(aes(x=tim1, y=y1))+
  geom_line(data=pd, aes(x=tim1, y=yp1),col="red")+
  ggtitle(" constant CC")+
  theme(plot.title = element_text(size=12, face="bold"))+
  labs(y="Cum. deaths", x = "Days")+
  geom_vline(xintercept = tim1[34], linetype="dotted", 
             color = "blue")
 


pp2 <- ggplot() +
  geom_line(data=pd, aes(x=tim1, y=y1)) +
  geom_point(aes(x=tim1, y=y1))+
  geom_line(data=pd, aes(x=tim1, y=yp3),colour="blue")+
  ggtitle("  CC with G=16608")+
  theme(plot.title = element_text(size=12, face="bold"))+
  labs(y="Cum. deaths", x = "Days")+
  geom_vline(xintercept = tim1[34], linetype="dotted", 
             color = "blue")

pp3 <- ggplot() +
  geom_line(data=pd, aes(x=tim1, y=y1)) +
  geom_point(aes(x=tim1, y=y1))+
  geom_line(data=pd, aes(x=tim1, y=yp31),colour="green")+
  ggtitle(" CC with G=5000")+
  theme(plot.title = element_text(size=12, face="bold"))+
  labs(y="Cum. deaths", x = "Days")+
  geom_vline(xintercept = tim1[34], linetype="dotted", 
             color = "blue")

pp4 <- ggplot() +
  geom_line(data=pd, aes(x=tim1, y=y1)) +
  geom_point(aes(x=tim1, y=y1))+
  geom_line(data=pd, aes(x=tim1, y=yp32),colour="yellowgreen")+
  ggtitle(" CC with G=2500")+
  theme(plot.title = element_text(size=12, face="bold"))+
  labs(y="Cum. deaths", x = "Days")+
  geom_vline(xintercept = tim1[34], linetype="dotted", 
             color = "blue")

pp5 <- ggplot() +
  geom_line(data=pd, aes(x=tim1, y=y1)) +
  geom_point(aes(x=tim1, y=y1))+
  geom_line(data=pd, aes(x=tim1, y=yp6),colour="orange3")+
  ggtitle(" CC with all estimated parameters")+
  theme(plot.title = element_text(size=12, face="bold"))+
  labs(y="Cum. deaths", x = "Days")+
  geom_vline(xintercept = tim1[34], linetype="dotted", 
             color = "blue")

pp6 <- ggplot() +
  geom_line(data=pd, aes(x=tim1, y=y1)) +
  geom_point(aes(x=tim1, y=y1))+
  geom_line(data=pd, aes(x=tim1, y=yp1),colour="red")+
  geom_line(data=pd, aes(x=tim1, y=yp3),colour="blue")+
  geom_line(data=pd, aes(x=tim1, y=yp31),colour="green")+
  geom_line(data=pd, aes(x=tim1, y=yp32),colour="yellowgreen")+
  geom_line(data=pd, aes(x=tim1, y=yp6),colour="orange3")+
  ggtitle("All curves together")+
  theme(plot.title = element_text(size=12, face="bold"))+
  labs(y="Cum. deaths", x = "Days")+
  geom_vline(xintercept = tim1[34], linetype="dotted", 
           color = "blue")
  
Don1 <- Don1 %>%
  mutate(date =mdy(Date))

dd <- diff(Don1$TotDeaths)
dd1 <- c(Don1$TotDeaths[1],dd)
newdata <- Don1
newdata$dif <- dd1

 Cum <- ggplot()+
  geom_line(data=newdata,aes(x=date,y=TotDeaths))+
  geom_point(aes(x=newdata$date, y=newdata$TotDeaths),col="blue")+
  ggtitle("Cumulative Deaths in Ontario")+
  labs(y="Cum. deaths", x = "Days")

DD <- ggplot()+
  geom_line(data=newdata,aes(date,dif))+
 geom_point(aes(x=newdata$date, y=newdata$dif),col="red")+
   ggtitle("Daily Deaths in Ontario")+
   labs(y="Deaths", x = "Days")




# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(p1, p2, p3, p4)

multiplot(pp1, pp2, pp3, pp4,pp5,pp6)
multiplot(Cum,DD)


###################################################################################
#####-----------Bootstrap with nlsBoot-------------------------------#########
##################################################################################

library(nlstools)
library(nlsMicrobio)


###########################################################################################
###-------Varying K = Incorporating HOsp as K  ------------------###
###########################################################################################

#data <- read.csv(file = "death.csv",header = T)
#library(useful)
#shiftData<-shift.column(data=data, columns="inHosp", len=6, up=FALSE)

#NewDat <- as.data.frame(shiftData)
#NewDat$times<- 1:24
SS2<-getInitial(TotHosp~SSlogis(time,xmid,scale,alpha),data=Don)
K_start2<-SS2["alpha"]
#K_start<-y[1]
R_start2<-1/SS2["scale"]
N0_start2<-SS2["alpha"]/(exp(SS2["xmid"]/SS2["scale"])+1)
#the formula for the model
formula2<-formula(TotHosp~K*N0*exp(R*time)/(K+N0*(exp(R*time)-1)))
#fit the model
#fomula <-  TotDeaths~TotHosp*N0*exp(R*time)/(TotHosp+N0*(exp(R*time)-1))
#start1=list(R=R_start2,N0=N0_start2,K=K_start2)
st3 <- list(R=0.15,N0=110,K=3000)
#mm2<-nls(formula,data=Don,start=start1)
mm2<-nls(formula2,data=Don,start=st3)
y.predi2 <- predict(mm2)
H.hat <- predict(mm2)
SS1<-getInitial(TotDeaths~SSlogis(time,xmid,scale,alpha),data=Don)
#R_start1<-1/SS1["scale"]
#N0_start1<-SS1["alpha"]/(exp(SS1["xmid"]/SS1["scale"])+1)
R_start1 <- 0.048
N0_start1 <- 30
#fomula <-  TotDeaths~H.hat*N0*exp(R*time)/(H.hat+N0*(exp(R*time)-1))
fomula <-  TotDeaths~TotHosp*N0*exp(R*time)/(TotHosp+N0*(exp(R*time)-1))
m2<-nls(fomula,data=Don,start=list(R=R_start1,N0=N0_start1))
y.predi1 <- predict(m2)
#new.data1 <- data.frame(inHosp.Shifted = c(1660,1701), times=c(31,32))
#predict(m2, new.data1)
y1<- Don$TotDeaths
plot(times,y1,xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario")
lines(times,predict(m2),col="red",lty=2,lwd=3)
legend("topleft", legend=c("DeathObs", "DeathsPredi"),
       lty=c(1,1), pch=c(1,5), col=c("black", "red"))

pred<- data.frame(TotalDeaths=NewDat$TotalDeaths, PredDeaths=round(as.numeric(m2$m$fitted()),2))
xtable(pred)


plot(times,y1,xlab = "days",ylab="Total deaths", main = "Total deaths  in Ontario")
lines(times,predict(m2),col="red",lty=2,lwd=3)
lines(times,predict(m),col="blue",lty=2,lwd=3)
legend("topleft", legend=c("Kconstant", "Kvarying"),
       lty=c(1,1), pch=c(1,5), col=c("blue", "red"))

confint(m2)
confint2(mm2)

##---------------------------------------------------------------------------------------------------
###################################################################################################
####################################################################################################
library(deSolve)

Logistic = function(t, y, parms) {
  # Pull state variables from y vector
  P=y[1]
  # Pull parameter values from parms vector
  K = parms["K"]
  L = parms["L"]
  # Define equations
  dP = K * P * (1-P/L)
  res = dP
  # Return list of gradients
  list(res)
}

Gompertz = function(t, y, parms) {
  # Pull state variables from y vector
  P=y[1]
  # Pull parameter values from parms vector
  K = parms["K"]
  L = parms["L"]
  # Define equations
  dP = K * P * log(L/P)
  res = dP
  # Return list of gradients
  list(res)
}

Von = function(t, y, parms) {
  # Pull state variables from y vector
  P=y[1]
  # Pull parameter values from parms vector
  K = parms["K"]
  L = parms["L"]
  a = parms["a"]
  # Define equations
  dP = K * P * (((P/L)^a)-1)
  res = dP
  # Return list of gradients
  list(res)
}

times = seq(0, 10, by = 1)
parms = c( K = 0.5, L = 200)
parmsv = c( K = 0.5, L = 200,a=-0.2)
start = 5

outL=ode(y=start, times=times, func=Logistic, parms=parms)
outG=ode(y=start, times=times, func=Gompertz, parms=parms)

outV=ode(y=start, times=times, func=Von, parms=parmsv)

plot(x=outL[,1], y=outL[,2], ylab="cases", xlab= "Time", type="l",lwd = 3)
lines(x=outL[,1], y=outG[,2], col="red",lwd = 3)
lines(x=outL[,1], y=outV[,2], col="green",lwd = 3)
legend("topleft", legend=c("Logistic", "Gompertz","Von"),
       lty=c(1,1), pch=c(1,NA), col=c("black", "red","green"))

########################################################################################
###-----------------Ontario-----------------------------------------

ft <- function(t,N0,lam){
  N <- N0* exp(lam*t)
  return(N)
}

#N0 <- 1
N0 <- 2741
#lam_ON <- 0.219364
lam_ON <- 1.999163
t <- 1:30
P <- ft(t=t,N0=N0,lam=lam_ON)
round(P)

## ---------------------Age group-------------------------------------
ft <- function(t,N0,lam){
  N<-matrix(0,nrow = length(N0),ncol = length(t))
  for(i in 1:dim(N)[1])
    N[i,] <- N0[i]* exp(lam[i]*t)
  return(N)
}

t <- 1:13
Nf <- c(62,811,985,701,180)  # the last totat cases
#N0 <- 2*Nf*exp(-lam*tau)     # lam: growth rate, tau: doubling time

#lam <- c(0.1861434,0.332403,0.2317963,0.3932461)
lam <- c(0.196898,0.1173082,0.08002679,0.154756,0.1600582)
PP <- ft(t=t,N0=Nf,lam=lam)
round(PP)

##----------------------Deaths-----------------------------------
Dt <- function(t,N0,lam){
  N<-matrix(0,nrow = length(N0),ncol = length(t))
  for(i in 1:dim(N)[1])
    N[i,] <- N0[i]* exp(lam[i]*t)
  return(N)
}

t <- 1:13

#Df <- c(153,150,43,26,380,631,150,1381)  # the last totat cases
#N0 <- 2*Nf*exp(-lam*tau)     # lam: growth rate, tau: doubling time

Df <- c(381,4825,2592,9387)
#lam <- c(0.1861434,0.332403,0.2317963,0.3932461)
#lamd <- c(0.2259312,0.2426562,0.09203202,0.213014,0.1861434,0.332403,
#          0.2317963,0.3932461)
lamd <- c(0.2262472,0.2794076,0.2372876,0.3055486)
D <- Dt(t=t,N0=Df,lam=lamd)
round(D)

###--------------------------------------------------------------
logi <- function(t,N0,K,lam){
  N <- K/(1+((K/N0)-1)* exp(-lam*t))
  return(N)
}

t <- 1:9
lam=0.1514
K=50000
N0=8612
P <- logi(t=t,N0=N0,K=K,lam=lam)
round(P)

##---------Modelling the spread of infectious disease-----------------------##

# Here I.t and I0 are fractions of population
infect <- function(t,I0,lam){
  I.t <- 1/(1+(1/I0-1)* exp(-lam*t))
  return(I.t)
}

N0 <- 8612
Ptot <- 40000000
I0 <- N0/Ptot
lam=0.1514
Infc <- infect(t=t,I0=I0,lam=lam)
n.infct <- Infc*Ptot

############################################################################

##-------------------Growth curves

library(easynls)

dgro <- data.frame(t=timDC,df=dailyD$cumdeath)

# fit the data to exponential model
model1 <- nlsfit(data=dgro,model=6,start = c(a=50000,b=0.23))
#plot data
nlsplot(data=dgro,model=6,start = c(a=50000,b=0.23),xlab="Days",ylab="Deaths",position=1)

yC=4.429873*exp(0.1486281*35)
#fit the data to Gompertz model
model2 <- nlsfit(data=CC[-1,],model=10,start = c(a=30,b=0.4,c=-500))
nlsplot(data=dgro,model=10,start = c(a=500,b=15,c=0.25),xlab="Days",ylab="cases",position=1)

yG = 1529251*exp(-15.121*exp(-0.02*38))

#Fit  the data to logistic model
model3 <- nlsfit(data=dgro,model=7,start = c(a=0.95,b=0.20,c=0.25))
nlsplot(data=dgro,model=10,start = c(a=500,b=15,c=0.25),xlab="Days",ylab="cases",position=1)

y = 1529251*exp(-15.121*exp(-0.02*38))


##-------------Linealization---------------------
f <- function(x,a,b) {a * exp(b * x)}

start <- list(a=exp(coef(model.0)[1]), b=coef(model.0)[2], c=c.0)

fm0 <- nls(log(df) ~ log(a * (1 + b * (exp(-c * t)))^-1),data=dgro,start = start)

model <- nls(df ~ a * (1 + b * (exp(-c * t)))^-1, data = dgro, start = start)
model


#find the parameters for the equation
SS<-getInitial(df~SSlogis(t,alpha,xmid,scale),data=dgro)
y <- dailyD$cumdeath
#we used a different parametrization
K_start<-SS["alpha"]
R_start<-1/SS["scale"]
N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
#the formula for the model
log_formula<-formula(y~K*N0*exp(R*timDC)/(K+N0*(exp(R*timDC)-1)))
#fit the model
m<-nls(log_formula,start=list(K=K_start,R=R_start,N0=N0_start))
#estimated parameters
summary(m)
plot(timDC,y)
lines(timDC,predict(m),col="red",lty=2,lwd=3)

new.data <- data.frame(timDC = c(39,40))
predict(m, new.data)
##-------------------------
ff <- function(t){
  y <- exp(-1.25873+0.09644*t)
  return(y)
}
t<-0:85

ff(t)
##------------------------------------------------------------------------
P <- daily_ON$cum
dP <- rep(0,85)
for(i in 1:85){
  dP[i]<- P[i+1]-P[i]
}

pDp <- dP/P[1:85]
dta <- data.frame(dP,P1=P[1:85])
plot(P[1:85],pDp)
abline(lm(pDp~P[1:85]))

########################################################################################
library(growthcurver)
model.wt <- SummarizeGrowth(Don$time, Don$TotalDeaths)
model.wt$vals
predict(model.wt$model)
predict(model.wt, new.data)


p1 <- ggplot(Don, aes(x=times,y=TotalDeaths)) + geom_point(alpha=0.5) + theme_bw()
df.predicted <- data.frame(times = Don$time, pred.wt = predict(model.wt$model))
p1 + geom_line(data=df.predicted, aes(y=df.predicted$pred.wt), color="red")


##-----------------------------------------------
dt <- read.csv(file = "death.csv",header = T)


death.pred <- nls(TotalDeaths ~ SSlogis(time, K, N0, r), data = dt)

summary(death.pred)


##-------------Making predictions-----------------------
predict.death <- predict(death.pred, data.frame(time = c(31,38,44,100)))
predict.death
death.pred$m$fitted()

##------------Plotting data----------------------------

alpha <- coef(death.pred)  #extracting coefficients
plot(TotalDeaths ~ time, data = dt, main = "Logistic Growth Model of Deaths in Ontario", 
     xlab = "Days", ylab = "Deaths", xlim = c(0, 40), ylim = c(0, 1000))  # Census data
curve(alpha[1]/(1 + exp(-(x - alpha[2])/alpha[3])), add = T, col = "red")  # Fitted model
points(31, 788.3499 , pch = "x", cex = 1.3)
points(38, 876.9210, pch = "x", cex = 1.3)


