# random walk Metropolis algorithm

RMMA <- function(n=1000, sigmas = 1, x0=0, beta = 1.5){
  x<- matrix(0, nrow=1, ncol =n)
  k=0

  x[1,1]<- x0
  for(i in 2:n){
    y<- x[1,i-1]+rnorm(1)*sigmas
    rmma <- exp(sum(abs(x[1,i-1])^beta-abs(y)^beta))
               if(runif(1)<rmma){
                 x[1,i]<- y
                 k <- k +1
                 }
               else{
                 x[1,i]<-x[,i-1]
               }
  }
  return (list(X=x, nb = k))
}
xmatrix11<- RMMA()

# run a sample with 1000 sample size to figure the trace of random walk
plot(seq(1,1000,1),xmatrix11$X, type = 's')


# produce a table which presents accept rate, mean, variance and scales to see impact of different
# on accept rate, autocorrelation rate, mean and variance
testassgn2 <- matrix(0, nrow =17, ncol= 5)

tablegenerate<- function(x){
  testassgnf<- as.matrix(x)
  listsigma <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 15, 20, 25, 30,35, 40,45,50)
  testassgnf[,1]<-listsigma
  for (i in 1:length(listsigma)){
    xmatrixn <- RMMA(n=10000, sigmas = listsigma[i], x0= 0, beta = 1.5)
    testassgnf[i,2]<- acf(xmatrixn$X[1,], plot = F)[1]$acf
    testassgnf[i,3]<- xmatrixn$nb/10000
    testassgnf[i,4]<- mean(xmatrixn$X[1,])
    testassgnf[i,5]<- var(xmatrixn$X[1,])
  }
  return(testassgnf)
}

dataframe1<-tablegenerate(testassgn2)
testtable<- data.frame("Scales"= dataframe1[,1], "ACF lag1"=dataframe1[,2], "Accept rate"= dataframe1[,3], "Mean"=dataframe1[,4], "Var"=dataframe1[,5] )

xmatrixn001 <- RMMA(n=5000, sigmas = 0.01, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn001$X, type = 's', ylab = "X(t)", main = "scale = 0.01")
xmatrixn005 <- RMMA(n=5000, sigmas = 0.05, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn005$X, type = 's', ylab = "X(t)",main = "scale = 0.05")
xmatrixn1 <- RMMA(n=5000, sigmas = 0.1, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn1$X, type = 's', ylab = "X(t)",main = "scale = 0.1")
xmatrixn2 <- RMMA(n=5000, sigmas = 0.5, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn2$X, type = 's',ylab = "X(t)",main = "scale = 0.5")
xmatrixn3 <- RMMA(n=5000, sigmas = 1, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn3$X, type = 's',ylab = "X(t)",main = "scale = 1")
xmatrixn4 <- RMMA(n=5000, sigmas = 5, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn4$X, type = 's',ylab = "X(t)",main = "scale = 5")
xmatrixn5 <- RMMA(n=5000, sigmas = 10, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn5$X, type = 's',ylab = "X(t)",main = "scale = 10")
xmatrixn6 <- RMMA(n=5000, sigmas =15, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn6$X, type = 's',ylab = "X(t)",main = "scale = 15")
xmatrixn7 <- RMMA(n=5000, sigmas =20, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn7$X, type = 's',ylab = "X(t)",main = "scale = 20")
xmatrixn8 <- RMMA(n=5000, sigmas =25, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn8$X, type = 's',ylab = "X(t)",main = "scale = 25")
xmatrixn9 <- RMMA(n=5000, sigmas =50, x0= 0, beta = 1.5)
plot(seq(1,5000,1),xmatrixn9$X, type = 's',ylab = "X(t)",main = "scale = 50")



# RMMAmiu5 <- function(n=1000, sigmas = 1, x0=0, beta = 1.5){
#   x<- matrix(0, nrow=1, ncol =n)
#   k=0
#   u<- runif(n)
#   x[,1]<- x0
#   for(i in 2:n){
#     y<- x[,i-1]+mvrnorm(mu = 0, Sigma = sigmas)
#     rmma <- exp(sum(abs(x[,i]-0)^beta-abs(y-0)^beta))
#     if(u[i]<rmma){
#       x[,i]<- y
#       k <- k +1
#     }
#     else{
#       x[,i]<-x[,i-1]
#     }
#   }
#   return (list(X=x, nb = k))
# }

plotgenerate<- function(sigmas1){
s1<- sigmas1
xmatrixn1 <- RMMA(n=10000, sigmas = s1, x0= 0, beta = 1.5)
consummean<- rep(0,10000)
for (i in 1:10000){
  consummean[i]<- sum(xmatrixn1$X[1:i])/i
}
plot(seq(1,10000,1),consummean, type = 's')
}

plotgenerate(0.1)
plotgenerate(0.5)
plotgenerate(1)
plotgenerate(5)
plotgenerate(10)
plotgenerate(20)

xmatrix700<- RMMA(n=700, sigmas = 1, x0= 0, beta = 1.5)
xmatrix500<- RMMA(n=500, sigmas = 1, x0= 0, beta = 1.5)
xmatrix900<- RMMA(n=900, sigmas = 1, x0= 0, beta = 1.5)
xmatrix1000<- RMMA(n=1000, sigmas = 1, x0= 0, beta = 1.5)
xmatrix2000<- RMMA(n=2000, sigmas = 1, x0= 0, beta = 1.5)
xmatrix3000<- RMMA(n=3000, sigmas = 1, x0= 0, beta = 1.5)
xmatrix5000<- RMMA(n=5000, sigmas = 1, x0= 0, beta = 1.5)
xmatrix10000<- RMMA(n=10000, sigmas = 1, x0= 0, beta = 1.5)

plotgenerateconv<- function(sigmas1, nb){
  s1<- sigmas1
  nb1<-nb
  xmatrixn1 <- RMMA(n=nb1, sigmas = s1, x0= 0, beta = 1.5)
  consummean<- rep(0,nb1)
  for (i in 1:nb1){
    consummean[i]<- sum(xmatrixn1$X[1:i])/i
  }
  plot(seq(1,nb1,1),consummean, type = 's')
}

plotgenerateconv(1,1000)
plotgenerateconv(1,2000)
plotgenerateconv(1,3000)
plotgenerateconv(1,5000)
plotgenerateconv(1,10000)


# plotgenerateconv<- function(sigmas1, nb){
#   s1<- sigmas1
#   nb1<- nb
#   xmatrixn1 <- RMMA(n=nb1, sigmas = s1, x0= 0, beta = 1.5)
#   consummean<- rep(0,10000)
#   for (i in 1:10000){
#     consummean[i]<- sum(xmatrixn1$X[1:i])/i
#   }
#   plot(consummean, type = 'l')
# }
# 
# plotgenerateconv(1, 10000)

par(mfrow=c(3,2))
# riemannsum(xmatrix10000$X, 10000, 1)
plot(cumsum(xmatrix10000$X) / seq_along(xmatrix10000$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")
# riemannsum(xmatrix5000$X, 5000, 1)
plot(cumsum(xmatrix5000$X) / seq_along(xmatrix5000$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")   
# riemannsum(xmatrix3000$X, 3000, 1)
plot(cumsum(xmatrix3000$X) / seq_along(xmatrix3000$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")   
# riemannsum(xmatrix2000$X, 2000, 1)
plot(cumsum(xmatrix2000$X) / seq_along(xmatrix2000$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")   
# riemannsum(xmatrix1000$X, 1000, 1)
plot(cumsum(xmatrix1000$X) / seq_along(xmatrix1000$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")   
plot(cumsum(xmatrix700$X) / seq_along(xmatrix700$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")   
plot(cumsum(xmatrix500$X) / seq_along(xmatrix500$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")   





# riemannsum(xmatrix900$X, 900, 1)
# riemannsum(xmatrix700$X, 700, 1)
# l1<-riemancumsumlist(xmatrix10000$X, 10000, 1)
# l2<-riemancumsumlist(xmatrix5000$X, 5000, 1)
# l3<-riemancumsumlist(xmatrix3000$X, 3000, 1)
# l4<-riemancumsumlist(xmatrix2000$X, 2000, 1)
# l5<-riemancumsumlist(xmatrix1000$X, 1000, 1)

riemansumlist()

mean(xmatrix1000$X^2)
mean(xmatrix2000$X^2)
mean(xmatrix3000$X^2)
mean(xmatrix5000$X^2)
mean(xmatrix10000$X^2)

# random walk metropoli to generate X(t), then calculate the second moment of f
# and variance of estimator 
RMMAreturnx <- function(n=1000, sigmas = 1, x0=0, beta = 1.5){
  x<- matrix(0, nrow=1, ncol =n)
  k=0
  
  x[1,1]<- x0
  for(i in 2:n){
    y<- x[1,i-1]+rnorm(1)*sigmas
    rmma <- exp(sum(abs(x[1,i-1])^beta-abs(y)^beta))
    if(runif(1)<rmma){
      x[1,i]<- y
    }
    else{
      x[,i]<-x[,i-1]
    }
  }
  return (x)
}
y<- matrix(0, ncol=10000, nrow =100)
meanlist<- rep(0, 100)
meanlist2<- rep(0, 100)
varlist2<- rep(0, 100)
for (i in 1:100){
  y[i,]<-RMMAreturnx(n=10000, sigmas = 1, x0=0, beta = 1.5)
  meanlist[i]<- mean(y[i,]^2)
  varlist2[i]<- var(y[i,])
  meanlist2[i]<-mean(y[i,])
}
var(meanlist)
var(meanlist2)
mean(meanlist)
mean(varlist2)

plot(cumsum(xmatrixburnin1000$X) / seq_along(xmatrixburnin1000$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")


par(mfrow = c(1,2))
xmatrixburnin001<- RMMA(n=5000, sigmas = 0.01, x0= 10, beta = 1.5)
plot(seq(1,5000,1),xmatrixburnin001$X, type = 'l')
plot(cumsum(xmatrixburnin001$X) / seq_along(xmatrixburnin001$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")
xmatrixburnin010<- RMMA(n=5000, sigmas = 0.1, x0= 10, beta = 1.5)
plot(seq(1,5000,1),xmatrixburnin010$X, type = 'l')
plot(cumsum(xmatrixburnin010$X) / seq_along(xmatrixburnin010$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")
xmatrixburnin050<- RMMA(n=5000, sigmas = 0.5, x0= 10, beta = 1.5)
plot(seq(1,5000,1),xmatrixburnin050$X, type = 'l')
plot(cumsum(xmatrixburnin050$X) / seq_along(xmatrixburnin050$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")
xmatrixburnin1<- RMMA(n=5000, sigmas = 1, x0= 10, beta = 1.5)
plot(seq(1,5000,1),xmatrixburnin1$X, type = 'l')
plot(cumsum(xmatrixburnin1$X) / seq_along(xmatrixburnin1$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")
xmatrixburnin5<- RMMA(n=5000, sigmas = 5, x0= 10, beta = 1.5)
plot(seq(1,5000,1),xmatrixburnin5$X, type = 'l')
plot(cumsum(xmatrixburnin5$X) / seq_along(xmatrixburnin5$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")
xmatrixburnin50<- RMMA(n=5000, sigmas = 50, x0= 10, beta = 1.5)
plot(seq(1,5000,1),xmatrixburnin50$X, type = 'l')
plot(cumsum(xmatrixburnin50$X) / seq_along(xmatrixburnin50$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means starting at (0,0)', col="red")


totalscalmatrix <- matrix(0,ncol = 1000, nrow = 51)
scalelist<- seq(1,2,0.02)

for(i in 1:51){
  testinterface<-RMMA(n=1000, sigmas = 1, x0= 0, beta = scalelist[i])
  for (j in 1:1000){
    totalscalmatrix[i,j]<- testinterface$X[j]
  }
}

riemannsum<- function(x, n, beta1){
  x1<-x
  be<- beta1
  x1<-sort(x1)
  nb<- n
  sum<- rep(0,nb)
  for (i in 2:length(x1)){
    sum[i]<- (x1[i]-x1[i-1])*exp(-abs(x1[i])^be)*(1/1.8)
      # producefx1(x1[i], miu1 = -1, miu2 = 2, sigma1= 0.2, sigma2= 0.3)
  }
  return(sum(sum))
}
riemansumlist<- rep(0, 51)

for (i in 1:51){
  riemansumlist[i]<- riemannsum(totalscalmatrix[i,], 1000, scalelist[i])
}
  
varlist<- function(x){
  matrixX<- as.matrix(x)
  varl<- rep(0,51)
  for (i in 1:51){
    varl[i]<- var(matrixX[i,])
  }
  return(varl)
}

varlist(totalscalmatrix)

riemancumsumlist<- function(x, n, beta2){
  x1<-x
  be<- beta2
  x1<-sort(x1)
  nb<- n
  sum<- rep(0,nb)
  for (i in 2:length(x1)){
    sum[i]<- (x1[i]-x1[i-1])*exp(-abs(x1[i])^be)*(1/1.8)
  }
  return(cumsum(sum))
}

plot(riemancumsumlist(totalscalmatrix[1,],1000,1), type = 'l')
plot(riemancumsumlist(totalscalmatrix[51,],1000,2), type = 'l')
xx<-rep(0, 51)
for (i in 1:51){
  x<- riemancumsumlist(totalscalmatrix[i,],1000,1)
    xx[i]<- x[1000]
}

gibbsmatrix2<- matrix(0, ncol = 1000, nrow = 2)

# fx function (density)
fx2<- function(x){
  if (x <=-3 && x>=-5){
    i1x<- 1}
  else{i1x<-0}
  if(x>=3 && x<=5){
    i2x<-1}
  else{i2x<-0}
  fx<- 3/8*(i1x*(1-(x+4)^2)+i2x*(1-(x-4)^2))
  return(fx)
}

xlist2<- seq(-10,10, 0.01)
xlist3<- rep(0,2000)
for (i in 1:length(xlist2)){
  xlist3[i]<- fx2(xlist2[i])
}

  
plot(xlist3, type = 'l')

# produce x appear in range of (-5, -3) and (3,5) equally

runif2<- function(x1, x2, x3, x4){
  
  if (runif(1)<0.5){
    runif2x1<- runif(1, x1, x2)
  }
  else{runif2x1<-runif(1, x3,x4)}
}


gibbsmatrixgenerate<- function(){
gibbsmatrix2<- matrix(0, ncol = 1000, nrow = 2)
gibbsmatrix2[1,1]<- runif(1,0,3/8)
gibbsmatrix2[2,1]<- frunifform2()
  # runif(1, -5, -3)
for (i in (2:1000)){
  gibbsmatrix2[1,i]<- runif(1, 0, fx2(gibbsmatrix2[2,i-1]))
  groots1<- sort(
    rootSolve::uniroot.all(
      function(x){fx2(x)-gibbsmatrix2[1,i]},
      lower = -5,
      upper = -3
    )
  )
  groots2<- sort(
    rootSolve::uniroot.all(
      function(x){fx2(x)-gibbsmatrix2[1,i]},
      lower = 3,
      upper = 5
    )
  )
  gibbsmatrix2[2,i]<- runif2(groots1[1],groots1[2], groots2[1],groots2[2])
}
return(gibbsmatrix2)
}
par(mfrow =c(1,1))
xmar<- gibbsmatrixgenerate()
plot(xmar[2,], xmar[1,], xlab = "X", ylab = "U")
plot(xmar[2,], xmar[1,],xlab = "X", ylab = "U", type = 's')

fxlistgibbsmean<- rep(0,100)
fx2list<- rep(0,100)
for (i in 1:100){
  xmatrixtest<- gibbsmatrixgenerate()
  fxlistgibbsmean[i]<-mean(xmatrixtest[2,])
  
}
  estimatefxmean_gibbs<- mean(fxlistgibbsmean)
  estimatefxvar_gibbs<- var(fxlistgibbsmean)
  
  
  
  fxlistgibbsvar<- rep(0,100)
  fx2listvar<- rep(0, 100)
for (i in 1:100){
    xmatrixtest<- gibbsmatrixgenerate()
    
    fxlistgibbsvar[i]<- var(xmatrixtest[2,])
}




frunifform2<- function(){
  x1<- runif(1, -5, -3)
  x2<- runif(1, 3,5)
  xi <- runif(1)
  if (xi <0.5){
    return(x1)
  }
  if (xi> 0.5){
    return(x2)
  }
}

fxu<- function(x){
  x1<-x
  u1<-runif(1, 0, fx2(x1))
  return(rbind(x1, u1))
}
  

RWMPSS1 <- function(n=1000, scaless = 1){
  xmatrixss<- matrix(0, nrow=1, ncol =n)
  k=0
  clist<-c()
  xmatrixss[1,1]<- frunifform2()
  # xmatrixss[2,1]<-runif(1, 0, fx2(xmatrixss[1,1]))
  for(i in 2:n){
    y<- xmatrixss[1,i-1]+rnorm(1)*scaless
    rwmpss<- fx2(y)/fx2(xmatrixss[1,i-1])
    if(runif(1)<rwmpss){
      xmatrixss[1,i]<- y
      k <- k +1
    }
    else{
      xmatrixss[1,i]<-xmatrixss[1,i-1]
    }
    # xmatrixss[2,i]<- runif(1,0,fx2(xmatrixss[1,i-1]))
  }
  return (list(X=xmatrixss, nb = k))
}


xmatrixx22y<- rep(0,1000)



RWMPSS <- function(n=1000, a){
  scale= a
  fxy<- 0
  fxyy<-0
  xmatrixss<- matrix(0, nrow=2, ncol =n)
  k=0
  clist<-c()
  xmatrixss[1,1]<- frunifform2()
  xmatrixss[2,1]<-runif(1, 0, fx2(xmatrixss[1,1]))
  for(i in 2:n){
    y<-rnorm(1, )
    if (fx2(xmatrixss[1,i-1])>=xmatrixss[2, i-1]){
    clist<- append(clist, xmatrixss[1,i-1])
    groots1<- sort(
      rootSolve::uniroot.all(
        function(x){fx2(x)-xmatrixss[2, i-1]},
        lower = -5,
        upper = -3
      )
    )
    groots2<- sort(
      rootSolve::uniroot.all(
        function(x){fx2(x)-xmatrixss[2, i-1]},
        lower = 3,
        upper = 5
      )
    )
    fxy<- abs(groots1[1]-groots1[2])+abs(groots2[1]-groots2[2])
    }
    
    # if (fx2(y)>=xmatrixss[2,i-1]){
    #   fxy<- sum(clist)+fx2(y)
    # }
    if(fx2(y)>=xmatrixss[2, i-1]){
      ygroots1<- sort(
        rootSolve::uniroot.all(
          function(x){fx2(x)-xmatrixss[2, i-1]},
          lower = -5,
          upper = -3
        )
      )
      ygroots2<- sort(
        rootSolve::uniroot.all(
          function(x){fx2(x)-xmatrixss[2, i-1]},
          lower = 3,
          upper = 5
        )
      )
      fxyy<- abs(ygroots1[1]-ygroots1[2])+abs(ygroots2[1]-ygroots2[2])
    }
    # else{fxy<-sum(clist)}
    # rmma <- sum(clist)/fxy
      rmma <- fxy/fxyy
    

#     
#     # groots1<- sort(
#     #   rootSolve::uniroot.all(
#     #     function(x){fx2(x)-fx2(y)},
#     #     lower = -5,
#     #     upper = -3
#     #   )
#     # )
#     # groots2<- sort(
#     #   rootSolve::uniroot.all(
#     #     function(x){fx2(x)-xmatrixss[2,i-1]},
#     #     lower = -5,
#     #     upper = -3
#     #   )
#     # )
#     # if (length(groots1)==0){
#     #   fxuy<- 1}
#     # if (length(groots1)!=0){
#     #   fxuy<- 1/(3/4*(-15*groots1[2]-groots1[2]^3/3-4*groots1[2]^2)-3/4*(-15*groots1[1]-groots1[1]^3/3-4*groots1[1]^2))
#     # }
#     # if(length(groots2)==0){
#     #   fxux<- 1}
#     # if (length(groots2)!=0){
#     #   fxux<- 1/(3/4*(-15*groots2[2]-groots2[2]^3/3-4*groots2[2]^2)-3/4*(-15*groots2[1]-groots2[1]^3/3-4*groots2[1]^2))
#     # }
#     # rmma <- fxuy/fxux
    if(runif(1)<rmma){
      xmatrixss[1,i]<- y
      k <- k +1
    }
    else{
      xmatrixss[1,i]<-xmatrixss[1,i-1]
    }
    xmatrixss[2,i]<- runif(1,0,fx2(xmatrixss[1,i-1]))
    }
  return (list(X=xmatrixss, nb = k))
}

list1<- RWMPSS(1000,5)
plot(list1$X[1,], list1$X[2,])


RWMPSStest001<-RWMPSS1(1000, 0.01)
RWMPSStest002<-RWMPSS1(1000, 0.02)
RWMPSStest005<-RWMPSS1(1000, 0.05)
RWMPSStest010<-RWMPSS1(1000, 0.1)
RWMPSStest015<-RWMPSS1(1000, 0.15)
RWMPSStest020<-RWMPSS1(1000, 0.2)
RWMPSStest025<-RWMPSS1(1000, 0.25)
RWMPSStest030<-RWMPSS1(1000, 0.3)
RWMPSStest035<-RWMPSS1(1000, 0.35)
RWMPSStest050<-RWMPSS1(1000, 0.5)
RWMPSStest075<-RWMPSS1(1000, 0.75)
RWMPSStest1<-RWMPSS1(1000, 1)
RWMPSStest1_5<-RWMPSS1(1000, 1.5)
RWMPSStest2<-RWMPSS1(1000, 2)
RWMPSStest2_5<-RWMPSS1(1000, 2.5)
RWMPSStest3<-RWMPSS1(1000, 3)
RWMPSStest5<-RWMPSS1(1000, 5)
RWMPSStest10<-RWMPSS1(1000, 10)
RWMPSStest15<-RWMPSS1(1000, 15)
RWMPSStest20<-RWMPSS1(1000, 20)
RWMPSStest25<-RWMPSS1(1000, 25)
RWMPSStest30<-RWMPSS1(1000, 30)
RWMPSStest35<-RWMPSS1(1000, 35)
RWMPSStest40<-RWMPSS1(1000, 40)
RWMPSStest45<-RWMPSS1(1000, 45)
RWMPSStest50<-RWMPSS1(1000, 50)
par(mfrow = c(3,3))
plot(RWMPSStest001$X[1,], type = "l", main = "scale = 0.01")
plot(RWMPSStest002$X[1,], type = "l", main = "scale = 0.02")
plot(RWMPSStest005$X[1,], type = "l", main = "scale = 0.05")
plot(RWMPSStest010$X[1,], type = "l", main = "scale = 0.1")
plot(RWMPSStest015$X[1,], type = "l", main = "scale = 0.15")
plot(RWMPSStest020$X[1,], type = "l", main = "scale = 0.2")
plot(RWMPSStest025$X[1,], type = "l", main = "scale = 0.25")
plot(RWMPSStest030$X[1,], type = "l", main = "scale = 0.3")
plot(RWMPSStest035$X[1,], type = "l", main = "scale = 0.35")
plot(RWMPSStest050$X[1,], type = "l", main = "scale = 0.5")
plot(RWMPSStest075$X[1,], type = "l", main = "scale = 0.75")
plot(RWMPSStest1$X[1,], type = "l",main = "scale = 1")
plot(RWMPSStest5$X[1,], type = "l",main = "scale = 5")
plot(RWMPSStest10$X[1,], type = "l", main = "scale = 10")
plot(RWMPSStest15$X[1,], type = "l", main = "scale = 15")
plot(RWMPSStest25$X[1,], type = "l", main = "scale = 25")
plot(RWMPSStest50$X[1,], type = "l", main = "scale = 50")
par(mfrow = c(1,3))

RWMPSStest<-RWMPSS(1000, 10)

scalelist2<- c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
fxlistMTPSS<- rep(0, 1000)
meanfxMTPSS<- rep(0, length(scalelist2))
for (i in 1:26){
  RWMPSSlist<- RWMPSS1(1000, scalelist2[i])
  
  meanfxMTPSS[i]<- mean(RWMPSSlist$X[1,])
}



fxlistMTPSSvar<- rep(0, 1000)
varfxMTPSS<- rep(0, length(scalelist2))
for (i in 1:26){
  RWMPSSlist<- RWMPSS1(1000, scalelist2[i])
  
  varfxMTPSS[i]<- var(RWMPSSlist$X[1,])
}
  
  matrixblank <- matrix(0, ncol = 5, nrow =26 )
tablegenerate2<- function(x){
  testassgnf<- as.matrix(x)
  listsigma <- c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  testassgnf[,1]<-listsigma
  for (i in 1:length(listsigma)){
    xmatrixn <- RWMPSS1(n=1000, listsigma[i])
    testassgnf[i,2]<- acf(xmatrixn$X[1,], plot = F)[1]$acf
    testassgnf[i,3]<- xmatrixn$nb/1000
    testassgnf[i,4]<- mean(xmatrixn$X[1,])
    testassgnf[i,5]<- var(xmatrixn$X[1,])
  }
  return(testassgnf)
}
tablesacaletotal<- tablegenerate2(matrixblank)
  scaletable2<- data.frame("Scales"= tablesacaletotal[,1], "ACF lag1"=tablesacaletotal[,2], "Accept rate"= tablesacaletotal[,3], "Mean"=tablesacaletotal[,4], "Var"=tablesacaletotal[,5] )
  

  par(mfrow= c(3,2))
  xmatrixtest<- gibbsmatrixgenerate()  
  plot(xmatrixtest[2,], type = 's')
  plot(RWMPSStest1$X[1,], type = "l",main = "scale = 1")
  plot(RWMPSStest5$X[1,], type = "l",main = "scale = 5")
  plot(cumsum(xmatrixtest[2,]) / seq_along(xmatrixtest[2,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means of slice sampler', col="red")
  plot(cumsum(RWMPSStest1$X) / seq_along(RWMPSStest1$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means of metropolis with scale 1', col="red")
  plot(cumsum(RWMPSStest5$X) / seq_along(RWMPSStest5$X), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means of metropolis with scale 5', col="red")
  
  for (i in 1:)
  
  for (i in 1:length(xmatrixx22$X)){
    xmatrixx22y[i]<- fx2(xmatrixx22$X[i])
  }
  
  plot(xmatrixx22$X, xmatrixx22y)
  
  # generate multivariate metropoli random walk
  RMMAmulti <- function(n=1000, sigmas = 1, x0=as.vector(c(0,1))){
    x<- matrix(0, nrow=2, ncol =n)
    k=0
    
    x[,1]<- x0
    for(i in 2:n){
      y<- x[,i-1]+mvrnorm(1, mu = c(0,0), Sigma = sigmas*diag(2))
      
      rmma <- exp(sum(abs(x[1,i-1])^x[2,i-1]-abs(y[1])^y[2]))
      if(runif(1)<rmma && (y[2]>1 && y[2]<2)){
        x[,i]<- y
        k <- k +1
      }
      else{
        x[,i]<-x[,i-1]
      }
      
      
    }
    return (list(X=x, nb = k))
  }
  
  xmatrixmulti001<- RMMAmulti(1000, 0.01)
  xmatrixmulti01<- RMMAmulti(1000, 0.1)
  xmatrixmulti05<- RMMAmulti(1000, 0.5)
  xmatrixmulti1<- RMMAmulti(1000, 1)
  xmatrixmulti2<- RMMAmulti(1000, 2)
  xmatrixmulti5<- RMMAmulti(1000, 5)
  xmatrixmulti10<- RMMAmulti(1000, 10)
  xmatrixmulti15<- RMMAmulti(1000, 15)
  xmatrixmulti20<- RMMAmulti(1000, 20)
 
  
  plot(xmatrixmulti001$X[1,], xmatrixmulti001$X[2,], type = "s", main = "scale = 0.01", ylab = "beta", xlab = "x")
  plot(xmatrixmulti01$X[1,], xmatrixmulti01$X[2,], type = "s", main = "scale = 0.1", ylab = "beta", xlab = "x")
  plot(xmatrixmulti05$X[1,], xmatrixmulti05$X[2,], type = "s", main = "scale = 0.5", ylab = "beta", xlab = "x")
  plot(xmatrixmulti1$X[1,], xmatrixmulti1$X[2,], type = "s", main = "scale = 1", ylab = "beta", xlab = "x")
  plot(xmatrixmulti2$X[1,], xmatrixmulti2$X[2,], type = "s", main = "scale = 2", ylab = "beta", xlab = "x")
  plot(xmatrixmulti5$X[1,], xmatrixmulti5$X[2,], type = "s", main = "scale = 5", ylab = "beta", xlab = "x")
  plot(xmatrixmulti10$X[1,], xmatrixmulti10$X[2,], type = "s", main = "scale = 10", ylab = "beta", xlab = "x")
  
  RMMAmulti1 <- function(n=1000, sigmas = 1, x0=as.vector(c(0,1))){
    x<- matrix(0, nrow=2, ncol =n)
    k=0
    
    x[,1]<- x0
    for(i in 2:n){
      y<- x[,i-1]+mvrnorm(1, mu = 0, Sigma = sigmas)
      
      rmma <- exp(sum(abs(x[1,i-1])^x[2,i-1]-abs(y[1])^y[2]))
      if(runif(1)<rmma && (y[2]>1 && y[2]<2)){
        x[,i]<- y
        k <- k +1
      }
      else{
        x[,i]<-x[,i-1]
      }
      
      
    }
    return (list(X=x, nb = k))
  }
  
  matrixblankmetro<- matrix(0, nrow = 15, ncol = 5)
  tablegenerate3<- function(x){
    testassgnf<- as.matrix(x)
    listsigma <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2, 3, 5, 10, 15, 20,25, 30)
    testassgnf[,1]<-listsigma
    for (i in 1:length(listsigma)){
      xmatrixn <- RMMAmulti(n=1000, sigmas = listsigma[i], x0= as.vector(c(0,1)))
      testassgnf[i,2]<- acf(xmatrixn$X[1,], plot = F)[1]$acf
      testassgnf[i,3]<- xmatrixn$nb/1000
      testassgnf[i,4]<- mean(xmatrixn$X[1,])
      testassgnf[i,5]<- var(xmatrixn$X[1,])
    }
    return(testassgnf)
  }
  dataframe2<- tablegenerate3(matrixblankmetro)
  multitesttable<- data.frame("Scales"= dataframe2[,1], "ACF lag1"=dataframe2[,2], "Accept rate"= dataframe2[,3], "Mean"=dataframe2[,4], "Var"=dataframe2[,5] )
  
  par(mfrow = c(3,3))
  plot(xmatrixmulti001$X[1,], type = "s", main = "scale = 0.01", ylab = "x")
  plot(xmatrixmulti01$X[1,], type = "s", main = "scale = 0.1", ylab = "x")
  plot(xmatrixmulti05$X[1,], type = "s", main = "scale = 0.05", ylab = "x")
  plot(xmatrixmulti1$X[1,], type = "s", main = "scale = 1", ylab = "x")
  plot(xmatrixmulti2$X[1,], type = "s", main = "scale = 2", ylab = "x")
  plot(xmatrixmulti5$X[1,], type = "s", main = "scale = 5", ylab = "x")
  plot(xmatrixmulti10$X[1,], type = "s", main = "scale = 10", ylab = "x")
  plot(xmatrixmulti001$X[2,], type = "s", main = "scale = 0.01", ylab = "beta")
  plot(xmatrixmulti01$X[2,], type = "s", main = "scale = 0.1", ylab = "beta")
  plot(xmatrixmulti05$X[2,], type = "s", main = "scale = 0.5", ylab = "beta")
  plot(xmatrixmulti1$X[2,], type = "s", main = "scale = 1", ylab = "beta")
  plot(xmatrixmulti2$X[2,], type = "s", main = "scale = 2", ylab = "beta")
  plot(xmatrixmulti5$X[2,], type = "s", main = "scale = 5", ylab = "beta")
  plot(xmatrixmulti10$X[2,], type = "s", main = "scale = 10", ylab = "beta")
  
  
  plot(cumsum(xmatrixmulti001$X[1,]) / seq_along(xmatrixmulti001$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 0.01', col="red")
  plot(cumsum(xmatrixmulti01$X[1,]) / seq_along(xmatrixmulti01$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 0.1', col="red")
  plot(cumsum(xmatrixmulti05$X[1,]) / seq_along(xmatrixmulti05$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 0.5', col="red")
  plot(cumsum(xmatrixmulti1$X[1,]) / seq_along(xmatrixmulti1$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 1', col="red")
  plot(cumsum(xmatrixmulti2$X[1,]) / seq_along(xmatrixmulti2$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 2', col="red")
  plot(cumsum(xmatrixmulti5$X[1,]) / seq_along(xmatrixmulti5$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 5', col="red")
  plot(cumsum(xmatrixmulti10$X[1,]) / seq_along(xmatrixmulti10$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 10', col="red")
  plot(cumsum(xmatrixmulti15$X[1,]) / seq_along(xmatrixmulti15$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 15', col="red")
  plot(cumsum(xmatrixmulti20$X[1,]) / seq_along(xmatrixmulti20$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 20', col="red")
  plot(cumsum(xmatrixmulti10$X[1,]) / seq_along(xmatrixmulti10$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 10', col="red")
  
  
  
  RWMPSS111<- function(n= 1000, scaless = 1){
    xmatrixss <- matrix(0, nrow = 2, ncol= n)
    k =0
    xmatrixss[1,1]<- frunifform2()
    xmatrixss[2,1]<- runif(1,0, fx2(xmatrixss[1,1]))
    for (i in 2:n){
      y<- xmatrixss[1,i-1]+rnorm(1)*scaless
      groots1<- sort(
        rootSolve::uniroot.all(
          function(x){fx2(y)-},
          lower = -5,
          upper = -3
        )
      )
      groots2<- sort(
        rootSolve::uniroot.all(
          function(x){fx2(x)-xmatrixss[2,i-1]},
          lower = -5,
          upper = -3
        )
      )
      if (length(groots1)==0){
        fxuy<- 1}
      if (length(groots1)!=0){
        fxuy<- 1/(3/4*(-15*groots1[2]-groots1[2]^3/3-4*groots1[2]^2)-3/4*(-15*groots1[1]-groots1[1]^3/3-4*groots1[1]^2))
      }
      if(length(groots2)==0){
        fxux<- 1
      }
      if (length(groots2)!=0){
        fxux<- 1/(3/4*(-15*groots2[2]-groots2[2]^3/3-4*groots2[2]^2)-3/4*(-15*groots2[1]-groots2[1]^3/3-4*groots2[1]^2))
      }
      rmma <- fxuy/fxux
      if(runif(1)<rmma){
        xmatrixss[,i]<- y
        k <- k +1
      }
      else{
        xmatrixss[,i]<-xmatrixss[,i-1]
      }
        
      
      
      RWMPSS1111 <- function(n=1000, scaless = 1){
        xmatrixss<- matrix(0, nrow=2, ncol =n)
        k=0
        clist<-c()
        xmatrixss[1,1]<- frunifform2()
        xmatrixss[2,1]<- runif(1, 0, fx2(xmatrixss[1,1]))
        # xmatrixss[2,1]<-runif(1, 0, fx2(xmatrixss[1,1]))
        for(i in 2:n){
          y<- xmatrixss[1,i-1]+rnorm(1)*scaless
          rwmpss<- fx2(y)/fx2(xmatrixss[1,i-1])
          if(runif(1)<rwmpss){
            xmatrixss[1,i]<- y
            k <- k +1
          }
          else{
            xmatrixss[1,i]<-xmatrixss[1,i-1]
          }
          xmatrixss[2,i]<- runif(1,0,fx2(xmatrixss[1,i-1]))
        }
        return (list(X=xmatrixss, nb = k))
      }
      
      metropoliceslicesampler1<- RWMPSS1111(1000, 0.1)
      metropoliceslicesampler2<- RWMPSS1111(1000, 0.5)
      metropoliceslicesampler3<- RWMPSS1111(1000, 1)
      metropoliceslicesampler4<- RWMPSS1111(1000, 2)
      metropoliceslicesampler5<- RWMPSS1111(1000, 3)
      metropoliceslicesampler6<- RWMPSS1111(1000, 5)
      par(mfrow =c(2,3))
      plot(metropoliceslicesampler1$X[1, ], metropoliceslicesampler1$X[2, ], ylim = c(0, 0.375), main ="scale =0.1")
      plot(metropoliceslicesampler2$X[1, ], metropoliceslicesampler2$X[2, ], ylim = c(0, 0.375), main ="scale =0.5")
      plot(metropoliceslicesampler3$X[1, ], metropoliceslicesampler3$X[2, ], ylim = c(0, 0.375), main ="scale =1")
      plot(metropoliceslicesampler4$X[1, ], metropoliceslicesampler4$X[2, ], ylim = c(0, 0.375), main ="scale =2")
      plot(metropoliceslicesampler5$X[1, ], metropoliceslicesampler5$X[2, ], ylim = c(0, 0.375), main ="scale =3")
      plot(metropoliceslicesampler6$X[1, ], metropoliceslicesampler6$X[2, ], ylim = c(0, 0.375), main ="scale =5")
      
      plot(cumsum(metropoliceslicesampler1$X[1,]) / seq_along(metropoliceslicesampler1$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 0.01', col="red")
      
      plot(cumsum(metropoliceslicesampler2$X[1,]) / seq_along(metropoliceslicesampler2$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 0.01', col="red")
      
      
      plot(cumsum(metropoliceslicesampler3$X[1,]) / seq_along(metropoliceslicesampler3$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 0.01', col="red")
      plot(cumsum(metropoliceslicesampler4$X[1,]) / seq_along(metropoliceslicesampler4$X[1,]), type = 'l', xlab = 'Iteration', ylab = 'Running Averaae',main = 'Trace of empirical means under scale 0.01', col="red")
      
      
  