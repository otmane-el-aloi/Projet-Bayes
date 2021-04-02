knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
#install.packages("truncnorm")
#install.packages("invgamma")
library(invgamma)
library(truncnorm)
# library(acp)
X <- c(1, 1.5, 1.5, 1.5, 2.5, 4, 5, 5, 7, 8, 8.5, 9, 9.5, 9.5, 10, 
       12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29, 31.5)

Y <- c(1.8, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 
       2.26, 2.4, 2.39, 2.41, 2.5, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 
       2.47, 2.64, 2.56, 2.7, 2.72, 2.57)

data <- data.frame(X,Y)

ggplot(data,mapping = aes(X,Y))+
  geom_point()+
  geom_smooth()+
  ggtitle("Evolution de la longueur des Dugongs en fonction de l'age") +
  xlab("age") + ylab("Longeur")
y <- density(data$Y)
plot(y)

dugong <- function(age, Y, nchain = 10^4, prop_sd = c(0.01,0.1,1,0.01)){
  
  init <- c(2,1,11,0.9)  #alpha, beta (alpha et beta superieur à 1), tau (standard deviation donc on commence à 1) et gamma (entre 0 et 1 donc on prend 0.5)   
  chain <- matrix(NA, nchain + 1, 4)
  chain[1,] <- init
  acc_rates <- rep(0, 4)
  n<-length(Y)
  sig<-1000#ecart type
  for (iter in 1:nchain){
    current <- chain[iter,]
    
    ## Mise a jour de alpha
    
    
    prop <- current
    prop[1] <- rlnorm(1 ,log(current[1]),prop_sd[1])
    kernel.ratio<-(prop[1]/current[1])
    top<-sum(dnorm(Y,prop[1]-current[2]*current[4]^age,1/sqrt(current[3]),log=TRUE))+dnorm(prop[1],0,sig,log=TRUE)
    #Bottom 
    bottom <- sum(dnorm(Y,current[1]-current[2]*current[4]^age,1/sqrt(current[3]),log=TRUE))+dnorm(current[1],0,sig,log=TRUE)
    #Alpha
    acc_prob <- exp(top - bottom)*kernel.ratio
    
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[1] <- acc_rates[1] + 1
    }
    #pour beta
    prop <- current
    prop[2] <- rlnorm(1 ,log(current[2]),prop_sd[2])
    kernel.ratio<-(prop[2]/current[2])
    top<-sum(dnorm(Y,current[1]-prop[2]*current[4]^age,1/sqrt(current[3]),log=TRUE))+dnorm(prop[2],0,sig,log=TRUE)
    #Bottom 
    bottom <- sum(dnorm(Y,current[1]-current[2]*current[4]^age,1/sqrt(current[3]),log=TRUE))+dnorm(current[2],0,sig,log=TRUE)
    #Alpha
    acc_prob <- exp(top - bottom)*kernel.ratio
    
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[2] <- acc_rates[2] + 1
    }
    
    
    
    
    ## Mise à jour de tau   
    prop <- current
    prop[3] <- rlnorm(1 ,log(current[3]),prop_sd[3])
    
    kernel.ratio<-(prop[3]/current[3])
    
    #Top
    top <- sum(dnorm(Y,current[1]-current[2]*current[4]^age,1/sqrt(prop[3]),log=TRUE))+dgamma(prop[3],shape=0.001,rate=0.001,log=TRUE)
    bottom <-sum(dnorm(Y,current[1]-current[2]*current[4]^age,1/sqrt(current[3]),log=TRUE))+dgamma(current[3],shape=0.001,rate=0.001,log=TRUE)   
    
    acc_prob <- exp(top - bottom)*kernel.ratio
    
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[3] <- acc_rates[3] + 1
    }
    
    ### Mise à jour de gamma
    prop <- current
    prop[4] <- rlnorm(1, log(current[4]), prop_sd[4])
    
    kernel.ratio<-(prop[4]/current[4])
    top <-sum(dnorm(Y,current[1]-current[2]*prop[4]^age,1/sqrt(current[3]),log=TRUE))+dunif(prop[4],min=0.5,max=1,log=TRUE)
    bottom <- sum(dnorm(Y,current[1]-current[2]*current[4]^age,1/sqrt(current[3]),log=TRUE))+dunif(current[4],min=0.5,max=1,log=TRUE)
    #Alpha
    acc_prob <- exp(top - bottom)*kernel.ratio
    
    if ((runif(1) < acc_prob)&(prop[4]<1)){
      current <- prop
      acc_rates[4] <- acc_rates[4] + 1
    } 
    
    ## Sauvegardons le nouvel etat
    chain[iter+1,] <- current
  }
  
  return(list(chain = chain, acc_rates = acc_rates / nchain))
}


nchain=10000
CM_params <- dugong(X, Y, nchain )


#CM_params
## remove burning period
CM_params$chain <- CM_params$chain[-(1:1000),]
expression<-c('alpha','beta','tau','gamma')
#plot
par(mfrow = c(4, 2), mar = c(4, 5, 0.5, 0.5))
for (j in 1:4){
  plot(CM_params$chain[,j], type = "l", ylab = expression[j])
  plot(density(CM_params$chain[,j]), type = "l", main = "")}

CM_params$acc_rates

moy_U3 <- mean(logit(CM_params$chain[8000:9001,4]))
moy_alpha <- mean(CM_params$chain[8000:9001,1])
moy_beta <- mean(CM_params$chain[8000:9001,2])
moy_sigma <- mean(1/(sqrt(CM_params$chain[8000:9001,3])))
moy_gamma <- mean(CM_params$chain[8000:9001,4])
sd_U3 <- sqrt(mean(logit(CM_params$chain[8000:9001,4])^2)-mean(logit(CM_params$chain[8000:9001,4]))^2)
sd_alpha <- sqrt(mean(CM_params$chain[8000:9001,1]^2)-mean(CM_params$chain[8000:9001,1])^2)
sd_beta <- sqrt(mean(CM_params$chain[8000:9001,2]^2)-mean(CM_params$chain[8000:9001,2])^2)
sd_sigma <- sqrt(mean((1/(sqrt(CM_params$chain[8000:9001,3])))^2)-mean((1/(sqrt(CM_params$chain[8000:9001,3]))))^2)
sd_gamma <- sqrt(mean(CM_params$chain[8000:9001,4]^2)-mean(CM_params$chain[8000:9001,4])^2)
quantile(logit(CM_params$chain[8000:9001,4]), prob=c(0.025,0.975))
quantile(CM_params$chain[8000:9001,1], prob=c(0.025,0.975))
quantile(CM_params$chain[8000:9001,2], prob=c(0.025,0.975))
quantile(1/(sqrt(CM_params$chain[8000:9001,3])), prob=c(0.025,0.975))
quantile(CM_params$chain[8000:9001,4], prob=c(0.025,0.975))

