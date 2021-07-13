########################################################################################
### MEMA: Measurement Error in Meta-Analysis
### contact: harlan.campbell@stat.ubc.ca

#####################
# Fresh start 

rm(list = ls(all = TRUE))
set.seed(87654321)




#####################
# Determine missing packages and load them:
required_packages <- c("MCMCvis", "latex2exp", "Matrix", "RCurl", "runjags", "rjags", "sjstats", "rlist", "tidyverse", "xtable", "MCMCpack")
not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed, type="source")                                           
invisible(suppressWarnings(lapply(required_packages, require, character.only = TRUE)))
ls()
rm("not_installed", "required_packages")
ls()


list_to_array<- function(mylist, dim1, dim2, dim3){
		aperm((array(as.numeric(unlist(mylist)), dim=c(dim1,dim2,dim3))))}


# Initial values for MCMC random seeds

	inits1<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(123))
	inits2<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(456))
	inits3<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(789))


# Data:



schooldata <- read.csv("https://raw.githubusercontent.com/harlanhappydog/MEMA/master/13schools.csv")

########## original schooldata with X = [reading, math]  (ADNELS88) 
### dataset :  NELS88_readmath

XtXinv <- lambda_hat <- mu_hat <- sigma_hat <- beta_hat  <- n <- X  <- COV <- y <- list()
nvec <- vector()
K <- length(table(schooldata$sch))
for(i in 1:K){
	nvec[i] <- dim(schooldata[schooldata$sch==i,])[1]	
	mydat 	<- schooldata[schooldata$sch==i,]
	n[[i]] 	<- dim(mydat)[1]
	y[[i]] 	<- mydat$sci
	X[[i]] 	<- cbind(mydat$rdg, mydat$math)
	mod_i 	<- lm(y[[i]] ~ X[[i]])
	beta_hat[[i]] 	<- coef(mod_i)
	sigma_hat[[i]] 	<- summary(mod_i)$sigma
	COV[[i]] 		<- solve(t(cbind(1,X[[i]]))%*%
						cbind(1,X[[i]]))*(sigma_hat[[i]]^2)
	XtXinv[[i]] 		<- solve(t(cbind(1,X[[i]]))%*%cbind(1,X[[i]]))
	lambda_hat[[i]] <- matrix(var(X[[i]]), dim(X[[i]])[2], dim(X[[i]])[2])	
	mu_hat[[i]]		<- colMeans(X[[i]])
	}

# all summary data in a list:
NELS88_readmath <-list()
NELS88_readmath$N 		<- sum(unlist(n))
NELS88_readmath$studyID <- schooldata$sch
NELS88_readmath$Q 		<- dim(X[[1]])[2]
NELS88_readmath$K 		<- length(table(schooldata$sch))
NELS88_readmath$y 		<- y
NELS88_readmath$X 		<- X
NELS88_readmath$nvec 	<- nvec
NELS88_readmath$beta_hat <- beta_hat
NELS88_readmath$sigma_hat <- sigma_hat
NELS88_readmath$COV <- COV
NELS88_readmath$XtXinv <- XtXinv
NELS88_readmath$lambda_hat <- lambda_hat
NELS88_readmath$mu_hat <- mu_hat





# cleaning up:
rm(mydat, mod_i, XtXinv, lambda_hat, mu_hat, sigma_hat, beta_hat, n, X, COV, y, nvec, K, i)
ls()

########



########## original schooldata with X = [reading, math]  (NELS88) 
### dataset :  NELS88star_readmath


schooldata <- read.csv("https://raw.githubusercontent.com/harlanhappydog/MEMA/master/13schools.csv")

kprime <- 5

XtXinv_me <- lambda_hat_me <- mu_hat_me <- sigma_hat_me <- beta_hat_me  <- n_me <- X_me  <- COV_me <- y_me <- list()
nvec_me <- vector()


schooldata_me <- schooldata

K <- length(table(schooldata_me$sch))

set.seed(123)

PHI_rdg <- rep(0,K)
PHI_rdg[1:kprime] <- 0
if(kprime<K){  PHI_rdg[(kprime+1):K] <- sample(c(4,5,6),length((kprime+1):K), replace=TRUE)  }

PHI_math <- rep(0,K)
PHI_math[1:kprime] <- 0
if(kprime<K){  PHI_math[(kprime+1):K] <- sample(c(8,10,12),length((kprime+1):K), replace=TRUE)  }

schooldata_me$math <- schooldata$math
for(i in (kprime+1):K){
schooldata_me[schooldata_me$sch==i,]$math <- schooldata[schooldata$sch==i,]$math +
				rnorm(length(schooldata[schooldata$sch==i,]$math), 0, PHI_math[i])
				}


schooldata_me$rdg <- schooldata$rdg
for(i in (kprime+1):K){
schooldata_me[schooldata_me$sch==i,]$rdg <- schooldata[schooldata$sch==i,]$rdg +
				rnorm(length(schooldata[schooldata$sch==i,]$rdg), 0, PHI_rdg[i])
				}


for(i in 1:K){
	nvec_me[i] <- dim(schooldata_me[schooldata_me$sch==i,])[1]	
	mydat_me 	<- schooldata_me[schooldata_me$sch==i,]
	n_me[[i]] 	<- dim(mydat_me)[1]
	y_me[[i]] 	<- mydat_me$sci
	X_me[[i]] 	<- cbind(mydat_me$rdg, mydat_me$math)
	mod_i_me 	<- lm(y_me[[i]] ~ X_me[[i]])
	beta_hat_me[[i]] 	<- coef(mod_i_me)
	sigma_hat_me[[i]] 	<- summary(mod_i_me)$sigma
	COV_me[[i]] 		<- solve(t(cbind(1,X_me[[i]]))%*%
						cbind(1,X_me[[i]]))*(sigma_hat_me[[i]]^2)
	XtXinv_me[[i]] 		<- solve(t(cbind(1,X_me[[i]]))%*%cbind(1,X_me[[i]]))
	lambda_hat_me[[i]] <- matrix(var(X_me[[i]]), dim(X_me[[i]])[2], dim(X_me[[i]])[2])	
	mu_hat_me[[i]]		<- colMeans(X_me[[i]])
	}



schooldata <- read.csv("https://raw.githubusercontent.com/harlanhappydog/MEMA/master/13schools.csv")

########## original schooldata with X = [reading, math]  (ADNELS88) 
### dataset :  NELS88_readmath

XtXinv <- lambda_hat <- mu_hat <- sigma_hat <- beta_hat  <- n <- X  <- COV <- y <- list()
nvec <- vector()
K <- length(table(schooldata$sch))
for(i in 1:K){
	nvec[i] <- dim(schooldata[schooldata$sch==i,])[1]	
	mydat 	<- schooldata[schooldata$sch==i,]
	n[[i]] 	<- dim(mydat)[1]
	y[[i]] 	<- mydat$sci
	X[[i]] 	<- cbind(mydat$rdg, mydat$math)}


cov(X_me[[1]]); cov(X_me[[8]]); cov(X_me[[13]])
cov(X[[1]]); cov(X[[8]]); cov(X[[13]])


# all summary data in a list:
NELS88star_readmath 		<- list()
NELS88star_readmath$N 		<- sum(unlist(n_me))
NELS88star_readmath$studyID   <- schooldata$sch
NELS88star_readmath$Q 		<- dim(X_me[[1]])[2]
NELS88star_readmath$K 		<- length(table(schooldata$sch))
NELS88star_readmath$y 		<- y_me
NELS88star_readmath$X 		<- X_me
NELS88star_readmath$nvec 	<- nvec_me
NELS88star_readmath$beta_hat 	<- beta_hat_me
NELS88star_readmath$sigma_hat <- sigma_hat_me
NELS88star_readmath$COV 		<- COV_me
NELS88star_readmath$XtXinv 	<- XtXinv_me
NELS88star_readmath$lambda_hat<- lambda_hat_me
NELS88star_readmath$mu_hat 	<- mu_hat_me


TableA<-data.frame(cbind(999, NELS88_readmath$nvec,cbind(do.call(rbind,NELS88_readmath$beta_hat), do.call(rbind,NELS88star_readmath$beta_hat))[,c(1,4,2,5,3,6)], PHI_rdg, PHI_math))
xtable(TableA)


# cleaning up:
rm(mydat_me, mod_i_me, XtXinv_me, lambda_hat_me, mu_hat_me, sigma_hat_me, beta_hat_me, n_me, X_me, COV_me, y_me, nvec_me, K, i)
ls()

#### #### #### #### #### #### #### #### #### #### #### #### 


library(rjags)

multiMA <- "model{

for(q in 1:(Q+1)){    
    tau[q]      ~   dt(0, 2.5, 1)I(0,);
	theta[q]    ~   dnorm(0, 1/100);  }

for (k in 1:NStudies){	  
for(q in 1:(Q+1)){   beta[k,q] ~ dnorm(theta[q], 1/(tau[q]^2));  }
  sigma[k] ~ dt(0, 2.5, 1)I(0,);
  mu[k,1:Q]  ~   dmnorm(m0, inverse(D100));
  invLambda[k,1:Q,1:Q] ~ dwish(D1[1:Q,1:Q], Q);
  Lambda[k,1:Q,1:Q] <- inverse(invLambda[k,1:Q,1:Q]);
}


for (j in 1:N){
X[j,2:(Q+1)] ~ dmnorm(mu[studyID[j],1:Q],  invLambda[studyID[j],1:Q,1:Q]);	
Y[j] ~ dnorm(inprod(X[j,], beta[studyID[j],]), 1/(sigma[studyID[j]]^2));
}

}"



multiMEMA <- "model{

for (k in 1:NStudies){
invLambda[k,1:Q,1:Q] ~ dwish(D1[1:Q,1:Q], Q+1);
Lambda[k,1:Q,1:Q] <- inverse(invLambda[k,1:Q,1:Q]);
invPhi[k,1:Q,1:Q] ~ dwish(2*zeta2*D1[1:Q,1:Q], 2*zeta1);
Phi[k,1:Q,1:Q] <- inverse(invPhi[k,1:Q,1:Q]);
}

  zeta1 ~ dexp(rho)I(Q/2,);
  zeta2 ~ dexp(rho);

for(q in 1:(Q+1)){    
    tau[q]      ~   dt(0, 2.5, 1)I(0,);
	theta[q]    ~   dnorm(0, 1/100);
	}



for (k in 1:NStudies){	  

for(q in 1:(Q+1)){   beta[k,q] ~ dnorm(theta[q], 1/(tau[q]^2));  }  
  sigma[k] ~ dt(0, 2.5, 1)I(0,);
  mu[k,1:Q]  ~   dmnorm(m0, inverse(D100));
}


for (j in 1:Nprime){
X[j,1] <- 1;
X[j,2:(Q+1)] ~ dmnorm(mu[studyID[j],1:Q], invLambda[studyID[j],1:Q,1:Q]);
Xstar[j,2:(Q+1)] ~ dmnorm(X[j,2:(Q+1)], inverse(zeromat));
Y[j] ~ dnorm(inprod(X[j,], beta[studyID[j],]), 1/(sigma[studyID[j]]^2));
}

for (j in (Nprime+1):N){
X[j,1] <- 1;
X[j,2:(Q+1)] ~ dmnorm(mu[studyID[j],1:Q], invLambda[studyID[j],1:Q,1:Q]);
Xstar[j,2:(Q+1)] ~ dmnorm(X[j,2:(Q+1)], invPhi[studyID[j],1:Q,1:Q]);
Y[j] ~ dnorm(inprod(X[j,], beta[studyID[j],]), 1/(sigma[studyID[j]]^2));
}

}"
                            
  
#######

# line 0:

QQ<-2
K<-13
N<-length(unlist(NELS88star_readmath$y))

jags.m <- jags.model(textConnection(multiMA), data = list (
  N = length(unlist(NELS88star_readmath$y)),
  Q = QQ,
  NStudies = K,
  studyID = NELS88star_readmath$studyID,
  Y = unlist(NELS88star_readmath$y),
  X = do.call(rbind,lapply(NELS88star_readmath$X, function(x) {(cbind(1,x))}))[,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ)
  ),
      inits =  list(inits1, inits2, inits3),  
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps1 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps1)
line1_multi <- summary(samps1)$quantiles[,c(1,3,5)]




#######

# line 1:

QQ<-2
kprime<-13
K<-13
N<-length(unlist(NELS88star_readmath$y))
Nprime <- suppressWarnings(min(N,min(c(1:N)[NELS88star_readmath$studyID>kprime])-1))
jags.m <- jags.model(textConnection(multiMEMA), data = list (
  N = length(unlist(NELS88star_readmath$y)),
  Nprime = Nprime,
  Q = QQ,
  NStudies = K,
  studyID = NELS88star_readmath$studyID,
  Y = unlist(NELS88star_readmath$y),
  Xstar = do.call(rbind,lapply(NELS88star_readmath$X, function(x) {(cbind(1,x))}))[,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ),
  zeromat=as.array(matrix(diag(rep(1,QQ)),nrow=QQ)/(1e+20)),
  rho=0.1
  ),
      inits =  list(inits1, inits2, inits3),  
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps1 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps1)
line1_multi <- summary(samps1)$quantiles[,c(1,3,5)]

#######

# line 2:

QQ<-2
kprime<-13
K<-13
N<-length(unlist(NELS88_readmath$y))
Nprime <- suppressWarnings(min(N,min(c(1:N)[NELS88_readmath$studyID>kprime])-1))
jags.m <- jags.model(textConnection(multiMEMA), data = list (
  N = length(unlist(NELS88_readmath$y)),
  Nprime = Nprime,
  Q = QQ,
  NStudies = K,
  studyID = NELS88_readmath$studyID,
  Y = unlist(NELS88_readmath$y),
  Xstar = do.call(rbind,lapply(NELS88_readmath$X, function(x) {(cbind(1,x))}))[,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ),
  zeromat=as.array(matrix(diag(rep(1,QQ)),nrow=QQ)/(1e+20)),
  rho=0.1
  ),
        inits =  list(inits1, inits2, inits3),  
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps2 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps2)
line2_multi <- summary(samps2)$quantiles[,c(1,3,5)]



#######

# line 3:

QQ<-2
kprime<-5
K<-13
N<-length(unlist(NELS88star_readmath$y))
Nprime <- suppressWarnings(min(N,min(c(1:N)[NELS88star_readmath$studyID>kprime])-1))
jags.m <- jags.model(textConnection(multiMEMA), data = list (
  N = length(unlist(NELS88star_readmath$y)),
  Nprime = Nprime,
  Q = QQ,
  NStudies = K,
  studyID = NELS88star_readmath$studyID,
  Y = unlist(NELS88star_readmath$y),
  Xstar = do.call(rbind,lapply(NELS88star_readmath$X, function(x) {(cbind(1,x))}))[,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ),
  zeromat=as.array(matrix(diag(rep(1,QQ)),nrow=QQ)/(1e+20)),
  rho=0.1
  ),
          inits =  list(inits1, inits2, inits3),  
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps3 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps3)
line3_multi <- summary(samps3)$quantiles[,c(1,3,5)]
line3_multi

MCMCtrace(samps3, params = c("theta"), 
 priors = cbind(rnorm(100000, 0, 10), rnorm(100000, 0, 10), rnorm(100000, 0, 10)) , main_den = c(			   
  TeX("Density $\\theta_{1}$"),
  TeX("Density $\\theta_{2}$"),
    TeX("Density $\\theta_{3}$")),
  main_tr = c(					   					                         
    TeX("Trace $\\theta_{1}$"),
    TeX("Trace $\\theta_{2}$"),
    TeX("Trace $\\theta_{3}$")),
  filename= "MCMCmulti_line3_repro.pdf")

#######

# line 4:

QQ<-2
kprime <- 5
K <- 5
N <- length(unlist(NELS88star_readmath$y)[NELS88star_readmath$studyID<=K])
Nprime <- min(c(N,suppressWarnings(min(N,min(c(1:N)[NELS88star_readmath$studyID>kprime])-1))), na.rm=TRUE)
jags.m <- jags.model(textConnection(multiMEMA), data = list (
  N = N,
  Nprime = Nprime,
  Q = QQ,
  NStudies = K,
  studyID = NELS88star_readmath$studyID[1:Nprime],
  Y = unlist(NELS88star_readmath$y)[1:Nprime],
  Xstar = do.call(rbind,lapply(NELS88star_readmath$X, function(x) {(cbind(1,x))}))[1:Nprime,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ),
  zeromat=as.array(matrix(diag(rep(1,QQ)),nrow=QQ)/(1e+20)),
  rho=0.1
  ),
          inits =  list(inits1, inits2, inits3),    
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps4 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps4)
line4_multi <- summary(samps4)$quantiles[,c(1,3,5)]




#######

# line 5:

QQ<-2
kprime<-0
K<-13
N<-length(unlist(NELS88star_readmath$y))
Nprime <- suppressWarnings(min(N,min(c(1:N)[NELS88star_readmath$studyID>kprime])-1))
jags.m <- jags.model(textConnection(multiMEMA), data = list (
  N = length(unlist(NELS88star_readmath$y)),
  Nprime = Nprime,
  Q = QQ,
  NStudies = K,
  studyID = NELS88star_readmath$studyID,
  Y = unlist(NELS88star_readmath$y),
  Xstar = do.call(rbind,lapply(NELS88star_readmath$X, function(x) {(cbind(1,x))}))[,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ),
  zeromat=as.array(matrix(diag(rep(1,QQ)),nrow=QQ)/(1e+20)),
  rho=0.1
  ),
          inits =  list(inits1, inits2, inits3),      
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps5 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps5)
line5_multi <- summary(samps5)$quantiles[,c(1,3,5)]




MCMCtrace(samps5, params = c("theta"), 
 priors = cbind(rnorm(100000, 0, 10), rnorm(100000, 0, 10), rnorm(100000, 0, 10)) , main_den = c(			   
  TeX("Density $\\theta_{1}$"),
  TeX("Density $\\theta_{2}$"),
    TeX("Density $\\theta_{3}$")),
  main_tr = c(					   					                         
    TeX("Trace $\\theta_{1}$"),
    TeX("Trace $\\theta_{2}$"),
    TeX("Trace $\\theta_{3}$")),
  filename= "MCMCmulti_line5_repro.pdf")


# params <- c("zeta1", "zeta2")		
# samps5zeta <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)

# MCMCtrace(samps5zeta, params = c("zeta1", "zeta2"), 
 # priors = cbind(rexp(100000, 0.1), rexp(100000, 0.1)) , main_den = c(			   
  # TeX("Density $\\zeta_{1}$"),
  # TeX("Density $\\zeta_{2}$")),
  # main_tr = c(					   					                         
    # TeX("Trace $\\zeta_{1}$"),
    # TeX("Trace $\\zeta_{2}$")),
  # filename= "MCMCmulti_line5_zeta.pdf")



#######

# line 6:


QQ<-2
kprime<-5
K<-13
N<-length(unlist(NELS88_readmath$y))

Nprime <- suppressWarnings(min(N,min(c(1:N)[NELS88_readmath$studyID>kprime])-1))
jags.m <- jags.model(textConnection(multiMEMA), data = list (
  N = length(unlist(NELS88_readmath$y)),
  Nprime = Nprime,
  Q = QQ,
  NStudies = K,
  studyID = NELS88_readmath$studyID,
  Y = unlist(NELS88_readmath$y),
  Xstar = do.call(rbind,lapply(NELS88_readmath$X, function(x) {(cbind(1,x))}))[,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ),
  zeromat=as.array(matrix(diag(rep(1,QQ)),nrow=QQ)/(1e+20)),
  rho = 0.1
  ),
          inits =  list(inits1, inits2, inits3),        
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps6 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps6)
line6_multi <- summary(samps6)$quantiles[,c(1,3,5)]



#######

# line 7:


QQ<-2
kprime<-0
K<-13
N<-length(unlist(NELS88_readmath$y))

Nprime <- suppressWarnings(min(N,min(c(1:N)[NELS88_readmath$studyID>kprime])-1))
jags.m <- jags.model(textConnection(multiMEMA), data = list (
  N = length(unlist(NELS88_readmath$y)),
  Nprime = Nprime,
  Q = QQ,
  NStudies = K,
  studyID = NELS88_readmath$studyID,
  Y = unlist(NELS88_readmath$y),
  Xstar = do.call(rbind,lapply(NELS88_readmath$X, function(x) {(cbind(1,x))}))[,1:(QQ+1)],
  D1 = matrix(diag(rep(1,QQ)),nrow=QQ)*1,
  D100 = matrix(diag(rep(1,QQ)),nrow=QQ)*100,
  m0 = rep(0,QQ),
  zeromat=as.array(matrix(diag(rep(1,QQ)),nrow=QQ)/(1e+20)),
  rho = 0.1
  ),
          inits =  list(inits1, inits2, inits3),        
n.chains = 3, n.adapt = 10000)


params <- c("theta", "tau")			
samps7 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps7)
line7_multi <- summary(samps7)$quantiles[,c(1,3,5)]

#####


makeTab<-function(linex_multi){c(
round(linex_multi["theta[1]",2],2),
paste(round(linex_multi["theta[1]",1],2),round(linex_multi["theta[1]",3],2),sep=", "),
round(linex_multi["theta[2]",2],2),
paste(round(linex_multi["theta[2]",1],2),round(linex_multi["theta[2]",3],2),sep=", "), round(linex_multi["theta[3]",2],2),
paste(round(linex_multi["theta[3]",1],2),round(linex_multi["theta[3]",3],2),sep=", "),
round(linex_multi["tau[2]","50%"],2),
round(linex_multi["tau[3]","50%"],2))
}

resultstab<-rbind(
makeTab(line1_multi),
makeTab(line2_multi),
makeTab(line3_multi),
makeTab(line4_multi),
makeTab(line5_multi),
makeTab(line6_multi),
makeTab(line7_multi))

xtable(resultstab)