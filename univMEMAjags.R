
########################################################################################
### MEMA: Measurement Error in Meta-Analysis
### contact: harlan.campbell@stat.ubc.ca



#####################
# Fresh start 

rm(list = ls(all = TRUE))
set.seed(87654321)

#####################
# Determine missing packages and load them:
required_packages <- c("MCMCvis", "xtable", "RCurl", "runjags", "rjags", 
                       "sjstats", "rlist", "tidyverse", "xtable", "MCMCpack",
                       "latex2exp")
not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed, type="source")                                           
invisible(suppressWarnings(lapply(required_packages, require, character.only = TRUE)))
ls()
rm("not_installed", "required_packages")
ls()


# Initial values for MCMC random seeds

	inits1<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(123))
	inits2<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(456))
	inits3<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(789))

#################################################################################
######################## load original schooldata  (NELS88) #################################
### dataset :  NELS88
schooldata <- read.csv("https://raw.githubusercontent.com/harlanhappydog/MEMA/master/13schools.csv")

n <- X <- y <- list()
se_alpha_hat <- se_beta_hat <- beta <- alpha <- alpha_hat <- beta_hat <- sigma_hat <- lambda_hat <- mu_hat <- w <- u <- vector()
K <- length(table(schooldata$sch))

for(i in 1:K){
  mydat 		    <- schooldata[schooldata$sch==i,]
  n[[i]] 		    <- dim(mydat)[1]
  y[[i]] 		    <- mydat$sci
  X[[i]] 		    <- mydat$rdg
  mod_i 		    <- lm(y[[i]] ~ X[[i]])
  
  alpha_hat[i] 	<- coef(mod_i)[1]
  beta_hat[i]     <- coef(mod_i)[2]
  se_alpha_hat[i] <- coefficients(summary(mod_i))[1,2]	
  se_beta_hat[i]  <- coefficients(summary(mod_i))[2,2]
  sigma_hat[i] 	<- summary(mod_i)$sigma
  lambda_hat[i] 	<- sqrt(var(X[[i]]))
  mu_hat[i] 	    <- mean(X[[i]])
}

#	c((sigma_hat[i])/((lambda_hat[i])* sqrt(n[[i]] -1)), se_beta_hat[i])


studyID <- schooldata$sch




NELS88_dataframe <- data.frame(
  "n_i"		 =	paste(as.numeric(table(studyID))),
  "alpha"		 = 	alpha_hat,
  "beta"		 = 	beta_hat, 
  "sigma2"	     = 	sigma_hat^2, 
  "mu"		     = 	mu_hat, 
  "lambda"	     = 	lambda_hat)	


print(xtable(data.frame(
  "n_i"		 =	paste(as.numeric(table(studyID))),
  "beta"		 = 	beta_hat, 
  "se_beta"	 = 	se_beta_hat, 
  "alpha"		 = 	alpha_hat,
  "se_alpha"	 = 	se_alpha_hat,
  "beta"		 = 	beta_hat, 
  "se_beta"	 = 	se_beta_hat, 	
  "sigma2"	     = 	sigma_hat^2), digits=2))	




NELS88IPD <-list (
  K             = length(table(schooldata$sch)),		
  n_per_study 	= as.numeric(table(schooldata$sch)),
  alpha_hat 	    = alpha_hat,
  beta_hat 	    = beta_hat,
  mu_hat        = mu_hat,
  sigma_hat 	    = sigma_hat,
  lambda_hat    = lambda_hat,
  y		        = y,
  X		        = X,
  n		        = n)

invtXX<-list()
for(k in 1:K){
  solve(t(cbind(1,X[[k]]))%*%cbind(1,X[[k]]))
  x11<-(lambda_hat[k]^2+mu_hat[k]^2)/(lambda_hat[k]^2 * (NELS88IPD$n_per_study[k]-1))
  x22<-1/(lambda_hat[k]^2* (NELS88IPD$n_per_study[k]-1))
  x12<-(-mu_hat[k])/(lambda_hat[k]^2* (NELS88IPD$n_per_study[k]-1))
  invtXX[[k]] <- matrix(c(x11,x12,x12,x22),2,2)}

NELS88 <-list (
  K 	           = length(table(schooldata$sch)),
  kprime       = length(table(schooldata$sch)),
  n_per_study  = as.numeric(table(schooldata$sch)),
  alpha_hat 	   = alpha_hat,
  beta_hat 	   = beta_hat,
  mu_hat 	   = mu_hat,
  sigma_hat 	   = sigma_hat,
  lambda_hat   = lambda_hat,
  se_alpha_hat = se_alpha_hat, 
  se_beta_hat  = se_beta_hat)

rm(alpha, alpha_hat, beta, beta_hat, i, K, lambda_hat, mod_i, 
   mu_hat, mydat, n, schooldata, sigma_hat, studyID, u, w, X, y, 
   se_alpha_hat, se_beta_hat)

#################################################################################
##### load schooldata with added measurement error in K-kprime schools (NELS88_star) ###
### dataset :  NELS88_star

kprime <- 5

schooldata <- read.csv("https://raw.githubusercontent.com/harlanhappydog/MEMA/master/13schools.csv")
schooldata_me 	<- schooldata
K <- length(table(schooldata_me$sch))

PHI_ <- rep(0,K)
PHI_[1:kprime] <- 0

if(kprime<K){  PHI_[(kprime+1):K] <- seq(1, 12, length.out=length((kprime+1):K))  }

PHI_squared 	<- PHI_^2
gamma <- ((1+(PHI_^2/NELS88$lambda_hat^2))^(-1))


for(i in 1:K){
  schooldata_me[schooldata_me$sch==i,]$rdg <- schooldata[schooldata$sch==i,]$rdg +
    rnorm(length(schooldata[schooldata$sch==i,]$rdg),0, PHI_[i])
}

n_me <- X_me <- y_me <- list()
se_alpha_hat_me <- se_beta_hat_me <- beta_me <- alpha_me <- alpha_hat_me <- beta_hat_me <- vector()
sigma_hat_me <- lambda_hat_me <- mu_hat_me <- w <- u <- vector()

for(i in 1:K){
  mydat_me 		<- schooldata_me[schooldata_me$sch==i,]
  n_me[[i]] 		<- dim(mydat_me)[1]
  y_me[[i]] 		<- mydat_me$sci
  X_me[[i]] 		<- mydat_me$rdg
  mod_i_me 		<- lm(y_me[[i]] ~ X_me[[i]])
  
  alpha_hat_me[i]  <- coef(mod_i_me)[1]
  beta_hat_me[i] 	 <- coef(mod_i_me)[2]
  sigma_hat_me[i]  <- summary(mod_i_me)$sigma
  lambda_hat_me[i] <- sqrt(var(X_me[[i]]))
  mu_hat_me[i] 	 <- mean(X_me[[i]])
  
  se_alpha_hat_me[i] <- coefficients(summary(mod_i_me))[1,2]	
  se_beta_hat_me[i]  <- coefficients(summary(mod_i_me))[2,2]
}

studyID <- schooldata_me$sch

NELS88star_dataframe <- data.frame(
  "n_i"		= paste(as.numeric(table(studyID))), 
  "PHI_"		= PHI_, 
  "gamma" 	    = gamma, 
  "alpha"		= alpha_hat_me,
  "beta"		= beta_hat_me, 
  "sigma"		= sigma_hat_me, 
  "mu"			= mu_hat_me, 
  "lambda"	    = lambda_hat_me)	


NELS88star <- list (
  K 	= length(table(schooldata_me$sch)),
  kprime       = kprime,		
  n_per_study 	= as.numeric(table(schooldata_me$sch)),
  alpha_hat 	= alpha_hat_me,
  beta_hat 	= beta_hat_me,
  mu_hat 		= mu_hat_me,
  sigma_hat 	= sigma_hat_me,
  lambda_hat 	= lambda_hat_me,
  PHI_		= PHI_,
  gamma		= gamma,
  se_alpha_hat = se_alpha_hat_me,
  se_beta_hat = se_beta_hat_me)	


print(xtable(data.frame(
  "n_i"		 =	paste(as.numeric(table(studyID))),
  "PHI_"	  	 =  PHI_,
  "gamma" 	 =  gamma,
  "beta"		 = 	beta_hat_me, 
  "se_beta"	 = 	se_beta_hat_me, 
  "alpha"		 = 	alpha_hat_me,
  "se_alpha"	 = 	se_alpha_hat_me,
  "beta"		 = 	beta_hat_me, 
  "se_beta"	 = 	se_beta_hat_me, 	
  "sigma2"	 = 	sigma_hat_me^2), digits=2))	



xtable(data.frame(
  n_per_study 	= as.numeric(table(schooldata_me$sch)),
  alpha_hat 	= alpha_hat_me,
  beta_hat 	= beta_hat_me,
  sigma_hat 	= sigma_hat_me,
  mu_hat 		= mu_hat_me,  
  lambda_hat 	= lambda_hat_me,
  PHI_		= PHI_,
  gamma		= gamma)	)

NELS88starIPD <- list (
  K 	= length(table(schooldata_me$sch)),		
  n_per_study 	= as.numeric(table(schooldata_me$sch)),
  alpha_hat 	= alpha_hat_me,
  beta_hat 	= beta_hat_me,
  mu_hat 		= mu_hat_me,
  sigma_hat 	= sigma_hat_me,
  lambda_hat 	= lambda_hat_me,
  PHI_		= PHI_,
  gamma		= gamma,
  y		     = y_me,
  X            = X_me,
  n            = n_me)

rm(se_alpha_hat_me , se_beta_hat_me , alpha_hat_me, alpha_me, beta_hat_me, beta_me, gamma, 
   i, k, K, lambda_hat_me, mod_i_me, mu_hat_me, mydat_me, 
   n_me, schooldata, schooldata_me, sigma_hat_me, studyID, 
   PHI_, PHI_squared, u, w, X_me, y_me)

round(mean(NELS88star$gamma),2)
round(var(NELS88star$gamma),2)

#################################################################################
######################## model code in JAGS for two models: BayesMA and MEMA ###########
############


BayesMA_onlybeta <- 'model{

  theta ~ dnorm(0, 1/100);  
  TAU ~ dt(0, 0.25, 1)I(0,);
  
 for (k in 1:K){
	beta_hat[k] ~ dnorm(theta, 1/((TAU^2) + (se_beta_hat[k]^2)));
	}

}'

#   For small n, we should perhaps use:
#	beta_hat[k] ~ dnorm(theta, 1/((TAU^2)+((n[k]-1)/n[k])*(se_beta_hat[k]^2)));

MEMA_onlybeta <- 'model{
 
  
 for (k in 1:kprime){
    gamma_hat[k] = 1;
  	} 

 for (k in (kprime+1):K){
	gamma_hat[k] ~ dunif(0,1);
    }

  theta ~ dnorm(0, 1/100);  
  TAU ~ dt(0, 0.25, 1)I(0,);

  for (k in 1:K){
    beta_star[k] ~ dnorm(gamma_hat[k]*beta[k], 1/(se_beta_hat_star[k]^2));
    beta[k] ~ dnorm(theta, 1/(TAU^2));  
    }
    
}'






BayesMA <- 'model{

  theta ~ dnorm(0, 1/100);  
  XI_ ~ dnorm(0, 1/100);  
  
  TAU ~ dt(0, 0.25, 1)I(0,);
  OMEGA_ ~ dt(0, 0.25, 1)I(0,);
  rho_ ~ dunif(-1, 1)

  TauMat[1,1] <- OMEGA_^2  
  TauMat[2,2] <- TAU^2
  TauMat[1, 2] <- rho_ * TAU * OMEGA_
  TauMat[2, 1] <- TauMat[1, 2]  
  

  for (k in 1:K){	
  coef_true[k,1:2] ~ dmnorm(c(XI_, theta), inverse(TauMat))
  }
  
 for (k in 1:K){

 mu_hat[k] = sqrt((se_alpha_hat[k]/se_beta_hat[k])^2 - (sigma_hat[k]^2)/(n[k]*se_beta_hat[k]^2))

lambda_hat[k] = sigma_hat[k]/(sqrt(n[k]-1)*se_beta_hat[k])

    cov[k,1,1] = ((lambda_hat[k]^2 + mu_hat[k]^2) * sigma_hat[k]^2) / ((lambda_hat[k]^2) * n_per_study[k]);
    cov[k,1,2] =  - mu_hat[k] * (sigma_hat[k]^2)/((lambda_hat[k]^2) * n_per_study[k]);    
    cov[k,2,1] =  - mu_hat[k] * (sigma_hat[k]^2)/((lambda_hat[k]^2) * n_per_study[k]);
    cov[k,2,2] =  (sigma_hat[k]^2)/((lambda_hat[k]^2) * n_per_study[k]);

	coef_hat[k, 1:2] ~ dmnorm(coef_true[k,1:2], inverse(cov[k, 1:2, 1:2]));
	}

}'


MEMA <- 'model{
 
 for (k in 1:kprime){
    PHI_squared[k] = 0;
    lambda_hat[k]  = sqrt(lambda_hat_star[k]^2 - PHI_squared[k]);
    gamma_hat[k] = 1;
  	} 

 for (k in (kprime+1):K){
    iPHI_squared[k] ~ dgamma(zeta1,zeta2)T(1/(lambda_hat_star[k]^2),);
    PHI_squared[k] <- 1/iPHI_squared[k]    

    lambda_hat[k] = sqrt(lambda_hat_star[k]^2 - PHI_squared[k]);
   	gamma_hat[k] = (lambda_hat[k]^2)/(lambda_hat[k]^2 + PHI_squared[k]);
    }
    
  for (k in 1:K){	
	mu_hat[k]      = mu_hat_star[k];
	sigma_hat[k]   = sqrt((sigma_hat_star[k])^2 - 
	         (1 - gamma_hat[k])*(coef_true[k,2]^2) * (lambda_hat[k]^2));
	}

  zeta1 ~ dexp(delta);
  zeta2 ~ dexp(delta);
  
  theta ~ dnorm(0, 1/100);  
  XI_ ~ dnorm(0, 1/100);  
  TAU ~ dt(0, 0.25, 1)I(0,);
  OMEGA_ ~ dt(0, 0.25, 1)I(0,);
  rho_ ~ dunif(-1, 1)
  TauMat[1,1] <- OMEGA_^2  
  TauMat[2,2] <- TAU^2
  TauMat[1, 2] <- rho_ * TAU * OMEGA_
  TauMat[2, 1] <- TauMat[1, 2]  
  
  for (k in 1:K){	

   coef_true[k,1:2] ~ dmnorm(c(XI_, theta), inverse(TauMat))


    cov[k,1,1] = ((lambda_hat_star[k]^2 + mu_hat_star[k]^2) * 
           sigma_hat_star[k]^2) / ((lambda_hat_star[k]^2) * n_per_study[k]);
    
    cov[k,1,2] =  - mu_hat_star[k] * 
           (sigma_hat_star[k]^2)/((lambda_hat_star[k]^2) * n_per_study[k]);
    
    cov[k,2,1] =  - mu_hat_star[k] * 
           (sigma_hat_star[k]^2)/((lambda_hat_star[k]^2) * n_per_study[k]);

    cov[k,2,2] =  (sigma_hat_star[k]^2)/((lambda_hat_star[k]^2) * n_per_study[k]);

  	alpha_star[k] = coef_true[k,1] + (1 - gamma_hat[k]) * coef_true[k,2] * mu_hat[k];
	
    beta_star[k] = gamma_hat[k] * coef_true[k,2];
    	
	coef_hat[k, 1:2] ~ dmnorm(c(alpha_star[k], beta_star[k]), inverse(cov[k, 1:2, 1:2]));
    }
    
}'


######################################################
#### run models with data:
######################################################




jags.m <- jags.model(textConnection(BayesMA_onlybeta), 
                     data = list (
                       K  =  length(NELS88$beta_hat),  				 			                   beta_hat =  NELS88$beta_hat,
                       se_beta_hat  =  NELS88$se_beta_hat),
  inits =  list(inits1, inits2, inits3),  
                     n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")			
samps <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
summary_samps <- summary(samps)$quantiles[,c(1, 3, 5)]
line0a <- round(summary_samps, 2)
line0a


jags.m <- jags.model(textConnection(BayesMA), data = list (
  K 	            =  length(NELS88$beta_hat),
  n_per_study 	=  NELS88$n_per_study,
  coef_hat 	    =  cbind(NELS88$alpha_hat, NELS88$beta_hat),
  se_alpha_hat    =  NELS88$se_alpha_hat,
  se_beta_hat     =  NELS88$se_beta_hat,  					
  sigma_hat       =  NELS88$sigma_hat,
  n               =  NELS88$n),
  inits =  list(inits1, inits2, inits3),  
  n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")			
samps <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
summary_samps <- summary(samps)$quantiles[,c(1, 3, 5)]
line0 <- round(summary_samps, 2)
line0




######  line 1.

######################################################



jags.m <- jags.model(textConnection(MEMA_onlybeta), data = list (
  K 	            =  length(NELS88star$beta_hat),
  beta_star 	     =  NELS88star$beta_hat,
  se_beta_hat_star =  NELS88star$se_beta_hat,
  kprime          =  13),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")			
samps1 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
#plot(samps1)
summary_samps <- summary(samps1)$quantiles[,c(1, 3, 5)]
line1a <- round(summary_samps, 2)
line1a


jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat),
  n_per_study 	=  NELS88star$n_per_study,
  coef_hat 	    =  cbind(NELS88star$alpha_hat, NELS88star$beta_hat),
  mu_hat_star 	=  NELS88star$mu_hat,
  sigma_hat_star  =  NELS88star$sigma_hat,
  lambda_hat_star =  NELS88star$lambda_hat,
  delta           = 0.1,
  kprime          =  13),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")			
samps1 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps1)
summary_samps <- summary(samps1)$quantiles[,c(1, 3, 5)]
line1 <- round(summary_samps, 2)
line1





jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat),
  n_per_study 	=  NELS88star$n_per_study,
  coef_hat 	    =  cbind(NELS88star$alpha_hat, NELS88star$beta_hat),
  mu_hat_star 	=  NELS88star$mu_hat,
  sigma_hat_star  =  NELS88star$sigma_hat,
  lambda_hat_star =  NELS88star$lambda_hat,
  delta           = 0.5,
  kprime          =  13),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")			
samps1 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps1)
summary_samps <- summary(samps1)$quantiles[,c(1, 3, 5)]
line1b <- round(summary_samps, 2)
line1b
######  line 2.

######################################################


jags.m <- jags.model(textConnection(MEMA_onlybeta), data = list (
  K 	            =  length(NELS88$beta_hat),
  beta_star 	     =  NELS88$beta_hat,
  se_beta_hat_star =  NELS88$se_beta_hat,
  kprime          =  13),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")			
samps2a <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps2a)
summary_samps <- summary(samps2a)$quantiles[,c(1, 3, 5)]
line2a <- round(summary_samps, 2)
line2a

jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88$beta_hat),
  n_per_study 	=  NELS88$n_per_study,
  coef_hat 	    =  cbind(NELS88$alpha_hat, NELS88$beta_hat),
  mu_hat_star 	=  NELS88$mu_hat,
  sigma_hat_star  =  NELS88$sigma_hat,
  lambda_hat_star =  NELS88$lambda_hat,
  delta             =  0.1,
  kprime          =  13),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps2 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps2)
summary_samps <- summary(samps2)$quantiles[,c(1, 3, 5)]
line2 <- round(summary_samps, 2)


jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88$beta_hat),
  n_per_study 	=  NELS88$n_per_study,
  coef_hat 	    =  cbind(NELS88$alpha_hat, NELS88$beta_hat),
  mu_hat_star 	=  NELS88$mu_hat,
  sigma_hat_star  =  NELS88$sigma_hat,
  lambda_hat_star =  NELS88$lambda_hat,
  delta             =  0.5,
  kprime          =  13),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps2 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps2)
summary_samps <- summary(samps2)$quantiles[,c(1, 3, 5)]
line2b <- round(summary_samps, 2)


######  line 3.

######################################################




jags.m <- jags.model(textConnection(MEMA_onlybeta), 
                     data = list (
                       K 	            =  length(NELS88star$beta_hat),
                       beta_star 	     =  NELS88star$beta_hat,
                       se_beta_hat_star =  NELS88star$se_beta_hat,
                       kprime          =  5),
 inits =  list(inits1, inits2, inits3),      
                     n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")		
samps3 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps3)
summary_samps <- summary(samps3)$quantiles[,c(1, 3, 5)]
line3a <- round(summary_samps, 2)
line3a


                  
                  
	
jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat),
  n_per_study 	=  NELS88star$n_per_study,
  coef_hat 	    =  cbind(NELS88star$alpha_hat, NELS88star$beta_hat),
  mu_hat_star 	=  NELS88star$mu_hat,
  sigma_hat_star  =  NELS88star$sigma_hat,
  lambda_hat_star =  NELS88star$lambda_hat,
  delta             =  0.1,
  kprime          =  5),
  inits =  list(inits1, inits2, inits3),  
  n.chains = 3, n.adapt = 10000)

params <- c("theta", "XI_", "TAU", "OMEGA_")		
samps3 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
head(samps3)
plot(samps3)
summary_samps <- summary(samps3)$quantiles[,c(1, 3, 5)]
line3 <- round(summary_samps, 2)
line3



# theta ~ dnorm(0, 1/100);  
# XI_ ~ dnorm(0, 1/100);  
# TAU ~ dt(0, 0.25, 1)I(0,);
# OMEGA_ ~ dt(0, 0.25, 1)I(0,);


# MCMCtrace(samps3, params = c("theta", "XI_", "TAU", "OMEGA_"), 
# priors = cbind(rnorm(100000, 0, 10), rnorm(100000, 0, 10), rcauchy(100000, 0, sqrt(1/0.25)), rcauchy(100000, 0, sqrt(1/0.25))) , main_den = c(			   
# TeX("Density $\\theta$"),
# TeX("Density $\\xi$"),
# TeX("Density $\\tau$"),
# TeX("Density $\\omega$")),
# main_tr = c(					   					                         
# TeX("Trace $\\theta$"),
# TeX("Trace $\\xi$"),
# TeX("Trace $\\tau$"),
# TeX("Trace $\\omega$")),
# filename= "MCMC_line3.pdf")

# uncomment:
MCMCtrace(samps3, params = c("theta", "XI_"), 
          priors = cbind(rnorm(100000, 0, 10), rnorm(100000, 0, 10)) , main_den = c(			   
            TeX("Density $\\theta$"),
            TeX("Density $\\xi$")),
          main_tr = c(					   					                         
            TeX("Trace $\\theta$"),
            TeX("Trace $\\xi$")),  
          filename= "MCMC_line3_repro_.pdf")


params <- c("zeta1", "zeta2")		
samps3zeta <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)

MCMCtrace(samps3zeta, params = c("zeta1", "zeta2"), 
          priors = cbind(rexp(100000, 0.1), rexp(100000, 0.1)) , main_den = c(			   
            TeX("Density $\\zeta_{1}$"),
            TeX("Density $\\zeta_{2}$")),
          main_tr = c(					   					                         
            TeX("Trace $\\zeta_{1}$"),
            TeX("Trace $\\zeta_{2}$")),
          filename= "MCMC_line3_zeta_repro_.pdf")



jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat),
  n_per_study 	=  NELS88star$n_per_study,
  coef_hat 	    =  cbind(NELS88star$alpha_hat, NELS88star$beta_hat),
  mu_hat_star 	=  NELS88star$mu_hat,
  sigma_hat_star  =  NELS88star$sigma_hat,
  lambda_hat_star =  NELS88star$lambda_hat,
  delta             =  0.5,
  kprime          =  5),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "XI_", "TAU", "OMEGA_")		
samps3 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps3)
summary_samps <- summary(samps3)$quantiles[,c(1, 3, 5)]
line3b <- round(summary_samps, 2)
line3b

######  line 4.

######################################################


jags.m <- jags.model(textConnection(MEMA_onlybeta), 
                     data = list (
                       K 	            =  length(NELS88star$beta_hat[1:5]),
                       beta_star 	     =  NELS88star$beta_hat[1:5],
                       se_beta_hat_star =  NELS88star$se_beta_hat[1:5],
                       kprime          =  5),
 inits =  list(inits1, inits2, inits3),
                     n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")		
samps4 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps4)
summary_samps <- summary(samps4)$quantiles[,c(1, 3, 5)]
line4a <- round(summary_samps, 2)
line4a



jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat[1:5]),
  n_per_study 	=  NELS88star$n_per_study[1:5],
  coef_hat 	    =  cbind(NELS88star$alpha_hat[1:5], NELS88star$beta_hat[1:5]),
  mu_hat_star 	=  NELS88star$mu_hat[1:5],
  sigma_hat_star  =  NELS88star$sigma_hat[1:5],
  lambda_hat_star =  NELS88star$lambda_hat[1:5],
  delta             =  0.1,
  kprime          =  5),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps4 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps4)
summary_samps <- summary(samps4)$quantiles[,c(1, 3, 5)]
line4 <- round(summary_samps, 2)
line4




jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat[1:5]),
  n_per_study 	=  NELS88star$n_per_study[1:5],
  coef_hat 	    =  cbind(NELS88star$alpha_hat[1:5], NELS88star$beta_hat[1:5]),
  mu_hat_star 	=  NELS88star$mu_hat[1:5],
  sigma_hat_star  =  NELS88star$sigma_hat[1:5],
  lambda_hat_star =  NELS88star$lambda_hat[1:5],
  delta             =  0.5,
  kprime          =  5),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps4 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps4)
summary_samps <- summary(samps4)$quantiles[,c(1, 3, 5)]
line4b <- round(summary_samps, 2)
line4b
######  line 5.

######################################################



jags.m <- jags.model(textConnection(MEMA_onlybeta), 
                     data = list (
                       K 	            =  length(NELS88star$beta_hat),
                       beta_star 	     =  NELS88star$beta_hat,
                       se_beta_hat_star =  NELS88star$se_beta_hat,
                       kprime          =  0),
 inits =  list(inits1, inits2, inits3),                   
                     n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")		
samps5 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps5)
summary_samps <- summary(samps5)$quantiles[,c(1, 3, 5)]
line5a <- round(summary_samps, 2)
line5a


jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat),
  n_per_study 	=  NELS88star$n_per_study,
  coef_hat 	    =  cbind(NELS88star$alpha_hat, NELS88star$beta_hat),
  mu_hat_star 	=  NELS88star$mu_hat,
  sigma_hat_star  =  NELS88star$sigma_hat,
  lambda_hat_star =  NELS88star$lambda_hat,
  delta             =  0.1,
  kprime          =  0),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps5 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps5)
summary_samps <- summary(samps5)$quantiles[,c(1, 3, 5)]
line5 <- round(summary_samps, 2)
line5


jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88star$beta_hat),
  n_per_study 	=  NELS88star$n_per_study,
  coef_hat 	    =  cbind(NELS88star$alpha_hat, NELS88star$beta_hat),
  mu_hat_star 	=  NELS88star$mu_hat,
  sigma_hat_star  =  NELS88star$sigma_hat,
  lambda_hat_star =  NELS88star$lambda_hat,
  delta             =  0.5,
  kprime          =  0),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps5 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps5)
summary_samps <- summary(samps5)$quantiles[,c(1, 3, 5)]
line5b <- round(summary_samps, 2)
line5b


######  line 6.

######################################################


jags.m <- jags.model(textConnection(MEMA_onlybeta), 
                     data = list (
                       K 	            =  length(NELS88$beta_hat),
                       beta_star 	     =  NELS88$beta_hat,
                       se_beta_hat_star =  NELS88$se_beta_hat,
                       kprime          =  5),
 inits =  list(inits1, inits2, inits3),                      
                     n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")		
samps6 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps6)
summary_samps <- summary(samps6)$quantiles[,c(1, 3, 5)]
line6a <- round(summary_samps, 2)
line6a


jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88$beta_hat),
  n_per_study 	=  NELS88$n_per_study,
  coef_hat 	    =  cbind(NELS88$alpha_hat, NELS88$beta_hat),
  mu_hat_star 	=  NELS88$mu_hat,
  sigma_hat_star  =  NELS88$sigma_hat,
  lambda_hat_star =  NELS88$lambda_hat,
  delta             =  0.1,
  kprime          =  5),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps6 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps6)
summary_samps <- summary(samps6)$quantiles[,c(1, 3, 5)]
line6 <- round(summary_samps, 2)


jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88$beta_hat),
  n_per_study 	=  NELS88$n_per_study,
  coef_hat 	    =  cbind(NELS88$alpha_hat, NELS88$beta_hat),
  mu_hat_star 	=  NELS88$mu_hat,
  sigma_hat_star  =  NELS88$sigma_hat,
  lambda_hat_star =  NELS88$lambda_hat,
  delta             =  0.5,
  kprime          =  5),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps6 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps6)
summary_samps <- summary(samps6)$quantiles[,c(1, 3, 5)]
line6b <- round(summary_samps, 2)

######  line 7.

######################################################


jags.m <- jags.model(textConnection(MEMA_onlybeta), 
                     data = list (
                       K 	            =  length(NELS88$beta_hat),
                       beta_star 	     =  NELS88$beta_hat,
                       se_beta_hat_star =  NELS88$se_beta_hat,
                       kprime          =  0),
 inits =  list(inits1, inits2, inits3),                   
                     n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU")		
samps7 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps7)
summary_samps <- summary(samps7)$quantiles[,c(1, 3, 5)]
line7a <- round(summary_samps, 2)
line7a



jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88$beta_hat),
  n_per_study 	=  NELS88$n_per_study,
  coef_hat 	    =  cbind(NELS88$alpha_hat, NELS88$beta_hat),
  mu_hat_star 	=  NELS88$mu_hat,
  sigma_hat_star  =  NELS88$sigma_hat,
  lambda_hat_star =  NELS88$lambda_hat,
  delta             =  0.1,
  kprime          =  0),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps7 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps7)
summary_samps <- summary(samps7)$quantiles[,c(1, 3, 5)]
line7 <- round(summary_samps, 2)



jags.m <- jags.model(textConnection(MEMA), data = list (
  K 	            =  length(NELS88$beta_hat),
  n_per_study 	=  NELS88$n_per_study,
  coef_hat 	    =  cbind(NELS88$alpha_hat, NELS88$beta_hat),
  mu_hat_star 	=  NELS88$mu_hat,
  sigma_hat_star  =  NELS88$sigma_hat,
  lambda_hat_star =  NELS88$lambda_hat,
  delta             =  0.5,
  kprime          =  0),
    inits =  list(inits1, inits2, inits3),    n.chains = 3, n.adapt = 10000)

params <- c("theta", "TAU", "XI_", "OMEGA_", "rho_")	
samps7 <- coda.samples(jags.m, params, n.iter = 100000, thin = 10)
plot(samps7)
summary_samps <- summary(samps7)$quantiles[,c(1, 3, 5)]
line7b <- round(summary_samps, 2)



tablesummary<-function(linex){c(theta=linex["theta",c("50%")], theta95=paste(linex["theta",c("2.5%")],linex["theta", c("97.5%")], sep=" , "), tau=linex["TAU","50%"])}

xtable(cbind(rep(c(0,1,2),7), sort(rep(1:7,3)), rbind(
  tablesummary(line1a),
  tablesummary(line1),
  tablesummary(line1b),
  tablesummary(line2a),
  tablesummary(line2),
  tablesummary(line2b),
  tablesummary(line3a),
  tablesummary(line3),
  tablesummary(line3b),
  tablesummary(line4a),
  tablesummary(line4),
  tablesummary(line4b),
  tablesummary(line5a),
  tablesummary(line5),
  tablesummary(line5b),
  tablesummary(line6a),
  tablesummary(line6),
  tablesummary(line6b),
  tablesummary(line7a),
  tablesummary(line7),
  tablesummary(line7b)
)))

### Table 1

schooldata <- read.csv("https://raw.githubusercontent.com/harlanhappydog/MEMA/master/13schools.csv")
schooldata_me 	<- schooldata

xtable(data.frame(
   n_per_study 	= as.numeric(table(schooldata_me$sch)),
   phi          = NELS88star$PHI_,
   gamma        = NELS88star$gamma,   
   beta_hat 	= NELS88star$beta_hat,
   se_beta_hat 	= NELS88star$se_beta_hat,
   alpha_hat 	= NELS88star$alpha_hat,
   se_alpha_hat	= NELS88star$se_alpha_hat,      
   beta_hat 	= NELS88star$beta_hat,
   se_beta_hat 	= NELS88star$se_beta_hat,
   sigma_hat 	= NELS88star$sigma_hat^2)	)




### Table 2
#linea = univariate
#line  = 0.1
#lineb = 0.5

xtable( rbind(
  tablesummary(line1a),
  tablesummary(line1),
  tablesummary(line2a),
  tablesummary(line2),
  tablesummary(line3a),
  tablesummary(line3),
  tablesummary(line3b),
  tablesummary(line4a),
  tablesummary(line4),
  tablesummary(line5a),
  tablesummary(line5),
  tablesummary(line5b),
  tablesummary(line6a),
  tablesummary(line6),
  tablesummary(line6b),
  tablesummary(line7a),
  tablesummary(line7),
  tablesummary(line7b)
))



