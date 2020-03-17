
##########################################################
######################## DATASETS ########################
##########################################################

schooldata <- read.csv("~/Desktop/UBC/RECODID_ZIKV/Rcode/13schools.csv")

########## original schooldata with X = [reading]  (ADNELS88) 
### dataset :  NELS88_read

XtXinv <- lambda_hat <- mu_hat <- sigma_hat <- beta_hat  <- n <- X  <- COV <- y <- list()
nvec <- vector()
K <- length(table(schooldata$sch))
for(i in 1:K){
	nvec[i] <- dim(schooldata[schooldata$sch==i,])[1]	
	mydat 	<- schooldata[schooldata$sch==i,]
	n[[i]] 	<- dim(mydat)[1]
	y[[i]] 	<- mydat$sci
	X[[i]] 	<- cbind(mydat$rdg)
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
NELS88_read<-list()
NELS88_read$N 		<- sum(unlist(n))
NELS88_read$studyID <- schooldata$sch
NELS88_read$Q 		<- dim(X[[1]])[2]
NELS88_read$K 		<- length(table(schooldata$sch))
NELS88_read$y 		<- y
NELS88_read$X 		<- X
NELS88_read$nvec 	<- nvec
NELS88_read$beta_hat <- beta_hat
NELS88_read$sigma_hat <- sigma_hat
NELS88_read$COV <- COV
NELS88_read$XtXinv <- XtXinv
NELS88_read$lambda_hat <- lambda_hat
NELS88_read$mu_hat <- mu_hat

# cleaning up:
rm(mydat, mod_i, XtXinv, lambda_hat, mu_hat, sigma_hat, beta_hat, n, X, COV, y, nvec, K, i)
ls()
#### #### #### #### #### #### #### #### #### #### #### #### 

schooldata <- read.csv("~/Desktop/UBC/RECODID_ZIKV/Rcode/13schools.csv")

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

#### #### #### #### #### #### #### #### #### #### #### #### 



########## original schooldata with X = [reading, math]  (NELS88) 
### dataset :  NELS88star_readmath

XtXinv_me <- lambda_hat_me <- mu_hat_me <- sigma_hat_me <- beta_hat_me  <- n_me <- X_me  <- COV_me <- y_me <- list()
nvec_me <- vector()


schooldata_me <- schooldata

K <- length(table(schooldata_me$sch))
Lambda22<-NULL
for(k in 1:K){Lambda22<-c(Lambda22, NELS88_readmath$lambda_hat[[k]][2,2]) }


set.seed(123)
sqrtPsi <- rep(0, NELS88_readmath$K)
for(k in 1:round(length(sqrtPsi)/3)){ sqrtPsi[k] <- 0.05*mean(sqrt(Lambda22))  }
for(k in (1+round(length(sqrtPsi)/3)):(2*round(length(sqrtPsi)/3))){ sqrtPsi[k] <- 0.5*mean(sqrt(Lambda22))  }
for(k in (1+2*round(length(sqrtPsi)/3)):length(sqrtPsi)){ sqrtPsi[k] <- 0.95*mean(sqrt(Lambda22)) }



K <- length(table(schooldata_me$sch))
	for(i in 1:K){
		schooldata_me[schooldata_me$sch==i,]$math <- schooldata[schooldata$sch==i,]$math +
				rnorm(length(schooldata[schooldata$sch==i,]$math),0, sqrtPsi[[i]])
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

# all summary data in a list:
NELS88star_readmath 		<- list()
NELS88star_readmath$N 		<- sum(unlist(n_me))
NELS88star_readmath$studyID <- schooldata$sch
NELS88star_readmath$Q 		<- dim(X_me[[1]])[2]
NELS88star_readmath$K 		<- length(table(schooldata$sch))
NELS88star_readmath$y 		<- y_me
NELS88star_readmath$X 		<- X_me
NELS88star_readmath$nvec 	<- nvec_me
NELS88star_readmath$beta_hat 	<- beta_hat_me
NELS88star_readmath$sigma_hat 	<- sigma_hat_me
NELS88star_readmath$COV 		<- COV_me
NELS88star_readmath$XtXinv 		<- XtXinv_me
NELS88star_readmath$lambda_hat  <- lambda_hat_me
NELS88star_readmath$mu_hat 		<- mu_hat_me
NELS88star_readmath$sqrtPsi <- sqrtPsi

# cleaning up:
rm(mydat_me, mod_i_me, XtXinv_me, lambda_hat_me, mu_hat_me, sigma_hat_me, beta_hat_me, n_me, X_me, COV_me, y_me, nvec_me, K, i, k, sqrtPsi, Lambda22)
ls()

#### #### #### #### #### #### #### #### #### #### #### #### 




########################################################################################
######################## MODELS ########################
########################################################################################

  

###########################################################################################

BayesMA_multi_AG <-
"data {
  int<lower=1> N;
  int<lower=0> Q;  			
  int<lower=1> NStudies;		
  int<lower=1,upper=NStudies> studyID[N];
  matrix[NStudies,Q+1] beta_hat;
  matrix[Q+1, Q+1] XtXinv[NStudies];
  vector<lower=0>[NStudies] sigma;
  vector[Q+1] thetaPrior_mu;
  matrix[Q+1, Q+1] thetaPrior_sigma; 
}


transformed data {
  matrix[Q+1, Q+1] COV[NStudies];
  
	for(k in 1:NStudies){
		COV[k] = XtXinv[k] *sigma[k]^2;
	}		 
 
 }


parameters {
  matrix[NStudies, Q+1] beta;
  vector[Q+1] theta;
  corr_matrix[Q+1] Omega_om;        // prior correlation
  vector<lower=0>[Q+1] tau_om;      // prior scale
}


transformed parameters {
  matrix[Q+1, Q+1] omega;	

  omega = quad_form_diag(Omega_om, tau_om); 
	
}

model {


    theta    ~   multi_normal(thetaPrior_mu, thetaPrior_sigma);  
    
    tau_om ~ cauchy(0, 2.5);
    
    Omega_om ~ lkj_corr(2);


  for (k in 1:NStudies)	  
	beta[k,] ~   multi_normal(theta, omega);

  for (k in 1:NStudies)	  
	beta_hat[k,] ~   multi_normal(beta[k,], COV[k]);
	
  	}"
  	
###################################################################################
##############################################################################
##############################################################################
#############


MA_multi_IPD <-
"data {
  int<lower=1> N;
  int<lower=1> Q;  			
  int<lower=1> NStudies;		
  int<lower=1,upper=NStudies> studyID[N];
  vector[N] Y;
  matrix[N, Q] X;
  int<lower=1> n_vec[NStudies];
}

parameters {
  vector<lower=0>[NStudies] sigma;
  matrix[NStudies, Q] mu;
  matrix[NStudies, Q+1] beta;
  vector[Q+1] theta;
  corr_matrix[Q+1] Omega_om;        // prior correlation for omega
  vector<lower=0>[Q+1] tau_om;      // prior scale for omega
  corr_matrix[Q] Omega_lambda[NStudies];        // prior correlation for lambda
  vector<lower=0>[Q] tau_lambda[NStudies];      // prior scale for lambda
}


transformed parameters {
  cov_matrix[Q+1] omega;
  cov_matrix[Q] lambda[NStudies];	
  
  omega = quad_form_diag(Omega_om, tau_om); 

  for(k in 1:NStudies){
	lambda[k] = quad_form_diag(Omega_lambda[k], tau_lambda[k]);
	}	
}

model {

  for(q in 1:(Q+1)){    theta[q]    ~   normal(0, 10);     }
    
    tau_om ~ cauchy(0, 2.5);
    Omega_om ~ lkj_corr(2);
	
  for (k in 1:NStudies){	  
	beta[k,] ~ multi_normal(theta, omega);  		
	sigma[k] ~ cauchy(0, 2.5);
	
    tau_lambda[k] ~ cauchy(0, 2.5);
    Omega_lambda[k] ~ lkj_corr(2);

	for(q in 1:(Q)){	mu[k,q] ~	normal(0, 10);	}		
	}
	
	
   for (j in 1:N){
		X[j,] ~   multi_normal(mu[studyID[j],], lambda[studyID[j]]);
		Y[j] ~   normal(append_col(1,X[j,])*(beta[studyID[j],]'), sigma[studyID[j]]);
		} 
	
  	}"
  	
  	
  	
  	
  	    	
  	
  	
MEMA1_multi_IPD_kappa <-
"data {
  int<lower=1> N;
  int<lower=1> Q;  			
  int<lower=1> NStudies;		
  int<lower=1,upper=NStudies> studyID[N];
  vector[N] Y;
  matrix[N, 1] X1;
  matrix[N, 1] Xstar2;
  int<lower=1> n_vec[NStudies];
  vector<lower=0>[NStudies] sqrtPsi; 
}

parameters {
  matrix[N, 1] X2;	
  vector<lower=0>[NStudies] sigma;
  matrix[NStudies, Q] mu;
  matrix[NStudies, Q+1] beta;
  vector[Q+1] theta;
  corr_matrix[Q+1] Omega_om;        // prior correlation for omega
  vector<lower=0>[Q+1] tau_om;      // prior scale for omega
  vector[NStudies] kappa0;
  vector[NStudies] kappa1;
  vector<lower=0>[NStudies] kappa2;
}

transformed parameters {
  cov_matrix[Q+1] omega;	
  matrix[N, 2] X;  	
  
  omega = quad_form_diag(Omega_om, tau_om); 
  X = append_col(X1, X2);
 
}

model {
    
    tau_om ~ cauchy(0, 2.5);
    Omega_om ~ lkj_corr(2);

	for(q in 1:(Q+1)){    theta[q]    ~   normal(0, 10);  }
	
    for (k in 1:NStudies){	  
		beta[k,] ~ multi_normal(theta, omega);  		
		sigma[k] ~ cauchy(0, 2.5);
		
		kappa0[k] ~ normal(0, 10);
		kappa1[k] ~ normal(0, 10);	 
		kappa2[k] ~ cauchy(0, 2.5);

		for(q in 1:(Q)){	mu[k,q] ~	normal(0, 10);	}		
	}
	
	
   for (j in 1:N){
   		Xstar2[j,1] ~ normal(X2[j,1], sqrtPsi[studyID[j]]);
   		X2[j,1] ~ normal(kappa0[studyID[j]] +kappa1[studyID[j]]*(X1[j,1]), kappa2[studyID[j]]);
		Y[j] ~   normal(append_col(1, X[j,])*(beta[studyID[j],]'), sigma[studyID[j]]);
	} 
	
  	}"


##############################################################################
#############



  	
MEMA2_multi_IPD_kappa <-
"data {
  int<lower=1> N;
  int<lower=1> Q;  			
  int<lower=1> NStudies;		
  int<lower=1,upper=NStudies> studyID[N];
  vector[N] Y;
  matrix[N, 1] X1;
  matrix[N, 1] Xstar2;
  int<lower=1> n_vec[NStudies];
  real<lower=0> a; 
  real<lower=0> b; 
}

parameters {
  matrix[N, 1] X2;	
  vector<lower=0>[NStudies] sigma;
  matrix[NStudies, Q] mu;
  matrix[NStudies, Q+1] beta;
  vector[Q+1] theta;
  corr_matrix[Q+1] Omega_om;        // prior correlation for omega
  vector<lower=0>[Q+1] tau_om;      // prior scale for omega
  vector[NStudies] kappa0;
  vector[NStudies] kappa1;
  vector<lower=0>[NStudies] kappa2;
  vector<lower=0>[NStudies] sqrtPsi;
}

transformed parameters {
  cov_matrix[Q+1] omega;	
  matrix[N, 2] X;  	
  
  omega = quad_form_diag(Omega_om, tau_om); 
  X = append_col(X1, X2);
 
}

model {
    
    tau_om ~ cauchy(0, 2.5);
    Omega_om ~ lkj_corr(2);

	for(q in 1:(Q+1)){    theta[q]    ~   normal(0, 10);  }
	
    for (k in 1:NStudies){	  
		sqrtPsi[k] ~ 	inv_gamma((a^2)/b + 2, (a^3)/b + a);	
		beta[k,] ~ multi_normal(theta, omega);  		
		sigma[k] ~ cauchy(0, 2.5);

		kappa0[k] ~ normal(0, 10);
		kappa1[k] ~ normal(0, 10);	 
		kappa2[k] ~ cauchy(0, 2.5);

		for(q in 1:(Q)){	mu[k,q] ~	normal(0, 10);	}		
	}
	
	
   for (j in 1:N){
   		Xstar2[j,1] ~ normal(X2[j,1], sqrtPsi[studyID[j]]);
   		X2[j,1] ~ normal(kappa0[studyID[j]] + kappa1[studyID[j]]*(X1[j,1]), kappa2[studyID[j]]);
		Y[j] ~   normal(append_col(1, X[j,])*(beta[studyID[j],]'), sigma[studyID[j]]);
	} 
	
  	}"


##############################################################################
#############
  	
MEMA3_multi_IPD_kappa <-
"data {
  int<lower=1> N;
  int<lower=1> Q;  			
  int<lower=1> NStudies;		
  int<lower=1,upper=NStudies> studyID[N];
  vector[N] Y;
  matrix[N, 1] X1;
  matrix[N, 1] Xstar2;
  int<lower=1> n_vec[NStudies];
  vector<lower=0>[NStudies] MEupperbounds;
}

parameters {
  matrix[N, 1] X2;	
  vector<lower=0>[NStudies] sigma;
  matrix[NStudies, Q] mu;
  matrix[NStudies, Q+1] beta;
  vector[Q+1] theta;
  corr_matrix[Q+1] Omega_om;        // prior correlation for omega
  vector<lower=0>[Q+1] tau_om;      // prior scale for omega
  vector[NStudies] kappa0;
  vector[NStudies] kappa1;
  vector<lower=0>[NStudies] kappa2;
  vector<lower=0>[NStudies] sqrtPsi; 
}

transformed parameters {
  cov_matrix[Q+1] omega;	
  matrix[N, 2] X;  	
  
  omega = quad_form_diag(Omega_om, tau_om); 
  X = append_col(X1, X2);
 
}

model {
    
    tau_om ~ cauchy(0, 2.5);
    Omega_om ~ lkj_corr(2);

	for(q in 1:(Q+1)){    theta[q]    ~   normal(0, 10);  }
	
    for (k in 1:NStudies){	  
		beta[k,] ~ multi_normal(theta, omega);  		
		sigma[k] ~ cauchy(0, 2.5);
		
		kappa0[k] ~ normal(0, 10);
		kappa1[k] ~ normal(0, 10);	 
		kappa2[k] ~ cauchy(0, 2.5);
		
		sqrtPsi[k] ~ 	uniform(0, MEupperbounds[k]);	

		for(q in 1:(Q)){	mu[k,q] ~	normal(0, 10);	}		
	}
	
	
   for (j in 1:N){
   		Xstar2[j,1] ~ normal(X2[j,1], sqrtPsi[studyID[j]]);
   		X2[j,1] ~ normal(kappa0[studyID[j]] +kappa1[studyID[j]]*(X1[j,1]), kappa2[studyID[j]]);
		Y[j] ~   normal(append_col(1, X[j,])*(beta[studyID[j],]'), sigma[studyID[j]]);
	} 
	
  	}"


##############################################################################
#############

## FREQUENTIST GLS:


#######################


GLSmulti <- function(n, y, X){

I <- length(unlist(n))

## Becker and Wu (2007)  (I think assuming that X is fixed not random)

# We stack the k sample slope vectors and make a blockwise diag matrix...
sigma_hat<-vector()
cov_b<-b_estimator<-list()
for(i in 1:I){
mod<-	lm(y[[i]]~X[[i]][,-1])
sigma_hat[i]<-summary(mod)$sigma
b_estimator[[i]] <- solve(t(X[[i]])%*%X[[i]])%*%t(X[[i]])%*%y[[i]]
cov_b[[i]] <- solve(t(X[[i]])%*%X[[i]])*(sigma_hat[i]^2)
}


b_estimator_vertical<-t(matrix(unlist(b_estimator),1,))

SIGMA_matrix <- bdiag(cov_b)

# Regardless of the components of W and beta, we estimate beta and its covariance as

I3by3<-matrix(c(1,0,0,0,1,0,0,0,1),3,3)
el <- (I3by3)
dups <- list(el)[rep(1,I)]
W <- do.call(rbind, dups)


beta_hat_star <- solve(t(W)%*%solve(SIGMA_matrix)%*%W)%*%t(W)%*%solve(SIGMA_matrix)%*% b_estimator_vertical
COV_beta_star_hat <- solve(t(W)%*%solve(SIGMA_matrix)%*%W)

# RESULTS

theta_2 <- (beta_hat_star)[2]

CI_theta_2 <- c((beta_hat_star[2] - qnorm(1-0.025)*sqrt(COV_beta_star_hat[2,2])),(beta_hat_star[2] + qnorm(1-0.025)*sqrt(COV_beta_star_hat[2,2])))

Tau_22 <- sqrt(COV_beta_star_hat[2,2])

theta_1 <- (beta_hat_star)[1]

Tau_11 <- sqrt(COV_beta_star_hat[1,1])


theta_3 <-(beta_hat_star)[3]

Tau_33 <- sqrt(COV_beta_star_hat[3,3])

return( round(c(theta_2=theta_2,
			CI_theta_2= CI_theta_2,
			sqrtTau_22= (Tau_22),
			theta_1= theta_1,
			sqrtTau_11= (Tau_11),
			theta_3= theta_3,
			sqrtTau_33= (Tau_33)),2))
}


