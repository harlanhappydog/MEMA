


######  linear_regression_MA  function ######

linear_regression_MA <- function(y, X, studyID, mcmc_iter=5000){



max_obs_per_study <- 100

Q <- dim(X)[2]
N <- length(y)
K <- length(unique(studyID))

the_data<-data.frame(cbind(y, X, studyID))

sigma_hat <- beta_hat  <- n <- X <- y <- list()

for(i in 1:K){
	
	mydat <- the_data[the_data$studyID==i,]
	n[[i]] <- dim(mydat)[1]
	y[[i]] <- mydat$y
	X[[i]] <- cbind(mydat[,-c(1, dim(the_data)[2])] )
	mod_i <- lm(y[[i]] ~ as.matrix(X[[i]]))
	beta_hat[[i]] <- coef(mod_i)
	sigma_hat[i] <- summary(mod_i)$sigma
	}
	
	
	
indicator <- matrix(0, K, max_obs_per_study)
for(k in 1:K){
ii <- (1:(N+1))[studyID == k]
indicator[k,1:length(ii)] <- ii
indicator[k,][indicator[k,]==0]<-(N+1)
}

multi_XAG <-list (
  N = N,	
  Q = Q,
  NStudies = K,
  studyID = studyID,
  beta_hat = do.call(rbind,beta_hat),
  X = rbind(cbind(1,do.call(rbind, X)),rep(0,Q+1)),
  indicator = indicator,
  sigma = unlist(sigma_hat),
  thetaPrior_mu = rep(0, Q+1),
  thetaPrior_sigma = diag(rep(10, Q+1)) ,
  omegaPrior_nu = Q+1,
  omegaPrior_sigma = diag(rep(25, Q+1)) 
					)
 
N <- NStudies <- studyID <- beta_hat  <- sigma_hat <- n <- X <- y <- indicator <- NULL

    
                 
MA_multivariable_H1 <-
"data {
  int<lower=1> N;
  int<lower=0> Q;  			
  int<lower=1> NStudies;		
  int<lower=1,upper=NStudies> studyID[N];
  matrix[NStudies,Q+1] beta_hat;
  matrix[N+1, Q+1] X;
  int indicator[NStudies, 100];
  vector<lower=0>[NStudies] sigma;
  vector[Q+1] thetaPrior_mu;
  matrix[Q+1, Q+1] thetaPrior_sigma;
  real<lower = Q> omegaPrior_nu;                                                                  
  cov_matrix[Q+1] omegaPrior_sigma;  
}


transformed data {
  matrix[Q+1, Q+1] COV[NStudies];
  
	for(k in 1:NStudies){
		COV[k] = (inverse_spd( (X[indicator[k,], ])' * (X[indicator[k,], ]) ))*sigma[k]^2;
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
  	


                 
MA_multivariable_H0 <-
"data {
  int<lower=1> N;
  int<lower=0> Q;  			
  int<lower=1> NStudies;		
  int<lower=1,upper=NStudies> studyID[N];
  matrix[NStudies,Q+1] beta_hat;
  matrix[N+1, Q+1] X;
  int indicator[NStudies, 100];
  vector<lower=0>[NStudies] sigma;
  vector[Q+1] thetaPrior_mu;
  matrix[Q+1, Q+1] thetaPrior_sigma;
  real<lower = Q> omegaPrior_nu;                                                                  
  cov_matrix[Q+1] omegaPrior_sigma;  
}


transformed data {
  matrix[Q+1, Q+1] COV[NStudies];
  
	for(k in 1:NStudies){
		COV[k] = (inverse_spd( (X[indicator[k,], ])' * (X[indicator[k,], ]) ))*sigma[k]^2;
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
  vector[2] theta_0;
  vector[1] zero;
  vector[Q+1] theta_all;

  zero[1] = 0.00;

  omega = quad_form_diag(Omega_om, tau_om); 
  theta_0 = append_row(theta[1], zero);
  theta_all = append_row(theta_0,theta[3:(Q+1)]);
}

model {


    theta    ~   multi_normal(thetaPrior_mu, thetaPrior_sigma);  
    
    tau_om ~ cauchy(0, 2.5);
    
    Omega_om ~ lkj_corr(2);


  for (k in 1:NStudies)	  
	beta[k,] ~   multi_normal(theta_all, omega);

  for (k in 1:NStudies)	  
	beta_hat[k,] ~   multi_normal(beta[k,], COV[k]);
	
  	}"
  	



init_fn <- function() {
  list(beta = cbind(multi_XAG$beta_hat),
	   theta = colMeans(cbind(multi_XAG$beta_hat)),
       omega = diag(rep(2, multi_XAG$Q+1))
       )
}


MAAD_XIPD_H0 <- stan(model_code=MA_multivariable_H0,
                 iter = mcmc_iter, 
                 warmup = round(mcmc_iter*0.10),
                 data = multi_XAG,
                 control=list(adapt_delta=0.95),
                 chains = 1, 
                 init = init_fn)


MAAD_XIPD_H1 <- stan(model_code=MA_multivariable_H1,
                 iter = mcmc_iter, 
                 warmup = round(mcmc_iter*0.10),
                 data = multi_XAG,
                 control=list(adapt_delta=0.95),
                 chains = 1, 
                 init = init_fn)



bridge_H0 <- bridge_sampler(MAAD_XIPD_H0); bridge_H1 <- bridge_sampler(MAAD_XIPD_H1)
BF_10_result <- bayes_factor(bridge_H1, bridge_H0)
BF_10 <- as.numeric(BF_10_result)[1]                 
                 
fit_MAAD_XIPD <- summary(MAAD_XIPD_H1)

print(Q)

theta_string <- rep("", Q+1)
for(q in 1:(Q+1)){theta_string[q] <- paste("theta[",q,"]", sep="")}

summary_theta <- matrix(0, Q+1,3)
for(q in 1:(Q+1)){
summary_theta[q,] <- c(fit_MAAD_XIPD$summary[theta_string[q],c("mean", "2.5%","97.5%")])
}
rownames(summary_theta) <- theta_string
colnames(summary_theta) <- c("mean", "2.5%","97.5%")
summary_theta


omega_string <- rep("", Q+1)
for(q in 1:(Q+1)){omega_string[q] <- paste("omega[", q, "," , q, "]", sep="")}

summary_omega <- matrix(0, Q+1, 3)
for(q in 1:(Q+1)){
summary_omega[q,] <- c(fit_MAAD_XIPD$summary[omega_string[q],c("mean", "2.5%","97.5%")])}

rownames(summary_omega) <- omega_string
colnames(summary_omega) <- c("mean", "2.5%","97.5%")
summary_omega

analysis_summary<- rbind(summary_theta, summary_omega)

return(list(summary=round(analysis_summary,5), BF_10_for_theta_2=BF_10)) }







##################



######  simple_linear_regression_MA  function ######

simple_linear_regression_MA <- function(y, X, studyID, mcmc_iter=5000){

library(bridgesampling)
library(rstan)

max_obs_per_study <- 100

Q <- dim(X)[2]
N <- length(y)
K <- length(unique(studyID))

the_data<-data.frame(cbind(y, X, studyID))

n <- X <- y <- list()
beta <- alpha <- alpha_hat <- beta_hat <- sigma_hat <- lambda_hat <- mu_hat <- w <- u <- vector()


for(i in 1:K){
	
	mydat <- the_data[the_data$studyID==i,]
	n[[i]] <- dim(mydat)[1]
	y[[i]] <- mydat$y
	X[[i]] <- cbind(mydat[,-c(1, dim(the_data)[2])] )
	mod_i <- lm(y[[i]] ~ as.matrix(X[[i]]))
	alpha_hat[i] <- coef(mod_i)[1]
	beta_hat[i] <- coef(mod_i)[2]
	sigma_hat[i] <- summary(mod_i)$sigma
	lambda_hat[i] <- sqrt(var(X[[i]]))
	mu_hat[i] <- mean(X[[i]])

	}
	
	
	
DATA_list<-list (
  NStudies = K,		
  n_per_study = unlist(n),
  alpha_hat = alpha_hat,
  beta_hat = beta_hat,
  mu_hat = mu_hat,
  sigma_hat = sigma_hat,
  lambda_hat = lambda_hat
)

 




MA_AD_H1 <- 
"data {
  int<lower=0> NStudies;		
  vector<lower=0>[NStudies]  n_per_study;
  vector[NStudies] alpha_hat;
  vector[NStudies] beta_hat;
  vector[NStudies] mu_hat;
  vector<lower=0>[NStudies] sigma_hat;
  vector<lower=0>[NStudies] lambda_hat;
}

transformed data {
  vector[NStudies] var_alphahat;
  vector[NStudies] var_betahat;	
  matrix[NStudies,2] coef_vec_hat;
  matrix[2, 2] SIG[NStudies];
  for (k in 1:NStudies)
	var_alphahat[k] = ((lambda_hat[k]^2 + mu_hat[k]^2)* sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);	
  for (k in 1:NStudies)
	var_betahat[k] = (sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);	

  for (k in 1:NStudies)
  	SIG[k][1,1] = var_alphahat[k];
  for (k in 1:NStudies)
  	SIG[k][2,2] = var_betahat[k];
  for (k in 1:NStudies)
  	SIG[k][1,2] = - mu_hat[k]*(sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);  	
  for (k in 1:NStudies)
  	SIG[k][2,1] = - mu_hat[k]*(sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);  	
print(SIG[1])  
  for (k in 1:NStudies)
coef_vec_hat[k,1] = alpha_hat[k];	
  for (k in 1:NStudies)
coef_vec_hat[k,2] = beta_hat[k];	

}	


parameters {
  vector[NStudies] alpha;
  vector[NStudies] beta;
  real theta;
  real<lower=0> omega;
  real alpha_mu;
  real<lower=0> alpha_sd;
}

transformed parameters {

matrix[NStudies, 2] alphabeta;	
alphabeta = append_col(alpha, beta);	

}


model {
  theta ~ normal(0, 10);  
  omega ~ cauchy(0, 25);
  alpha_mu ~ normal(0, 10);  
  alpha_sd ~ cauchy(0, 25);
  for (k in 1:NStudies)	
	coef_vec_hat[k,] ~ multi_normal(alphabeta[k,], SIG[k]);
  for (k in 1:NStudies)
	alpha[k] ~ normal(alpha_mu, alpha_sd);
  for (k in 1:NStudies)
	beta[k] ~ normal(theta, omega);

  	}
"  	

MA_AD_H0 <- 
"data {
  int<lower=0> NStudies;		
  vector<lower=0>[NStudies]  n_per_study;
  vector[NStudies] alpha_hat;
  vector[NStudies] beta_hat;
  vector[NStudies] mu_hat;
  vector<lower=0>[NStudies] sigma_hat;
  vector<lower=0>[NStudies] lambda_hat;
}

transformed data {
  vector[NStudies] var_alphahat;
  vector[NStudies] var_betahat;	
  matrix[NStudies,2] coef_vec_hat;
  matrix[2, 2] SIG[NStudies];
  for (k in 1:NStudies)
	var_alphahat[k] = ((lambda_hat[k]^2 + mu_hat[k]^2)* sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);	
  for (k in 1:NStudies)
	var_betahat[k] = (sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);	

  for (k in 1:NStudies)
  	SIG[k][1,1] = var_alphahat[k];
  for (k in 1:NStudies)
  	SIG[k][2,2] = var_betahat[k];
  for (k in 1:NStudies)
  	SIG[k][1,2] = - mu_hat[k]*(sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);  	
  for (k in 1:NStudies)
  	SIG[k][2,1] = - mu_hat[k]*(sigma_hat[k]^2)/((lambda_hat[k]^2)* n_per_study[k]);  	
print(SIG[1])  
  for (k in 1:NStudies)
coef_vec_hat[k,1] = alpha_hat[k];	
  for (k in 1:NStudies)
coef_vec_hat[k,2] = beta_hat[k];	

}	


parameters {
  vector[NStudies] alpha;
  vector[NStudies] beta;

  real<lower=0> omega;
  real alpha_mu;
  real<lower=0> alpha_sd;
}

transformed parameters {

matrix[NStudies, 2] alphabeta;	
alphabeta = append_col(alpha, beta);	

}


model {

  omega ~ cauchy(0, 25);
  alpha_mu ~ normal(0, 10);  
  alpha_sd ~ cauchy(0, 25);
  for (k in 1:NStudies)	
	coef_vec_hat[k,] ~ multi_normal(alphabeta[k,], SIG[k]);
  for (k in 1:NStudies)
	alpha[k] ~ normal(alpha_mu, alpha_sd);
  for (k in 1:NStudies)
	beta[k] ~ normal(0, omega);

  	}
"  	


init_fn <- function() {
  list(alpha = DATA_list$alpha_hat,
       beta = DATA_list$beta_hat,
       theta = mean(DATA_list$beta_hat),
       omega = sd(DATA_list$beta_hat),
       alpha_mu = mean(DATA_list$alpha_hat),
       alpha_sd = sd(DATA_list$alpha_hat)
       )
}

mod_H0 <- stan(model_code=MA_AD_H0,
                 iter = mcmc_iter, 
                 warmup = round(mcmc_iter*0.10),
                 data = DATA_list,
                 control=list(adapt_delta=0.95),
                 chains = 1, init = init_fn)


mod_H1 <- stan(model_code=MA_AD_H1,
                 iter = mcmc_iter, 
                 warmup = round(mcmc_iter*0.10),
                 data = DATA_list,
                 control=list(adapt_delta=0.95),
                 chains = 1, init = init_fn)
                 

library("bridgesampling")
bridge_H0 <- bridge_sampler(mod_H0); 
bridge_H1 <- bridge_sampler(mod_H1)


BF_10_result <- bf(bridge_H1, bridge_H0)
BF_10 <- as.numeric(BF_10_result)[1]     

fit_MAAD <- summary(mod_H1)

summary_MAAD <- rbind(
fit_MAAD$summary["alpha_mu",c("mean", "2.5%","97.5%")], 
fit_MAAD$summary["theta",c("mean", "2.5%","97.5%")], 
fit_MAAD$summary["alpha_sd",c("mean", "2.5%","97.5%")],
fit_MAAD$summary["omega",c("mean", "2.5%","97.5%")]
)


rownames(summary_MAAD)<- c("alpha_mu", "theta", "alpha_sd", "omega") 
colnames(summary_MAAD) <- c("mean", "2.5%","97.5%")


analysis_summary<- summary_MAAD

return(list(summary=round(analysis_summary,5), BF_10_for_theta=BF_10)) }

