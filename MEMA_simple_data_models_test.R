
########################################################################################
########################################################################################


########################################################################################
######################## DATASETS ########################
########################################################################################


######################## original schooldata  (NELS88) #################################
### dataset :  NELS88
schooldata <- read.csv("~/Desktop/UBC/RECODID_ZIKV/Rcode/13schools.csv")

n <- X <- y <- list()
beta <- alpha <- alpha_hat <- beta_hat <- sigma_hat <- lambda_hat <- mu_hat <- w <- u <- vector()
K <- length(table(schooldata$sch))

for(i in 1:K){
	mydat 			<- schooldata[schooldata$sch==i,]
	n[[i]] 			<- dim(mydat)[1]
	y[[i]] 			<- mydat$sci
	X[[i]] 			<- mydat$rdg
	mod_i 			<- lm(y[[i]] ~ X[[i]])

	alpha_hat[i] 	<- coef(mod_i)[1]
	beta_hat[i] 		<- coef(mod_i)[2]
	sigma_hat[i] 	<- summary(mod_i)$sigma
	lambda_hat[i] 	<- sqrt(var(X[[i]]))
	mu_hat[i] 		<- mean(X[[i]])
	}
	
studyID <- schooldata$sch

NELS88_dataframe <- data.frame(
	"n_i"		=	paste(as.numeric(table(studyID))),
	"alpha"		= 	alpha_hat,
	"beta"		= 	beta_hat, 
	"sigma2"	= 	sigma_hat^2, 
	"mu"			= 	mu_hat, 
	"lambda"	= 	lambda_hat)	
	

NELS88 <-list (
  NStudies 		= length(table(schooldata$sch)),		
  n_per_study 	= as.numeric(table(schooldata$sch)),
  alpha_hat 	= alpha_hat,
  beta_hat 		= beta_hat,
  mu_hat 		= mu_hat,
  sigma_hat 	= sigma_hat,
  lambda_hat 	= lambda_hat)

rm(alpha, alpha_hat, beta, beta_hat, i, K, lambda_hat, mod_i, mu_hat, mydat, n, NELS88_dataframe, schooldata, sigma_hat, studyID, u, w, X, y)


######################## tainted schooldata  (NELS88_star) #################################
### dataset :  NELS88_star
schooldata 		<- read.csv("~/Desktop/UBC/RECODID_ZIKV/Rcode/13schools.csv")
schooldata_me 	<- schooldata

set.seed(123)
tau <- rep(0,NELS88$NStudies)
for(k in 1:3){
tau[k] <- 2.75+k/10
}

for(k in 4:13){
tau[k] <- k/100
}


tausquared 	<- tau^2
gamma 		<- ((1+(tau^2/NELS88$lambda_hat^2))^(-1))

## average, variance, and range attenuation factor is:
# round(range(gamma), 2)
# round(mean(gamma), 2)
# round(var(gamma), 2)

true_theta <- 0.57
true_omega <- 0.05


# predicted value of theta_star:
mean(gamma)* true_theta

# predicted value of omega_star:
((mean(gamma))^2)*true_omega^2 + var(gamma)*((true_omega^2) + true_theta^2)


K <- length(table(schooldata_me$sch))
	for(i in 1:K){
		schooldata_me[schooldata_me$sch==i,]$rdg <- schooldata[schooldata$sch==i,]$rdg +
				rnorm(length(schooldata[schooldata$sch==i,]$rdg),0, tau[[i]])
				}

n_me <- X_me <- y_me <- list()
beta_me <- alpha_me <- alpha_hat_me <- beta_hat_me <- sigma_hat_me <- lambda_hat_me <- mu_hat_me <- w <- u <- vector()

for(i in 1:K){
	mydat_me 		<- schooldata_me[schooldata_me$sch==i,]
	n_me[[i]] 		<- dim(mydat_me)[1]
	y_me[[i]] 		<- mydat_me$sci
	X_me[[i]] 		<- mydat_me$rdg
	mod_i_me 		<- lm(y_me[[i]] ~ X_me[[i]])

	alpha_hat_me[i] <- coef(mod_i_me)[1]
	beta_hat_me[i] 	<- coef(mod_i_me)[2]
	sigma_hat_me[i] <- summary(mod_i_me)$sigma
	lambda_hat_me[i]<- sqrt(var(X_me[[i]]))
	mu_hat_me[i] 	<- mean(X_me[[i]])
	}
	
studyID <- schooldata_me$sch

NELS88star_dataframe <- data.frame(
	"n_i"		= paste(as.numeric(table(studyID))), 
	"tau"		= tau, 
	"gamma" 	= gamma, 
	"alpha"		= alpha_hat_me,
	"beta"		= beta_hat_me, 
	"sigma"		= sigma_hat_me, 
	"mu"			= mu_hat_me, 
	"lambda"	= lambda_hat_me)	
	

NELS88star <-list (
  NStudies 		= length(table(schooldata_me$sch)),		
  n_per_study 	= as.numeric(table(schooldata_me$sch)),
  alpha_hat 	= alpha_hat_me,
  beta_hat 		= beta_hat_me,
  mu_hat 		= mu_hat_me,
  sigma_hat 	= sigma_hat_me,
  lambda_hat 	= lambda_hat_me,
  tau		= tau,
  gamma		= gamma)



rm(alpha_hat_me, alpha_me, beta_hat_me, beta_me, gamma, i, k, K, lambda_hat_me,
mod_i_me, mu_hat_me, mydat_me, n_me, NELS88star_dataframe, schooldata, schooldata_me,
sigma_hat_me, studyID, tau, tausquared, true_omega, true_theta, u, w, X_me, y_me)


# a <- mean(tausquared); b <- var(tausquared) 
# myalpha <- (a^2)/b + 2
# mybeta <- (a^3)/b + a
# rinvg <- rinvgamma(9000000, myalpha, mybeta)
# mean(rinvg); mybeta/(myalpha-1); a
# var(rinvg); (mybeta^2)/(((myalpha-1)^2)*(myalpha-2)); b

#########################################################




########################################################################################
######################## MODELS ########################
########################################################################################



# FREQUENTIST:


REML_model<-function(alpha_hat, beta_hat, mu_hat, sigma_hat, lambda_hat, n, iter=1000, alpha_sig=0.05){


n<-as.list(n)

# intial values:
omega_hat <- sqrt(var(beta_hat)) 
theta_hat <- mean(beta_hat)
alpha_mu_hat <- mean(alpha_hat)
alpha_sd_hat <- sqrt(var(alpha_hat))

# iterative approach to estimate theta and omega with REML estimators
# see Chung, Rabe-Hesketh and Choi (SiM, 2012) :

I<-length(alpha_hat)
print(I)
w <- u <- rep(0, I)
for(jjj in 1:iter){

	for(i in 1:I){
	w[i] <- (omega_hat^2 + 
				(sigma_hat[i]^2/(n[[i]]*lambda_hat[i]^2) ) )^(-1)
	
	u[i] <- ( alpha_sd_hat^2 +
			(((lambda_hat[i]^2+mu_hat[i]^2)*sigma_hat[i]^2)/
			(n[[i]]*lambda_hat[i]^2) )  )^(-1)			
				
				
				}

alpha_sdvarB <- ( ((lambda_hat^2+mu_hat^2)*sigma_hat^2) /	(unlist(n)*lambda_hat^2) ) 
	
alpha_mu_hat <- sum(u*alpha_hat)/sum(u)

alpha_sd_hat <- sqrt(max(c(0, (
					sum((u^2)*
					( (alpha_hat-alpha_mu_hat)^2 - alpha_sdvarB)
						)/ (sum(u^2))	+				
					   (1/sum(u)) ) )))
	
	
	
omegavarB <-	(sigma_hat^2/(unlist(n)*lambda_hat^2))

theta_hat <- sum(w*beta_hat)/sum(w)	
	
omega_hat <- sqrt(max(c(0, (
				sum((w^2)*
				( (beta_hat-theta_hat)^2 - omegavarB)
				         )/  (sum(w^2)) +
				 		(1/sum(w)) ) )))

}

theta_CImargin <-  qnorm(1-alpha_sig/2) * ( sum((omega_hat^2 + ((sigma_hat^2)/(unlist(n)*lambda_hat^2)) )^(-1)) )^(-0.5)


Zstat <- theta_hat / (( sum((omega_hat^2 + ((sigma_hat^2)/(unlist(n)*lambda_hat^2)) )^(-1)) )^(-0.5))

pval <- 2*pnorm(-abs(Zstat))

summarydat <- c(theta_hat, theta_hat-theta_CImargin, theta_hat+theta_CImargin, pval, omega_hat, alpha_mu_hat, alpha_sd_hat) 
names(summarydat)<- c("theta", "theta_CI_2.5%", "theta_CI_97.5%", "H0", "omega", "alpha_mu","alpha_sd") 
return(summarydat)
}



############
BayesMA <- 
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
############

############
BayesMA_H0 <- 
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
########


BayesMEMA1 <- 
"data {
  int<lower=0> NStudies;		
  vector<lower=0>[NStudies]  n_per_study;
  vector[NStudies] alpha_hat_star;
  vector[NStudies] beta_hat_star;
  vector[NStudies] mu_hat_star;
  vector<lower=0>[NStudies] sigma_hat_star;
  vector<lower=0>[NStudies] lambda_hat_star;
  vector<lower=0>[NStudies] tau;  
}

transformed data {      
  vector[NStudies] gamma_hat;
  matrix[NStudies,2] coef_vec_hat;
  vector<lower=0>[NStudies] tausquared;  
  
  for (k in 1:NStudies)
	tausquared[k] = tau[k]^2;
  for (k in 1:NStudies)
  	gamma_hat[k] = ( 1+ tausquared[k]/((lambda_hat_star[k]^2) - tausquared[k]))^(-1);
  for (k in 1:NStudies)
    coef_vec_hat[k,1] = alpha_hat_star[k];	
  for (k in 1:NStudies)
    coef_vec_hat[k,2] = beta_hat_star[k];	
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
  vector[NStudies] mu_hat;
  vector[NStudies] sigma_hat;
  vector[NStudies] lambda_hat;
  vector[NStudies] var_alphahat_star;
  vector[NStudies] var_betahat_star;	 
  vector[NStudies] alpha_star;
  vector[NStudies] beta_star;
  matrix[NStudies, 2] alphabeta_star;	
  matrix[2, 2] SIG_star[NStudies];

  for (k in 1:NStudies)  
	mu_hat[k] = mu_hat_star[k];

  for (k in 1:NStudies)  
     lambda_hat[k] = sqrt(lambda_hat_star[k]^2 - tausquared[k]);
  
  for (k in 1:NStudies)  
	sigma_hat[k] = sqrt((sigma_hat_star[k])^2 - (1 - gamma_hat[k])*(beta[k]^2)*(lambda_hat[k]^2))  ;
    
  for (k in 1:NStudies)
	var_alphahat_star[k] = ((lambda_hat_star[k]^2 + mu_hat_star[k]^2) * sigma_hat_star[k]^2) / ((lambda_hat_star[k]^2) * n_per_study[k]);	

  for (k in 1:NStudies)
	var_betahat_star[k] = (sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);	

  for (k in 1:NStudies)
  	SIG_star[k][1,1] = var_alphahat_star[k];

  for (k in 1:NStudies)
  	SIG_star[k][2,2] = var_betahat_star[k];

  for (k in 1:NStudies)
  	SIG_star[k][1,2] = - mu_hat_star[k]*(sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);  	

  for (k in 1:NStudies)
  	SIG_star[k][2,1] = - mu_hat_star[k]*(sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);  	

  for (k in 1:NStudies)
    alpha_star[k] = alpha[k] + (1-gamma_hat[k])*beta[k]*mu_hat[k];

  for (k in 1:NStudies)
    beta_star[k] = gamma_hat[k]*beta[k];

  alphabeta_star = append_col(alpha_star, beta_star);	

}


model {
  theta ~ normal(0, 10);  
  omega ~ cauchy(0, 25);
  alpha_mu ~ normal(0, 10);  
  alpha_sd ~ cauchy(0, 25);
  for (k in 1:NStudies)	
	coef_vec_hat[k,] ~ multi_normal(alphabeta_star[k,], SIG_star[k]);
  for (k in 1:NStudies)
	alpha[k] ~ normal(alpha_mu, alpha_sd);
  for (k in 1:NStudies)
	beta[k] ~ normal(theta, omega);

  	}
"  	

#######

BayesMEMA2 <- 
"data {
  int<lower=0> NStudies;		
  vector<lower=0>[NStudies]  n_per_study;
  vector[NStudies] alpha_hat_star;
  vector[NStudies] beta_hat_star;
  vector[NStudies] mu_hat_star;
  vector<lower=0>[NStudies] sigma_hat_star;
  vector<lower=0>[NStudies] lambda_hat_star;
  real<lower=0> a;
  real<lower=0> b;
}

transformed data {      
  matrix[NStudies,2] coef_vec_hat;

  for (k in 1:NStudies)
    coef_vec_hat[k,1] = alpha_hat_star[k];	
  for (k in 1:NStudies)
    coef_vec_hat[k,2] = beta_hat_star[k];	
}



parameters {
  vector<lower=0>[NStudies] tau;		
  vector[NStudies] alpha;
  vector[NStudies] beta;
  real theta;
  real<lower=0> omega;
  real alpha_mu;
  real<lower=0> alpha_sd;
}

transformed parameters {
  vector[NStudies]    gamma;	
  vector[NStudies]    mu_hat;
  vector[NStudies]    sigma_hat;
  vector[NStudies]    lambda_hat;
  vector[NStudies]    var_alphahat_star;
  vector[NStudies]    var_betahat_star;	
  matrix[2, 2]        SIG_star[NStudies];  
  vector[NStudies]    alpha_star;
  vector[NStudies]    beta_star;
  matrix[NStudies, 2] alphabeta_star;	
  vector<lower=0>[NStudies] tausquared;  
  
  for (k in 1:NStudies)
	tausquared[k] = tau[k]^2;
  	
  for (k in 1:NStudies)
  	gamma[k] = ( 1+ tausquared[k]/((lambda_hat_star[k]^2) - tausquared[k]))^(-1);
 
   for (k in 1:NStudies)  
	mu_hat[k] = mu_hat_star[k];

  for (k in 1:NStudies)  
     lambda_hat[k] = sqrt(lambda_hat_star[k]^2 - tausquared[k]);
  
  for (k in 1:NStudies)  
	sigma_hat[k] = sqrt((sigma_hat_star[k])^2 - (1 - gamma[k])*(beta[k]^2)*(lambda_hat[k]^2))  ;
  
  for (k in 1:NStudies)
	var_alphahat_star[k] = ((lambda_hat_star[k]^2 + mu_hat_star[k]^2)* sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);	

  for (k in 1:NStudies)
	var_betahat_star[k] = (sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);	

  for (k in 1:NStudies)
  	SIG_star[k][1,1] = var_alphahat_star[k];

  for (k in 1:NStudies)
  	SIG_star[k][2,2] = var_betahat_star[k];

  for (k in 1:NStudies)
  	SIG_star[k][1,2] = - mu_hat_star[k]*(sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);  	

  for (k in 1:NStudies)
  	SIG_star[k][2,1] = - mu_hat_star[k]*(sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);  	


  for (k in 1:NStudies)
    alpha_star[k] = alpha[k] + (1-gamma[k])*beta[k]*mu_hat[k];

  for (k in 1:NStudies)
    beta_star[k] = gamma[k]*beta[k];

  alphabeta_star = append_col(alpha_star, beta_star);	

}


model {
  theta ~ normal(0, 10);  
  omega ~ cauchy(0, 25);
  alpha_mu ~ normal(0, 10);  
  alpha_sd ~ cauchy(0, 25);

  for (k in 1:NStudies)
	tau[k] ~ inv_gamma((a^2)/b + 2, (a^3)/b + a);
  
  for (k in 1:NStudies)	
	coef_vec_hat[k,] ~ multi_normal(alphabeta_star[k,], SIG_star[k]);
  for (k in 1:NStudies)
	alpha[k] ~ normal(alpha_mu, alpha_sd);
  for (k in 1:NStudies)
	beta[k] ~ normal(theta, omega);

  	}
"  	



#######

BayesMEMA3 <- 
"data {
  int<lower=0> NStudies;		
  vector<lower=0>[NStudies]  n_per_study;
  vector[NStudies] alpha_hat_star;
  vector[NStudies] beta_hat_star;
  vector[NStudies] mu_hat_star;
  vector<lower=0>[NStudies] sigma_hat_star;
  vector<lower=0>[NStudies] lambda_hat_star;
  real<lower=0> a_min;
  real<lower=0> a_max;
  real<lower=0> b_min;
  real<lower=0> b_max;

}

transformed data {
  matrix[NStudies,2] coef_vec_hat;

  for (k in 1:NStudies)
    coef_vec_hat[k,1] = alpha_hat_star[k];	
  for (k in 1:NStudies)
    coef_vec_hat[k,2] = beta_hat_star[k];	
}



parameters {
  real<lower=0> a;
  real<lower=0> b;
  vector<lower=0>[NStudies] tau;		
  vector[NStudies] alpha;
  vector[NStudies] beta;
  real theta;
  real<lower=0> omega;
  real alpha_mu;
  real<lower=0> alpha_sd;
}

transformed parameters {
  vector[NStudies]    gamma;	
  vector[NStudies]    mu_hat;
  vector[NStudies]    sigma_hat;
  vector[NStudies]    lambda_hat;
  vector[NStudies]    var_alphahat_star;
  vector[NStudies]    var_betahat_star;	
  matrix[2, 2]        SIG_star[NStudies];  
  vector[NStudies]    alpha_star;
  vector[NStudies]    beta_star;
  matrix[NStudies, 2] alphabeta_star;	
  vector<lower=0>[NStudies] tausquared;  
  
  for (k in 1:NStudies)
	tausquared[k] = tau[k]^2;

  for (k in 1:NStudies)
  	gamma[k] = ( 1+ tausquared[k]/((lambda_hat_star[k]^2) - tausquared[k]))^(-1);

  for (k in 1:NStudies)  
	mu_hat[k] = mu_hat_star[k];
  
   for (k in 1:NStudies)  
     lambda_hat[k] = sqrt(lambda_hat_star[k]^2 - tausquared[k]);
  
  for (k in 1:NStudies)  
	sigma_hat[k] = sqrt((sigma_hat_star[k])^2 - (1 - gamma[k])*(beta[k]^2)*(lambda_hat[k]^2))  ;
  
  for (k in 1:NStudies)
	var_alphahat_star[k] = ((lambda_hat_star[k]^2 + mu_hat_star[k]^2)* sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);	

  for (k in 1:NStudies)
	var_betahat_star[k] = (sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);	

  for (k in 1:NStudies)
  	SIG_star[k][1,1] = var_alphahat_star[k];

  for (k in 1:NStudies)
  	SIG_star[k][2,2] = var_betahat_star[k];

  for (k in 1:NStudies)
  	SIG_star[k][1,2] = - mu_hat_star[k]*(sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);  	

  for (k in 1:NStudies)
  	SIG_star[k][2,1] = - mu_hat_star[k]*(sigma_hat_star[k]^2)/((lambda_hat_star[k]^2)* n_per_study[k]);  	


  for (k in 1:NStudies)
    alpha_star[k] = alpha[k] + (1-gamma[k])*beta[k]*mu_hat[k];

  for (k in 1:NStudies)
    beta_star[k] = gamma[k]*beta[k];

  alphabeta_star = append_col(alpha_star, beta_star);	



}


model {

  a ~ uniform(a_min, a_max);
  b ~ uniform(b_min, b_max);
  
  theta ~ normal(0, 10);  
  omega ~ cauchy(0, 25);
  alpha_mu ~ normal(0, 10);  
  alpha_sd ~ cauchy(0, 25);   
    
  for (k in 1:NStudies)
	tau[k] ~ inv_gamma((a^2)/b + 2, (a^3)/b + a);

  for (k in 1:NStudies)	
	coef_vec_hat[k,] ~ multi_normal(alphabeta_star[k,], SIG_star[k]);
  for (k in 1:NStudies)
	alpha[k] ~ normal(alpha_mu, alpha_sd);
  for (k in 1:NStudies)
	beta[k] ~ normal(theta, omega);

  	}
"  	
