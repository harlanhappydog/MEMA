

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



