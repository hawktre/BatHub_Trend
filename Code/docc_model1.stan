//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions{
  real docc_lpmf(int[] dets,
                 matrix[] logit_psi,
                 vector logit_p,
                 int n_sites,
                 int n_years,
                 int[] n_visits,
                 int[] naive_ind){
    vector[n_sites] out;
    int occ_states[2];
      
    int pos1 = 1;
    int pos2 = 1;
    
    occ_states[1] = 1;
    occ_states[2] = 0;
    
    for(n in 1:n_sites){
      vector[2] forward[n_years];
      vector[n_years] log_p;
      int temp_length[n_years];
      
      int n_visits_temp[n_years] = segment(n_visits, pos1, n_years);
      int naive_temp[n_years] = segment(naive_ind, pos1, n_years);
      
      for(t in 1:n_years){
        int dets_temp[n_visits_temp[t]] = segment(dets, pos2, n_visits_temp[t]);
        vector[n_visits_temp[t]] logit_p_temp = segment(logit_p, pos2, n_visits_temp[t]);
        
        log_p[t] = bernoulli_logit_lpmf(dets_temp | logit_p_temp);
        temp_length[t] = 2 - naive_temp[t];
        
        pos2 += n_visits_temp[t];
      }
      
      for(k in 1:temp_length[1]){
        forward[1, k] = bernoulli_logit_lpmf(occ_states[k] | logit_psi[n][1, k]) +
                          (k == 1) * log_p[1];
      }
      
      for(t in 2:n_years){
        for(k in 1:temp_length[t]){
          real acc[temp_length[(t-1)]];
          for(j in 1:temp_length[(t-1)]){
            acc[j] = forward[(t-1), j] +
                       bernoulli_logit_lpmf(occ_states[k] | logit_psi[n][t, j]) +
                       (k == 1) * log_p[t];
          }
          forward[t, k] = log_sum_exp(acc);
        }
      }
      
      if(naive_temp[n_years] == 1){
        out[n] = forward[n_years, 1];
      } else{
        out[n] = log_sum_exp(forward[n_years]);
      }
      
      pos1 += n_years;
    }
    
    return(sum(out));
    
  }
  
}

data{
  int<lower=1> n_sites_total;  //total number of sites
  int<lower=1> n_site_years;   //number of sites by year
  int<lower=1> n_years;        //number of years
  int<lower=1> n_obs;          //number of observations total (total visits)
  
  int<lower=0, upper=1> dets[n_obs];      //detection histories by species
  int<lower=0> n_visits[n_site_years];    //number of visits per site
  
  //Indicator for whether each site was naively occupied or not, by year
  int<lower=0, upper=1> naive_ind[n_site_years];
  
  //note that the covariate matrices do not include intercepts and therefore
  //assume at least one covariate for both occupancy and detection.
  int<lower=1> n_covs1;                //number of covariates for occupancy
  matrix[n_sites_total, n_covs1] xmat; //(scaled) covariate matrix for occupancy
  int<lower=1> n_covs2;                //number of covariates for detection
  matrix[n_obs, n_covs2] vmat;         //(scaled) covariate matrix for detection
}

parameters{
  //Occupancy parameters
  //Intercepts, different for the first year
  real alpha01;
  real alpha02;
  
  //Coefficients for environmental variables
  vector[n_covs1] alphas;
  
  //Autocorrelation parameter (Markov-structured covariate)
  real alpha_auto;
  
  //Detection parameters, intercept and coefficients
  real beta0;
  vector[n_covs2] betas;
}

transformed parameters{
  //Vectors of logit probabilities, using approximate probit function
  matrix[n_years, 2] logit_psi[n_sites_total];
  vector[n_obs] logit_p;
  
  //Calculate probit and logit probabilities
  for(n in 1:n_sites_total){
    logit_psi[n][1, 1] = alpha01 + xmat[n] * alphas;
    logit_psi[n][1, 2] = alpha01 + xmat[n] * alphas;
    for(t in 2:n_years){
      logit_psi[n][t, 1] = alpha02 + xmat[n] * alphas + alpha_auto;
      logit_psi[n][t, 2] = alpha02 + xmat[n] * alphas;
    }
  }
  
  logit_p = beta0 + vmat * betas;
}

model{
  dets ~ docc(logit_psi, logit_p,
              n_sites_total, n_years, n_visits, naive_ind);
  
  //Priors     
  alpha01 ~ normal(0, 10);
  alpha02 ~ normal(0, 10);
  alphas ~ normal(0, 5);
  alpha_auto ~ normal(0, 5);
  
  beta0 ~ normal(0, 10);
  betas ~ normal(0, 5);
}
