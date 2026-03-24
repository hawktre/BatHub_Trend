// Dynamic Multi-Year Occupancy Model
// Following Kery and Royle parameterization
// phi and gamma are direct probabilities (no covariates)
// Flexible mean structure for occupancy and detection via design matrices
// Latent state z is marginalized out using the forward algorithm

data {
  int<lower=1> n_sites;                                                     // number of sites
  int<lower=1> n_years;                                                     // number of years
  int<lower=1> n_visits_max;                                                // maximum number of visits per site per year
  int<lower=1> n_xcovs;                                                     // number of occupancy covariates (including intercept)
  int<lower=1> n_pcovs;                                                     // number of detection covariates (including intercept)

  array[n_sites, n_years] int<lower=0> n_visits;                           // number of visits per site per year
  array[n_sites, n_visits_max, n_years] int<lower=0, upper=1> dets;        // observed detections (0/1)

  matrix[n_sites, n_xcovs] xmat;                                            // occupancy design matrix (time-stable)
  array[n_sites, n_visits_max, n_years] row_vector[n_pcovs] pmat;           // detection design matrix (site x visit x year)

  // Prior hyperparameters
  real alpha_prior_mean;                                                     // prior mean for occupancy coefficients
  real<lower=0> alpha_prior_sigma;                                           // prior sigma for occupancy coefficients
  real beta_prior_mean;                                                      // prior mean for detection coefficients
  real<lower=0> beta_prior_sigma;                                            // prior sigma for detection coefficients
}

parameters {
  vector[n_xcovs] alphas;                                                    // occupancy covariate coefficients (includes intercept)
  vector[n_pcovs] betas;                                                     // detection covariate coefficients (includes intercept)
  vector<lower=0, upper=1>[n_years - 1] phi;                               // survival probability for each year transition
  vector<lower=0, upper=1>[n_years - 1] gamma;                             // colonization probability for each year transition
}

model {
  // Priors
  alphas ~ normal(alpha_prior_mean, alpha_prior_sigma);                      // occupancy coefficients
  betas  ~ normal(beta_prior_mean,  beta_prior_sigma);                       // detection coefficients
  phi    ~ beta(1, 1);                                                       // uniform prior on survival probability
  gamma  ~ beta(1, 1);                                                       // uniform prior on colonization probability

  // Likelihood via forward algorithm
  for (i in 1:n_sites) {

    // --- Year 1 ---
    // Occupancy probability in year 1 from occurrence submodel
    real psi_i1 = inv_logit(xmat[i] * alphas);

    // Log-likelihood of detections in year 1 given site is occupied (z=1)
    real log_p_det_1 = 0;
    for (j in 1:n_visits[i, 1]) {
      log_p_det_1 += bernoulli_lpmf(dets[i, j, 1] | inv_logit(pmat[i, j, 1] * betas));
    }

    // Marginal likelihood for year 1, integrating over z_i1 in {0, 1}:
    // P(y_i1) = psi_i1 * P(y_i1 | z=1) + (1 - psi_i1) * P(y_i1 | z=0)
    // If any detections observed, site must be occupied so P(y|z=0) = 0
    real log_lik_1;
    if (sum(dets[i, 1:n_visits_max, 1]) > 0) {
      log_lik_1 = log(psi_i1) + log_p_det_1;
    } else {
      log_lik_1 = log_sum_exp(log(psi_i1) + log_p_det_1, log1m(psi_i1));
    }

    target += log_lik_1;

    // Initialize forward variable a_i1 = P(z_i1 = 1 | y_i1)
    // This is the posterior probability of occupancy after observing year 1 data
    real log_a_occ   = log(psi_i1)   + log_p_det_1 - log_lik_1;            // log P(z_i1 = 1 | y_i1)
    real log_a_unocc = log1m(psi_i1)               - log_lik_1;            // log P(z_i1 = 0 | y_i1)

    // --- Years 2 to n_years ---
    for (t in 2:n_years) {

      // Kery and Royle transition:
      // P(z_it = 1) = phi_{t-1} * P(z_i,t-1 = 1) + gamma_{t-1} * P(z_i,t-1 = 0)
      // P(z_it = 0) = (1 - phi_{t-1}) * P(z_i,t-1 = 1) + (1 - gamma_{t-1}) * P(z_i,t-1 = 0)
      real log_psi_it    = log_sum_exp(log_a_occ   + log(phi[t-1]),
                                       log_a_unocc  + log(gamma[t-1]));
      real log_1m_psi_it = log_sum_exp(log_a_occ   + log1m(phi[t-1]),
                                       log_a_unocc  + log1m(gamma[t-1]));

      // Log-likelihood of detections in year t given site is occupied (z=1)
      real log_p_det_t = 0;
      for (j in 1:n_visits[i, t]) {
        log_p_det_t += bernoulli_lpmf(dets[i, j, t] | inv_logit(pmat[i, j, t] * betas));
      }

      // Marginal likelihood for year t, integrating over z_it in {0, 1}
      real log_lik_t;
      if (sum(dets[i, 1:n_visits_max, t]) > 0) {
        log_lik_t = log_psi_it + log_p_det_t;
      } else {
        log_lik_t = log_sum_exp(log_psi_it + log_p_det_t, log_1m_psi_it);
      }

      target += log_lik_t;

      // Update forward variable: a_it = P(z_it = 1 | y_i1, ..., y_it)
      log_a_occ   = log_psi_it    + log_p_det_t - log_lik_t;
      log_a_unocc = log_1m_psi_it               - log_lik_t;
    }
  }
}

generated quantities {
  // Site-level marginal occupancy probability per site and year
  matrix[n_sites, n_years] psi;

  // Mean occupancy probability per year averaged across sites
  vector[n_years] mean_psi;

  // Mean detection probability per year averaged across sites and visits
  vector[n_years] mean_p;

  // Year 1: occupancy from occurrence submodel
  {
    real sum_psi = 0;
    real sum_p   = 0;
    int  n_obs   = 0;
    for (i in 1:n_sites) {
      psi[i, 1]  = inv_logit(xmat[i] * alphas);
      sum_psi   += psi[i, 1];
      for (j in 1:n_visits[i, 1]) {
        sum_p  += inv_logit(pmat[i, j, 1] * betas);
        n_obs  += 1;
      }
    }
    mean_psi[1] = sum_psi / n_sites;
    mean_p[1]   = sum_p   / n_obs;
  }

  // Years 2 to n_years: occupancy from Kery and Royle transition
  for (t in 2:n_years) {
    real sum_psi = 0;
    real sum_p   = 0;
    int  n_obs   = 0;
    for (i in 1:n_sites) {
      // Marginal occupancy: P(z_it=1) = phi * P(z_i,t-1=1) + gamma * P(z_i,t-1=0)
      psi[i, t]  = phi[t-1] * psi[i, t-1] + gamma[t-1] * (1 - psi[i, t-1]);
      sum_psi   += psi[i, t];
      for (j in 1:n_visits[i, t]) {
        sum_p  += inv_logit(pmat[i, j, t] * betas);
        n_obs  += 1;
      }
    }
    mean_psi[t] = sum_psi / n_sites;
    mean_p[t]   = sum_p   / n_obs;
  }
}