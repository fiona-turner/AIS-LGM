functions {
  matrix L_cov_exp_quad_ARD(vector[] x,
                            real alpha,
                            vector rho,
                            real delta) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha
                      * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta + 3; //add fixed variance term for generated quantities - measurement error, spatial error
    return cholesky_decompose(K);
  }
}
data {
  int<lower=0> N; //47
  int<lower=0> M; //10
  matrix[N,M] y1; // 47 * 10
  vector[M] yobs;
  vector[5] x1[N]; //the prior data as a set of vectors
  matrix[N,5] xx1; //the prior data as a matrix
}
transformed data{ //combine the HadCM3 output and ice core obs
  real delta = 1e-9;
  matrix[N+1,M] y;
  for (n1 in 1:N) y[n1,] = y1[n1,];
  y[N + 1,] = yobs';
}

parameters {
  real intercept1;
  real intercept2; // intercept for the regression at each of 10 sites
  real intercept3;
  real intercept4;
  real intercept5;
  real intercept6;
  real intercept7;
  real intercept8;
  real intercept9;
  real intercept10;
  
  vector[5] beta1;
  vector[5] beta2; // regression coefficients at each of 10 sites.
  vector[5] beta3;
  vector[5] beta4;
  vector[5] beta5;
  vector[5] beta6;
  vector[5] beta7;
  vector[5] beta8;
  vector[5] beta9;
  vector[5] beta10;
  
  vector<lower=0>[5] rho;
  real<lower=0> alpha;
  vector<lower=0>[M] sigma1; // nugget std dev on the emulator
  
  vector[5] xobs; //empty vector for prior var values for yobs
  vector[N+1] eta; //N(0,1) variable, so f = L*eta
}

transformed parameters {
  vector[5] x[N+1];
  vector[N+1] f;
  matrix[N+1,5] xx;
  {
    for (n in 1:N) x[n] = x1[n]; //combine x1 and xobs in to one matrix for emulator
    x[N + 1] = xobs;
    
    for (n in 1:N) xx[n] = xx1[n]; //now combine x1 and xobs for prior mean
    xx[N + 1,] = to_row_vector(xobs);
  }
  {
    matrix[N+1, N+1] L_K = L_cov_exp_quad_ARD(x, alpha, rho, delta);
    f = L_K * eta;
  }
}
model{
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma1 ~ inv_gamma(1,1);
  eta ~ std_normal();
  
  intercept1 ~ normal(0, 100);
  intercept2 ~ normal(0, 100);
  intercept3 ~ normal(0, 100);
  intercept4 ~ normal(0, 100);
  intercept5 ~ normal(0, 100);
  intercept6 ~ normal(0, 100);
  intercept7 ~ normal(0, 100);
  intercept8 ~ normal(0, 100);
  intercept9 ~ normal(0, 100);
  intercept10 ~ normal(0, 100);
  
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  beta3 ~ normal(0, 100);
  beta4 ~ normal(0, 100);
  beta5 ~ normal(0, 100);
  beta6 ~ normal(0, 100);
  beta7 ~ normal(0, 100);
  beta8 ~ normal(0, 100);
  beta9 ~ normal(0, 100);
  beta10 ~ normal(0, 100);
  
  xobs[1] ~ normal(0, 0.5); // prior for xobs - the 5 dimensional coordinate for the observed data
  xobs[2] ~ normal(0, 0.5);
  xobs[3] ~ normal(0, 0.6);
  xobs[4] ~ normal(0, 0.5);
  xobs[5] ~ normal(0, 1);
  
  y[,1] ~ normal(intercept1 + xx*beta1 + f, sigma1[1]);
  y[,2] ~ normal(intercept2 + xx*beta2 + f, sigma1[2]);
  y[,3] ~ normal(intercept3 + xx*beta3 + f, sigma1[3]);
  y[,4] ~ normal(intercept4 + xx*beta4 + f, sigma1[4]);
  y[,5] ~ normal(intercept5 + xx*beta5 + f, sigma1[5]);
  y[,6] ~ normal(intercept6 + xx*beta6 + f, sigma1[6]);
  y[,7] ~ normal(intercept7 + xx*beta7 + f, sigma1[7]);
  y[,8] ~ normal(intercept8 + xx*beta8 + f, sigma1[8]);
  y[,9] ~ normal(intercept9 + xx*beta9 + f, sigma1[9]);
  y[,10] ~ normal(intercept10 + xx*beta10 + f, sigma1[10]);
}
generated quantities {
  vector[M] ypred;
  ypred[1] = normal_rng(intercept1 + xobs'*beta1 + f[N+1], sigma1[1]);
  ypred[2] = normal_rng(intercept2 + xobs'*beta2 + f[N+1], sigma1[2]);
  ypred[3] = normal_rng(intercept3 + xobs'*beta3 + f[N+1], sigma1[3]);
  ypred[4] = normal_rng(intercept4 + xobs'*beta4 + f[N+1], sigma1[4]);
  ypred[5] = normal_rng(intercept5 + xobs'*beta5 + f[N+1], sigma1[5]);
  ypred[6] = normal_rng(intercept6 + xobs'*beta6 + f[N+1], sigma1[6]);
  ypred[7] = normal_rng(intercept7 + xobs'*beta7 + f[N+1], sigma1[7]);
  ypred[8] = normal_rng(intercept8 + xobs'*beta8 + f[N+1], sigma1[8]);
  ypred[9] = normal_rng(intercept9 + xobs'*beta9 + f[N+1], sigma1[9]);
  ypred[10] = normal_rng(intercept10 + xobs'*beta10 + f[N+1], sigma1[10]);
}
