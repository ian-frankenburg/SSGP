functions{
  vector predator_prey(real t, vector y, vector theta) {
    vector[2] dudvdt;
    real u = y[1];
    real v = y[2];
    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];
    dudvdt[1] = (alpha - beta * v) * u;
    dudvdt[2] = (-gamma + delta * u) * v;
    return dudvdt;
  }
}

data {
  int<lower = 0> N; // num measurements
  array[N] real ts; // measurement times > 0
  vector[2] y0; // initial measured population
  real t0;
  array[N] vector<lower = 0>[2] y; // measured population at measurement times
}

parameters {
  vector<lower=0>[4] theta; // theta = { alpha, beta, gamma, delta }
  vector<lower=0>[2] z0; // initial population
  vector<lower=0>[2] sigma; // error scale
}

transformed parameters {
  // ode_rk45(sho, y0, t0, ts, theta)
  array[N] vector[2] z = ode_rk45(predator_prey, z0, t0, ts, theta);
  // real z[N, 2]
  //   = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
  //                        rep_array(0.0, 0), rep_array(0, 0),
  //                        1e-6, 1e-5, 1e3);
}

model {
  // priors on ODE parameters
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  // observation noise parameter
  sigma ~ lognormal(-1, 1);
  // prior on initial ODE state
  z0 ~ lognormal(log(10), 1);
  for (k in 1:2) {
    y0[k] ~ lognormal(log(z0[k]), sigma[k]);
    y[ , k] ~ lognormal(log(z[, k]), sigma[k]);
  }
}

generated quantities {
  vector[2] y0_rep;
  array[N] vector[2] y_rep;
  for (k in 1:2) {
    y0_rep[k] = lognormal_rng(log(z0[k]), sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(z[n, k]), sigma[k]);
  }
}
