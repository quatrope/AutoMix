#include "automix.h"
#include "float.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// loggamma is in automix.c but not in the header
extern double loggamma(double x);

double logp_truncnormal_sampler(int model_k, int mdim, double *xp);
double logp_normal_sampler(int model_k, int mdim, double *xp);
double logp_beta_sampler(int model_k, int mdim, double *xp);
double logp_normal_params(int model_k, int mdim, double *params);
double logp_beta_params(int model_k, int mdim, double *params);
double logp_gamma_params(int model_k, int mdim, double *params);
double logp_gamma_beta(int model_k, int mdim, double *params);
double logp_normal_beta(int model_k, int mdim, double *params);
double logp_normal_gamma(int model_k, int mdim, double *params);
void init_normal_sampler(int model_k, int mdim, double *xp);
void init_truncnormal_sampler(int model_k, int mdim, double *xp);
void init_beta_sampler(int model_k, int mdim, double *xp);
void init_normal_params(int model_k, int mdim, double *xp);
void init_beta_params(int model_k, int mdim, double *xp);
void init_gamma_params(int model_k, int mdim, double *xp);
void init_gamma_beta(int model_k, int mdim, double *xp);
void init_normal_gamma(int model_k, int mdim, double *xp);
void init_normal_beta(int model_k, int mdim, double *xp);

int test_setUp(amSampler *am, int models, int *model_dims,
               targetFunc logposterior, double **initRWM);
int test_sampler(amSampler *am, double true_mean, double true_sigma,
                 double lower, double upper);
int test_dist_params(amSampler *am, double true_param1, double true_param2);
int test_two_models(amSampler *am, double true_k1_p1, double true_k1_p2,
                    double true_k2_p1, double true_k2_p2, double true_k1_frac);

int nsamples = 10;
double data_samples[] = {0.50613293, 0.70961096, 0.28166951, 0.12532996,
                         0.46374168, 0.58337466, 0.52458217, 0.56052633,
                         0.57215576, 0.68698825};

int main(int argc, char *argv[]) {
  int pass = EXIT_SUCCESS;
  int model_dim;
  int model_dims[2];
  amSampler am;
  double **initRWM;
  initRWM = (double **)malloc(2 * sizeof(double *));
  initRWM[0] = (double *)malloc(2 * sizeof(double));
  initRWM[1] = (double *)malloc(2 * sizeof(double));

  printf("Test Normal Sampler: . . .");
  model_dim = 1;
  initRWM[0][0] = 0.5;
  test_setUp(&am, 1, &model_dim, logp_normal_sampler, initRWM);
  pass |= test_sampler(&am, 0.5, 1.0, -DBL_MAX, DBL_MAX);
  freeAMSampler(&am);

  printf("Test Truncated Normal Sampler: . . .");
  model_dim = 1;
  initRWM[0][0] = 1.0;
  test_setUp(&am, 1, &model_dim, logp_truncnormal_sampler, initRWM);
  pass |= test_sampler(&am, 1.3, 1.5, 0.0, 10.0);
  freeAMSampler(&am);

  printf("Test Beta Sampler: . . .");
  model_dim = 1;
  initRWM[0][0] = 0.5;
  test_setUp(&am, 1, &model_dim, logp_beta_sampler, initRWM);
  pass |= test_sampler(&am, 0.5, 0.5, 0.0, 1.0);
  freeAMSampler(&am);

  printf("Test Normal Param Estimation: . . .");
  model_dim = 2;
  initRWM[0][0] = 0.5;
  initRWM[0][1] = 0.5;
  test_setUp(&am, 1, &model_dim, logp_normal_params, initRWM);
  pass |= test_dist_params(&am, 0.2, 0.5);
  freeAMSampler(&am);

  printf("Test Beta Param Estimation: . . .");
  model_dim = 2;
  initRWM[0][0] = 2.0;
  initRWM[0][1] = 2.0;
  test_setUp(&am, 1, &model_dim, logp_beta_params, initRWM);
  pass |= test_dist_params(&am, 4.5, 5.0);
  freeAMSampler(&am);

  printf("Test Gamma Param Estimation: . . .");
  model_dim = 2;
  initRWM[0][0] = 9.0;
  initRWM[0][1] = 2.0;
  test_setUp(&am, 1, &model_dim, logp_gamma_params, initRWM);
  pass |= test_dist_params(&am, 7.0, 14.5);
  freeAMSampler(&am);

  printf("Test Gamma-Beta Model Selection: . . .");
  model_dims[0] = 2;
  model_dims[1] = 2;
  initRWM[0][0] = 9.0;
  initRWM[0][1] = 2.0;
  initRWM[1][0] = 2.0;
  initRWM[1][1] = 2.0;
  test_setUp(&am, 2, model_dims, logp_gamma_beta, initRWM);
  pass |= test_two_models(&am, 7.0, 14.5, 4.7, 4.8, 0.37);
  freeAMSampler(&am);

  printf("Test Normal-Beta Model Selection: . . .");
  model_dims[0] = 2;
  model_dims[1] = 2;
  initRWM[0][0] = 0.5;
  initRWM[0][1] = 0.5;
  initRWM[1][0] = 2.0;
  initRWM[1][1] = 2.0;
  test_setUp(&am, 2, model_dims, logp_normal_beta, initRWM);
  pass |= test_two_models(&am, 0.2, 0.5, 4.7, 4.8, 0.95);
  freeAMSampler(&am);

  printf("Test Normal-Gamma Model Selection: . . .");
  model_dims[0] = 2;
  model_dims[1] = 2;
  initRWM[0][0] = 0.5;
  initRWM[0][1] = 0.5;
  initRWM[1][0] = 9.0;
  initRWM[1][1] = 2.0;
  test_setUp(&am, 2, model_dims, logp_normal_gamma, initRWM);
  pass |= test_two_models(&am, 0.2, 0.5, 7.1, 14.5, 0.97);
  freeAMSampler(&am);

  free(initRWM[0]);
  free(initRWM[1]);
  free(initRWM);

  return pass;
}

int test_setUp(amSampler *amp, int nmodels, int *model_dims,
               targetFunc logposterior, double **initRWM) {
  initAMSampler(amp, nmodels, model_dims, logposterior, initRWM);
  int ncond_prob_sweeps = 100000; // 1E5
  estimate_conditional_probs(amp, ncond_prob_sweeps);
  int nburn_sweeps = 10000; // 1E4
  burn_samples(amp, nburn_sweeps);
  int nsweeps = 100000; // 1E5
  rjmcmc_samples(amp, nsweeps);
  return EXIT_SUCCESS;
}

/************************************************
 *                                               *
 *                 Test Checks                   *
 *                                               *
 ************************************************/

int test_sampler(amSampler *am, double true_mean, double true_sigma,
                 double lower, double upper) {
  int ndraws = 100000;
  double mean = 0.0;
  double sumsq = 0.0;
  double **s_samples = (am->st).theta_summary[0];
  for (int i = 0; i < ndraws; i++) {
    double datum = *(s_samples[i]);
    mean += datum;
    sumsq += datum * datum;
    if (datum >= upper || datum <= lower) {
      int pass = 0; // didn't pass the test.
      printf("FAIL\nValue outside range encountered: %lf\n", datum);
      return !pass;
    }
  }
  mean /= ndraws;
  double sigma = sqrt((sumsq - mean * mean) / (ndraws - 1));
  double tol = 0.5;
  int pass = fabs(mean - true_mean) < tol && fabs(sigma - true_sigma) < tol;
  if (!pass) {
    printf("FAIL\nmean=%lf, sigma = %lf\n", mean, sigma);
  } else {
    printf("OK\n");
  }
  return !pass;
}

int test_dist_params(amSampler *am, double true_param1, double true_param2) {
  int ndraws = 100000;
  double p1_mean = 0.0;
  double p2_mean = 0.0;
  double **s_samples = (am->st).theta_summary[0];
  for (int i = 0; i < ndraws; i++) {
    p1_mean += s_samples[i][0];
    p2_mean += s_samples[i][1];
  }
  p1_mean /= ndraws;
  p2_mean /= ndraws;

  double tol = 0.5;
  int pass =
      fabs(p1_mean - true_param1) < tol && fabs(p2_mean - true_param2) < tol;
  if (!pass) {
    printf("FAIL\nparam1=%lf, param2 = %lf\n", p1_mean, p2_mean);
  } else {
    printf("OK\n");
  }
  return !pass;
}

int test_two_models(amSampler *am, double true_k1_p1, double true_k1_p2,
                    double true_k2_p1, double true_k2_p2, double true_k1_frac) {
  int *k_count = (am->st).ksummary;
  double p1_mean[2], p2_mean[2];
  for (int model_k = 0; model_k < 2; ++model_k) {
    double **s_samples = (am->st).theta_summary[model_k];
    int ndraws = (am->st).theta_summary_len[model_k];
    p1_mean[model_k] = 0.0;
    for (int i = 0; i < ndraws; ++i) {
      p1_mean[model_k] += s_samples[i][0];
      p2_mean[model_k] += s_samples[i][1];
    }
    p1_mean[model_k] /= k_count[model_k];
    p2_mean[model_k] /= k_count[model_k];
  }
  double tol = 0.5;
  int pass1 = fabs(p1_mean[0] - true_k1_p1) < tol &&
              fabs(p2_mean[0] - true_k1_p2) < tol;
  if (!pass1) {
    printf("FAIL\nModel 1: (p1=%lf, p2=%lf)\n", p1_mean[0], p2_mean[0]);
  }
  int pass2 = fabs(p1_mean[1] - true_k2_p1) < tol &&
              fabs(p2_mean[1] - true_k2_p2) < tol;
  if (!pass2) {
    printf("FAIL\nModel 2: (p1=%lf, p2=%lf)\n", p1_mean[1], p2_mean[1]);
  }
  double k1_frac = (double)k_count[0] / (k_count[0] + k_count[1]);
  int pass3 = fabs(k1_frac - true_k1_frac) < tol;
  if (!pass3) {
    printf("FAIL\nModel 1 probability: %lf\n", k1_frac);
  }
  if (pass1 && pass2 && pass3) {
    printf("OK\n");
  }
  return !(pass1 && pass2 && pass3);
}

/************************************************
 *                                               *
 *         Distributions to sample from          *
 *                                               *
 ************************************************/

double logp_truncnormal_sampler(int model_k, int mdim, double *xp) {
  double x = *xp;
  double a = 0.0;
  double b = 10.0;
  if (x <= a || x >= b) {
    return -DBL_MAX;
  }
  double prob;
  double x0 = 1.0;
  double sigma = 1.0;
  prob = -(x - x0) * (x - x0) / (2.0 * sigma * sigma);
  return prob;
}

double logp_normal_sampler(int model_k, int mdim, double *xp) {
  double x = *xp;
  double prob;
  double x0 = 0.5;
  double sigma = 1.0;
  prob = -(x - x0) * (x - x0) / (2.0 * sigma * sigma);
  return prob;
}

double logp_beta_sampler(int model_k, int mdim, double *xp) {
  double x = *xp;
  if (x <= 0.0 || x >= 1.0) {
    return -DBL_MAX;
  }
  double alpha = 2.0;
  double beta = 2.0;
  double prod = (alpha - 1.0) * log(x) + (beta - 1.0) * log(1.0 - x);
  prod += loggamma(alpha + beta) - loggamma(alpha) - loggamma(beta);
  return prod;
}

/************************************************
 *                                               *
 *  Distributions with indeterminate parameters  *
 *                                               *
 ************************************************/

double logp_normal_params(int model_k, int mdim, double *params) {
  double sigma = params[0];
  double x0 = params[1];
  double prod = 0;
  for (int i = 0; i < nsamples; i++) {
    double x = data_samples[i];
    prod += -(x - x0) * (x - x0);
  }
  prod = -nsamples * log(sigma) + prod / (2.0 * sigma * sigma);
  return prod;
}

double logp_beta_params(int model_k, int mdim, double *params) {
  double alpha = params[0];
  double beta = params[1];
  if (alpha <= 0.0 || beta <= 0.0) {
    return -DBL_MAX;
  }
  double prod = 0.0;
  for (int i = 0; i < nsamples; i++) {
    double x = data_samples[i];
    prod += (alpha - 1.0) * log(x) + (beta - 1.0) * log(1.0 - x);
  }
  prod +=
      nsamples * (loggamma(alpha + beta) - loggamma(alpha) - loggamma(beta));
  return prod;
}

double logp_gamma_params(int model_k, int mdim, double *params) {
  double alpha = params[0];
  double beta = params[1];
  double prod = 0.0;
  for (int i = 0; i < nsamples; ++i) {
    double x = data_samples[i];
    prod += (alpha - 1.0) * log(x) - beta * x;
  }
  prod += nsamples * (alpha * log(beta) - loggamma(alpha));
  return prod;
}

/************************************************
 *                                               *
 *              Two-model RJMCMC                 *
 *                                               *
 ************************************************/

double logp_gamma_beta(int model_k, int mdim, double *params) {
  if (model_k == 0) {
    return logp_gamma_params(0, 2, params);
  } else if (model_k == 1) {
    return logp_beta_params(0, 2, params);
  }
  return 0.0;
}

double logp_normal_beta(int model_k, int mdim, double *params) {
  if (model_k == 0) {
    return logp_normal_params(0, 2, params);
  } else if (model_k == 1) {
    return logp_beta_params(0, 2, params);
  }
  return 0.0;
}

double logp_normal_gamma(int model_k, int mdim, double *params) {
  if (model_k == 0) {
    return logp_normal_params(0, 2, params);
  } else if (model_k == 1) {
    return logp_gamma_params(0, 2, params);
  }
  return 0.0;
}
