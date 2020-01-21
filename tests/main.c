#include "automix.h"
#include "logwrite.h"
#include "utils.h"
#include "float.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double logp_truncnormal_sampler(int model_k, int mdim, double *xp);
double logp_normal_sampler(int model_k, int mdim, double *xp);
double logp_beta_sampler(int model_k, int mdim, double *xp);
double logp_normal_params(double sigma, double x0);
double logp_beta_params(double alpha, double beta);
double logp_gamma_params(double alpha, double beta);
void init_normal_sampler(int model_k, int mdim, double *xp);
int test_setUp(int models, int *model_dims, targetFunc logposterior,
               rwmInitFunc initRWM);
int test_tearDown(char *filename, int nmodels);
int test_normal_sampler();

int nsamples = 10;
double data_samples[] = {0.50613293, 0.70961096, 0.28166951, 0.12532996,
                         0.46374168, 0.58337466, 0.52458217, 0.56052633,
                         0.57215576, 0.68698825};

int main(int argc, char *argv[]) {
  int pass = EXIT_SUCCESS;

  printf("Test Normal Sampler: . . .\n");
  int model_dims = 1;
  test_setUp(1, &model_dims, logp_normal_sampler, init_normal_sampler);
  pass |= test_normal_sampler();
  test_tearDown("test", 1);

  return pass;
}

int test_setUp(int nmodels, int *model_dims, targetFunc logposterior,
               rwmInitFunc initRWM) {
  amSampler am;
  initAMSampler(&am, nmodels, model_dims, logposterior, initRWM);
  int ncond_prob_sweeps = 100000; // 1E5
  estimate_conditional_probs(&am, ncond_prob_sweeps);
  report_cond_prob_estimation("test", am);
  int nburn_sweeps = 10000; // 1E4
  burn_samples(&am, nburn_sweeps);
  int nsweeps = 100000; // 1E5
  rjmcmc_samples(&am, nsweeps);
  report_rjmcmc_run("test", am, 0, ncond_prob_sweeps, nsweeps);
  freeAMSampler(&am);
  return EXIT_SUCCESS;
}

int test_tearDown(char *filename, int nmodels) {
  char fname_buffer[50];
  sprintf(fname_buffer, "%s_ac.data", filename);
  remove(fname_buffer);
  sprintf(fname_buffer, "%s_cf.data", filename);
  remove(fname_buffer);
  sprintf(fname_buffer, "%s_log.data", filename);
  remove(fname_buffer);
  sprintf(fname_buffer, "%s_mix.data", filename);
  remove(fname_buffer);
  sprintf(fname_buffer, "%s_adapt.data", filename);
  remove(fname_buffer);
  sprintf(fname_buffer, "%s_k.data", filename);
  remove(fname_buffer);
  sprintf(fname_buffer, "%s_lp.data", filename);
  remove(fname_buffer);
  sprintf(fname_buffer, "%s_pk.data", filename);
  remove(fname_buffer);
  for (int i = 0; i < nmodels; ++i) {
    sprintf(fname_buffer, "%s_theta%d.data", filename, i + 1);
    remove(fname_buffer);
  }
  return EXIT_SUCCESS;
}

/************************************************
 *                                               *
 *                 Test Checks                   *
 *                                               *
 ************************************************/

int test_normal_sampler() {
  FILE *fp = fopen("test_theta1.data", "r");
  if (fp == NULL) {
      return EXIT_FAILURE;
  }
  int ndraws = 100000;
  double mean = 0.0;
  double sumsq = 0.0;
  for (int i = 0; i < ndraws; i++) {
    double datum;
    fscanf(fp, "%lf", &datum);
    mean += datum;
    sumsq += datum * datum;
  }
  fclose(fp);
  mean /= ndraws;
  double sigma = sqrt((sumsq - mean * mean) / (ndraws - 1));
  double true_mean = 0.5;
  double true_sigma = 1.0;
  double tol = 0.5;
  int pass = fabs(mean - true_mean) < tol && fabs(sigma - true_sigma) < tol;
  printf("Test Normal Sampler:......");
  if (!pass) {
    printf("Test didn't pass.\nmean=%lf, sigma = %lf\n", mean, sigma);
  } else {
    printf("OK\n");
  }
  return !pass;
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

double logp_normal_params(double sigma, double x0) {
  double prod = 0;
  for (int i = 0; i < nsamples; i++) {
    double x = data_samples[i];
    prod += -(x - x0) * (x - x0);
  }
  prod = -nsamples * log(sigma) + prod / (2.0 * sigma * sigma);
  return prod;
}

double logp_beta_params(double alpha, double beta) {
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

double logp_gamma_params(double alpha, double beta) {
  double prod = 0.0;
  for (int i = 0; i < nsamples; ++i) {
    double x = data_samples[i];
    prod += (alpha - 1.0) * log(x) - beta * x;
  }
  prod += nsamples * (alpha * log(beta) - loggamma(alpha));
  return prod;
}

void init_normal_sampler(int model_k, int mdim, double *xp) { *xp = 0.5; }