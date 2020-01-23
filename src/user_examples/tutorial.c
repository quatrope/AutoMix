#include "automix.h"
#include "float.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int nsamples = 10;
double data_samples[] = {0.2,  0.13, 0.35, 0.17, 0.89,
                         0.33, 0.78, 0.23, 0.54, 0.16};

double logposterior(int model_k, int mdim, double *theta);
void get_rwm_init(int k, int mdim, double *rwm);

int main() {

  // Initialize the AutoMix sampler
  int nmodels = 3;
  int model_dims[] = {2, 2, 2};
  amSampler am;
  double **initRWM = (double **)malloc(3 * sizeof(double *));
  double initRWM1[2] = {0.5, 0.5};
  double initRWM2[2] = {2.0, 2.0};
  double initRWM3[2] = {9.0, 2.0};
  initRWM[0] = initRWM1;
  initRWM[1] = initRWM2;
  initRWM[2] = initRWM3;
  initAMSampler(&am, nmodels, model_dims, logposterior, initRWM);
  free(initRWM);
  int nsweeps_cpest = 100000; // 100,000
  estimate_conditional_probs(&am, nsweeps_cpest);
  int nburn = 10000; // 10,000
  burn_samples(&am, nburn);
  int nsweeps = 100000; // 100,000
  rjmcmc_samples(&am, nsweeps);

  printf("Probability of model 1: %lf\n", (double)am.st.ksummary[0] / nsweeps);
  printf("Probability of model 2: %lf\n", (double)am.st.ksummary[1] / nsweeps);
  printf("Probability of model 3: %lf\n", (double)am.st.ksummary[2] / nsweeps);

  freeAMSampler(&am);
  return EXIT_SUCCESS;
}

double logp_normal(double sigma, double x0) {
  double prod = 0;
  for (int i = 0; i < nsamples; i++) {
    double x = data_samples[i];
    prod += -(x - x0) * (x - x0);
  }
  prod = -nsamples * log(sigma) + prod / (2.0 * sigma * sigma);
  return prod;
}

double logp_beta(double alpha, double beta) {
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

double logp_gamma(double alpha, double beta) {
  if (alpha <= 0.0 || beta <= 0.0) {
    return -DBL_MAX;
  }
  double prod = 0.0;
  for (int i = 0; i < nsamples; ++i) {
    double x = data_samples[i];
    prod += (alpha - 1.0) * log(x) - beta * x;
  }
  prod += nsamples * (alpha * log(beta) - loggamma(alpha));
  return prod;
}

double logposterior(int model_k, int mdim, double *theta) {
  if (model_k == 0) {
    return logp_normal(theta[0], theta[1]);
  } else if (model_k == 1) {
    return logp_beta(theta[0], theta[1]);
  } else if (model_k == 2) {
    return logp_gamma(theta[0], theta[1]);
  }
  return 0.0;
}
