#include "automix.h"
#include "float.h"
#include <math.h>
#include <stdio.h>

// loggamma is in automix.c but not in the header
extern double loggamma(double x);

int nsamples = 10;
double data_samples[] = {0.2,  0.13, 0.35, 0.17, 0.89,
                         0.33, 0.78, 0.23, 0.54, 0.16};

double logposterior(int model_k, double *theta);
void get_rwm_init(int k, int mdim, double *rwm);

int main() {
  int nmodels = 3;
  int model_dims[] = {2, 2, 2};
  double initRWM[] = {0.5, 0.5, 2.0, 2.0, 9.0, 2.0};
  amSampler am;
  initAMSampler(&am, nmodels, model_dims, logposterior, initRWM);
  burn_samples(&am, 10000);
  int nsweeps = 100000;
  rjmcmc_samples(&am, nsweeps);

  printf("p(M=1|E) = %lf\n", (double)am.st.ksummary[0] / nsweeps);
  printf("p(M=2|E) = %lf\n", (double)am.st.ksummary[1] / nsweeps);
  printf("p(M=3|E) = %lf\n", (double)am.st.ksummary[2] / nsweeps);

  freeAMSampler(&am);
  return 0;
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

double logposterior(int model_k, double *theta) {
  if (model_k == 0) {
    return logp_normal(theta[0], theta[1]);
  } else if (model_k == 1) {
    return logp_beta(theta[0], theta[1]);
  } else if (model_k == 2) {
    return logp_gamma(theta[0], theta[1]);
  }
  return 0.0;
}
