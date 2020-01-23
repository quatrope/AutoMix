.. _tutorial:

Tutorial
========

Example
-------

Suppose we have 10 samples from a random process X that we know not much about.

.. math::

    X = \left\{ 0.2, 0.13, 0.35, 0.17, 0.89, 0.33, 0.78, 0.23, 0.54, 0.16 \right\}

We will assume the samples are independent, identically distributed (IID) from some unknown probability density distribution (PDF).
The nature of the problem indicates that it could come from a Normal distribution, a Gamma distribution or perhaps a Beta distribution.
We don't have any other information on the parameters,

So we propose three models to account for what we know:

* **Model 1**: Normal distribution with arbitrary parameters.

  .. math::

      p_1(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left(-\frac{(x-x_0)^2}{2 \sigma^2}\right)

  This model has undertermined parameters :math:`\sigma` and :math:`x_0`.

* **Model 2**: Beta distribution

  .. math::

      p_2(x) = \frac{1}{B(\alpha, \beta)} x^{\alpha - 1} (1 - x)^{\beta - 1}

  This model has undetermined parameters :math:`\alpha` and :math:`\beta`
  under the restriction :math:`\alpha > 0` and :math:`\beta > 0`.

* **Model 3**: Gamma distribution

  .. math::

      p_3(x) = \frac{\beta^{\alpha} x^{\alpha - 1} e^{-\beta x}}{\Gamma(\alpha)}

  This model has undetermined parameters :math:`\alpha` and :math:`\beta`
  under the restriction :math:`\alpha > 0` and :math:`\beta > 0`.

.. note::

    The dimensions for the problems are 2 for all 3 models.

We would like to find the probability of each model given the evidence (called posterior in a Bayesian analysis.)

Assuming independent, identically distributed (IID) samples, the probability of the evidence given the model (or likelikhood) can be easily calculated

  .. math::

      p(X|M) = \prod_i p_{M}(x_i)

The probability of the model given the evidence is then

  .. math::

    p(M|X) \propto \prod_i p_{M}(x_i) p(M)

In what follows we will assume all models are equally likely and set :math:`p(M)=1/3 \forall M`.

The AutoMix sampler needs to be provided with a function that returns the target distribution (the posterior in our example),
for each model. The function prototype must be in the form:

.. code-block:: c

  double targetDistribution(int model, int model_dim, double *params);

The log-posteriors need to be implemented up to an additive constant.
A quick implementation of the log-PDFs could be:

.. code-block:: c

  #include "float.h"
  #include "utils.h"
  #include <math.h>
  
  int nsamples = 10;
  double data_samples[] = {0.2,  0.13, 0.35, 0.17, 0.89,
                           0.33, 0.78, 0.23, 0.54, 0.16};
  
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
    double prod = 0.0;
    for (int i = 0; i < nsamples; ++i) {
      double x = data_samples[i];
      prod += (alpha - 1.0) * log(x) - beta * x;
    }
    prod += nsamples * (alpha * log(beta) - loggamma(alpha));
    return prod;
  }
  
  double logposterior(int model, int model_dim, double *params) {
    if (model == 1) {
      return logp_normal(params[0], params[1]);
    } else if (model == 2) {
      return logp_beta(params[0], params[1]);
    } else if (model == 3) {
      return logp_gamma(params[0], params[1]);
    }
    return 0.0;
  }

Now that we have defined our problem and log-posteriors, we can set up AutoMix to generate samples from the posteriors of each model.


AutoMix Sampler
---------------

Below is an exanple of a minimal call to generate 1,000 samples from a given problem:

.. code-block:: c

  #include "automix.h"
  double logposterior(int model, int mdim, double* theta);
  void get_rwm_init(int model, int mdim, double *theta);

  int main() {
    AutoMixSampler am;
    am.setTargetDistribution(logposterior);
    am.setRWMInitFunction(get_rwm_init);
    am.setNumberOfModels(3);
    int model_dims[] = {4, 3, 2};
    am.setModelDimensions(model_dims);
    am.rjmcmc_samples(1000);
    return 0;
  } 

This is the most simple set-up for a call to AutoMixSampler.
Let's analyze it by parts.

The first line includes the AutoMixSampler header file, where the AutoMixSampler class is defined::

    #include "AutoMixSampler.h"

Next we have the declaration of two functions provided by the user::

    double logposterior(int model, int mdim, double* theta);
    void get_rwm_init(int model, int mdim, double *theta);

The first one is the **logarithm** of the `probability density function <https://en.wikipedia.org/wiki/Probability_density_function>`_ (PDF),
(or the log-posterior distribution in the context of a Bayesian analysis) for each of the models in our problem.

The second function simply gives initial values (could be random) in the parameter space of the model indicated by the `model` index.
These values will be used to initialize the Random Walk Metropolis (RWM) part of the program.
For our problem we can set:

.. code-block:: c

  void get_rwm_init(int model, int mdim, double *x) {
    if (model == 1) {
      x[0] = 1.0; // sigma
      x[1] = 0.5; // x0
      x[2] = 0.0; // a
      x[3] = 1.0; // b
      return;
    } else if (model == 2) {
      x[0] = 1.0; // sigma
      x[1] = 0.5; // x0
      x[2] = 1.0; // b
      return;
    } else if (model == 3) {
      x[0] = 2.0; // alpha
      x[1] = 2.0; // beta
      return;
    }
  }

The first line in the main function is the construction of an AutoMixSampler object::

      AutoMixSampler am;

This is the default constructor with all the default vales. See :ref:`constructor` for more detail on the available constructors.

The next two lines, we let `am` know of the log-posterior and initial values for the RWM as explained above::

      am.setTargetDistribution(logposterior);
      am.setRWMInitFunction(get_rwm_init);

The next few lines we define a problem where we have two models, the first one with a parameter space of dimension 1 (one single free parameter), and the other model with two free parameters::

      am.setNumberOfModels(2);
      int model_dims[] = {1, 2};
      am.setModelDimensions(model_dims);

The following line creates the 1,000 samples::

      am.rjmcmc_samples(1000);

Run Statistics
--------------

Our previous main program lacked one fundamental problems: the lack of output.

AutoMix saves the run statistics in different C data structures.

The most important is `rjmcmcStats`.
To ask AutoMix to save `rjmcmcStats` to file, use the helper logwrite functions::

.. code-block:: c

  #include "logwrite.h"
  ...

  int main() {
  ...

  am.rjmcmc_samples(1000);
  report_rjmcmc_run("output", am.st);
  return 0;
  }


Running Conditional Probability Estimation
------------------------------------------

Setting-up Optional Arguments
-----------------------------

Optional Burning of Samples

.. _constructor:

Constructor
-----------

Default.

Reading from file

Setting model and blah blah.
