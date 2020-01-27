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

Now that we have defined our problem and log-posteriors, we can set up AutoMix to generate samples from the posteriors of each model.


AutoMix Sampler
---------------

Below is an exanple of a minimal call to generate 1,000 samples from a given problem:

.. code-block:: c

  #include "automix.h"
  
  int main() {
    int nmodels = 3;
    int model_dims[] = {2, 2, 2};
    double initRWM[] = {0.5, 0.5, 2.0, 2.0, 9.0, 2.0};
    amSampler am;
    initAMSampler(&am, nmodels, model_dims, logposterior, initRWM);
    estimate_conditional_probs(&am, 100000);
    burn_samples(&am, 10000);
    int nsweeps = 100000;
    rjmcmc_samples(&am, nsweeps);
    freeAMSampler(&am);
    return 0;
  }

This is a simple set-up for an AutoMix program.
Let's analyze it by parts.

The first line includes the AutoMix header file, where the :c:type:`amSampler` structure and automix functions are defined:

.. code-block:: c

  #include "automix.h"

To initiate the :c:type:`amSampler` struct, we need to set 4 things:

  * The number of models we will use (3 in our example):

    .. code-block:: c

      int nmodels = 3;

  * The dimensions of each model:

    .. code-block:: c

      int model_dims[] = {2, 2, 2};

  * A function that returns the **logarithm** of the `probability density function <https://en.wikipedia.org/wiki/Probability_density_function>`_ (PDF),
    (or the log-posterior distribution in the context of a Bayesian analysis) for each of the models in our problem.
    It must have the prototype:

    .. code-block:: c

      double logposterior(int model, double *params);

    And an implementation:

    .. code-block:: c

      double logposterior(int model, double *params) {
        if (model == 1) {
          return logp_normal(params[0], params[1]);
        } else if (model == 2) {
          return logp_beta(params[0], params[1]);
        } else if (model == 3) {
          return logp_gamma(params[0], params[1]);
        }
        return 0.0;
      }

  * The initial value of the parameters for every model. This should be given as a 1d continuous array with an appropriate initial value
    for each model. For example the Beta distribution is bounded to [0, 1] and any initial value has to be in that interval as well.
    This is to avoid starting the chain in a forbidden region of the parameter space or very far from the average values of the parameters:

    .. code-block:: c

      double initRWM[] = {0.5, 0.5, 2.0, 2.0, 9.0, 2.0};

    In our case we start our chain with values :math:`\sigma` = 0.5, :math:`x_0` = 0.5 for the Gaussian model;
    :math:`\alpha` = 2.0, :math:`\beta` = 2.0 for the Beta model;
    and :math:`\alpha` = 9.0, :math:`\beta` = 2.0 for the Gamma model.

Once all of the above is defined for our problem we can init our :c:type:`amSampler`:

.. code-block:: c

  amSampler am;
  initAMSampler(&am, nmodels, model_dims, logposterior, initRWM);

The next step is to estimate the conditional probabilities for our posterior model to
create a Proposal Distribution with a multi-modal Normal mixture:

.. code-block:: c

  estimate_conditional_probs(&am, 100000);

We just need to pass the automix sampler struct and the number of sweeps we want for the estimation.

It is recommended to "burn" some initial samples to let the MCMC chain achieve convergence.
We do this with the following line:

.. code-block:: c

  burn_samples(&am, 10000);

Finally we can create however many RJMCMC samples as we want:

.. code-block:: c

  int nsweeps = 100000;
  rjmcmc_samples(&am, nsweeps);


Run Statistics
--------------

Our previous main program has one fundamental problem: the lack of output.

AutoMix saves the run statistics in different C data structures.

The most important is the struct :c:type:`runStats`.
In :c:type:`amSampler` is the member :c:member:`st`. There, we will find several statistics, including the parameters sampled during the RJMCMC run.
The two most important ones are :c:member:`st.k_summary` and :c:member:`st.theta_summary`.
For a full description of all the parameters, see :ref:`api`.

Just as an example, let's print out the relative probability of each model:

.. code-block:: c

  #include <stdio.h>

  int main() {
    ...
    rjmcmc_samples(&am, nsweeps);

    printf("p(M=1|E) = %lf\n", (double)am.st.ksummary[0] / nsweeps);
    printf("p(M=2|E) = %lf\n", (double)am.st.ksummary[1] / nsweeps);
    printf("p(M=3|E) = %lf\n", (double)am.st.ksummary[2] / nsweeps);

    freeAMSampler(&am);
    return 0;
  }

This should print something like the following on screen::

  p(M=1|E) = 0.792750
  p(M=2|E) = 0.023890
  p(M=3|E) = 0.183360
