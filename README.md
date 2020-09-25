[![Build Status](https://travis-ci.org/quatrope/AutoMix.svg?branch=master)](https://travis-ci.org/quatrope/AutoMix)
[![Documentation Status](https://readthedocs.org/projects/automix/badge/?version=latest)](https://automix.readthedocs.io/en/latest/?badge=latest)
![Style Format](https://img.shields.io/badge/style-clang--format-success)
[![CodeFactor](https://www.codefactor.io/repository/github/quatrope/automix/badge/master)](https://www.codefactor.io/repository/github/quatrope/automix/overview/master)

# AutoMix

## The AutoMix Sampler
 
_This README is a modified version of the original README.txt file._
_Check the original README.txt for license information._

The AutoMix package is a C program for Unix-like systems, implementing the automatic Reversible Jump MCMC sampler of the same name described in Chapters 4, 5, and 6 of David Hastie's Ph.D. thesis (included in `docs/thesis`).

While the original AutoMix is highly useful, the fact that it can only be used as an executable can limit its applicability.
**LibAutoMix makes the core algorithms of Source Extractor available as a library of stand-alone functions and data structures.**
The code is completely derived from the original AutoMix code base and aims to produce results compatible with it whenever possible.
LibAutoMix is a C library with no dependencies outside the standard library.

### Warning:

> Potential users should carefully understand the limitations of using the AutoMix sampler.
> The reliability of results from this sampler depends on many factors,
> including the scaling of the parameters and the degree of multimodality of the within-model conditionals of the target distribution.

## What is Reversible Jump MCMC?

Reversible Jump Markov Chain Monte Carlo (RJMCMC) extends the standard MCMC sampler to include a
discrete random variable `k` that represents the index of a *model*.
So, instead of sampling from the usual parameter space of a given distribution,
RJMCMC will also sample across different models (distributions).

The final samples of a RJMCMC reflect the probabilities of the parameters of each model,
but also the relative probabilities of the models themselves.

## Main advantages of AutoMix:

* Reversible Jump MCMC allows sampling from several distributions (models) simultaneously.
* The different models may have different number of free parameters (dimension).
* The relative frequency of sampling for different models is proportional to the probability of the model.
  That means that AutoMix can be used as a model selection sampler.
* AutoMix requires minimum input from the user.
* AutoMix automatically adapts proposal distributions with a multi-modal Normal mixture.

## Simple example program

AutoMix is most useful in the situation where multiple models are tested simultaneously,
Nevertheless, below is an example for a simple 1D Normal distribution sampler.

    // main.c file
    #include "automix.h"
    double logp_normal_sampler(int model, double* xp);
    int main() {
      int nmodels = 1;
      int model_dims[] = {1};
      double initRWM[] = {0.5, 0.5};
      amSampler am;
      initAMSampler(&am, nmodels, model_dims, logp_normal_sampler, initRWM);
      burn_samples(&am, 10000);
      int nsweeps = 100000;
      rjmcmc_samples(&am, nsweeps);
      freeAMSampler(&am);
      return 0;
    }
    
    double logp_normal_sampler(int model_k, double *xp) {
      double x = *xp;
      double prob;
      double x0 = 0.5;
      double sigma = 1.0;
      prob = -(x - x0) * (x - x0) / (2.0 * sigma * sigma);
      return prob;
    }

The previous C file will output samples from a Normal distribution.
To compile:

    $ make; make install
    $ sudo ldconfig # only for Linux
    $ cc main.c -lautomix -o normal
    $ ./normal

For full documentation, examples and a quick tour of the library, check out the [documentation](https://automix.readthedocs.io/).
