.. AutoMix documentation master file, created by
   sphinx-quickstart on Fri Jan  3 09:38:43 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index:

AutoMix
=======

.. toctree::
   :maxdepth: 2
   :hidden:

   tutorial
   compilation
   api
   license

The AutoMix package is a C program for Unix-like systems, implementing the automatic Reversible Jump MCMC sampler of the same name described in Chapters 4, 5, and 6 of David Hastie's Ph.D. thesis (included in ``docs/thesis``).

While the original AutoMix is highly useful, the fact that it can only be used as an executable can limit its applicability.
**LibAutoMix makes the core algorithms of Source Extractor available as a library of stand-alone functions and data structures.**
The code is completely derived from the original AutoMix code base and aims to produce results compatible with it whenever possible.
LibAutoMix is a C library with no dependencies outside the standard library.

What is Reversible Jump MCMC?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reversible Jump Markov Chain Monte Carlo (RJMCMC) extends the standard MCMC sampler to include a
discrete random variable ``k`` that represents the index of a *model*.
So, instead of sampling from the usual parameter space of a given distribution,
RJMCMC will also sample across different models (distributions).

The final samples of a RJMCMC reflect the probabilities of the parameters of each model,
but also the relative probabilities of the models themselves.

Main advantages of AutoMix:
^^^^^^^^^^^^^^^^^^^^^^^^^^^

  * Reversible Jump MCMC allows sampling from several distributions (models) simultaneously.
  * The different models may have different number of free parameters (dimension).
  * The relative frequency of sampling for different models is proportional to the probability of the model.
    That means that AutoMix can be used as a model selection sampler.
  * AutoMix requires minimum input from the user.
  * AutoMix automatically adapts proposal distributions with a multi-modal Normal mixture.

.. caution::

    Potential users should carefully understand the limitations of using the AutoMix sampler.
    The reliability of results from this sampler depends on many factors, including the scaling of the parameters and the degree of multimodality of the within-model conditionals of the target distribution.

Where to go from here?
^^^^^^^^^^^^^^^^^^^^^^

For a quick tour of the library, check out the :ref:`tutorial`.

To compile your program against the AutoMix library see :ref:`compile`.

For a description of the full public API with statistics information, see :ref:`api`.
