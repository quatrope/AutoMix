.. AutoMix documentation master file, created by
   sphinx-quickstart on Fri Jan  3 09:38:43 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index:

AutoMix
=======

The AutoMix package is a C program for Unix-like systems, implementing the automatic reversible jump MCMC sampler of the same name described in Chapters 4, 5, and 6 of David Hastie's Ph.D. thesis included in `docs/thesis`.

.. note::

    Potential users should carefully understand the limitations of using the AutoMix sampler.
    The reliability of results from this sampler depends on many factors, including the scaling of the parameters and the degree of multimodality of the within-model conditionals of the target distribution.

Usage
-----

To run the sampler for a particular problem the program must be compiled with a user-provided file containing four C functions that define the problem in question.
Examples of such files for problems considered within the aforementioned thesis are bundled with the archived AutoMix software.
More details of the functions are provided below.

The output of the sampler is in the form of several text files summarizing the MCMC sampler.
We discuss these files below.
Subsequent analysis may be performed on the output with the use of statistical packages.
We recommend `R <http://www.stats.bris.ac.uk/R/>`_ as a good free choice of such a package.

Contents:
^^^^^^^^^

.. toctree::
   :maxdepth: 2

   compilation
   running
   user
   output
   license
