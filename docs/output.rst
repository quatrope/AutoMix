.. _output:

Output Files
------------

The following files are returned by the AutoMix sampler (assuming the filestem is "output")

`output_ac.data`
^^^^^^^^^^^^^^^^

  A file containing the estimated autocorrelation coefficient for the chain. The coefficients go up to the order that was used to compute the autocorrelation time using Sokal's method (see Green and Han, 1992)

`output_adapt.data`
^^^^^^^^^^^^^^^^^^^

  A file summarising the AAP adaptation process for each model - adapting the scale parameter of the single component RWM. For a model with nk components, there are `2*nk` columns, with the odd columns showing the evolution of the scale parameters and the even columns showing the evolution of the average RWM  acceptance rate for the corresponding scale parameter.

`output_cf.data`
^^^^^^^^^^^^^^^^

  A file summarising the evolution of the mixture fitting for each model. In particular, for each model the first column records the number of mixture components at the current iteration, the third column summarises the cost function (to be minimised to find the best mixture) and the final column is an indicator function which takes the value 1 if a component is annihilated naturally at that iteration, 2 if a component is removed by forced annihiliation. (The second column is a term that contibutes to the cost function, see Figueiredo and Jain, 2002).

`output_k.data`
^^^^^^^^^^^^^^^

  A file recording the model index at each sweep of the RJ sampler (after burn-in).

`output_log.data`
^^^^^^^^^^^^^^^^^

  A log file, containing important run statistics, including run-time options, run seed, mixture and RWM parameters for each model, and acceptance rates and autocorrelation statistics for the RJ stage 3. Model probabilities are also returned in this file.

`output_lp.data`
^^^^^^^^^^^^^^^^

  A file recording the log posterior (1st column) and log likelihood (2nd column) at each sweep of the RJ sampler (after burn-in).

`output_mix.data`
^^^^^^^^^^^^^^^^^

  A file containing the fitted mixture parameters for the each model. The file is for use by the AutoMix program and is not easily readable by the user. The file contains the number of models and the dimension of each model. Then for each model in turn, the file records the adapted random walk scaling parameters, the number of components of the mixture fitted to that model, and for each component in turn the weight, the mean and the lower triangle of the B matrix. It is this file that is required by the AutoMix sampler if it is run in mode 1.

`output_pk.data`
^^^^^^^^^^^^^^^^

  A file containing the evolution of the model jumping proposal parameters for the RJ stage 3. Column k is the probability of jumping to model k.

and for each model k...

`output_thetak.data`
^^^^^^^^^^^^^^^^^^^^

  A file recording the parameters conditional on model k at each sweep.
