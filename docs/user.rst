.. _user:

Writing User Functions
----------------------

For each example to which the user wishes to apply the AutoMix sampler the user must supply a file (which is linked at compile time, see Makefile for examples) containing four user functions written in C.

.. warning::

    A familiarity with the C programming language is assumed.

These functions must have names and arguments detailed below and return the following information:

1. A function that returns the number of models:

  .. code-block:: c

      int get_nmodels(void)

2. A function to load the dimension of each model:

  .. code-block:: c

      void load_model_dims(int nmodels, int *model_dims)

  Given the number of models `nmodels`, on exit from the function the vector `model_dims` should contain the dimensions of the models.

3. A function to get initial conditions (possibly random) for the RWM for a given model:

  .. code-block:: c

      void get_rwm_init(int model_k, int mdim, double *rwm)

  Given the model index k, and the dimension `model_k` of that model, on exit from the function the vector `rwm` should contain the initial conditions.
  In particular, rwm[j-1] should contain the (possibly random) intial state for component j of the parameter vector associated with model k, (j=1,...,mdim).
    
  If random initial states are used, the user must also declare the random number functions in the file (see the example files).

4. A function to return the log of the target function pi evaluated at a given point in the state space, up to an additive constant:

  .. code-block:: c

      void logpost(int model_k, int mdim, double *theta, double* lp, double *llh)

  Given the model index `model_k`, and parameter vector `theta` (of dimension `mdim`), the function must load in `lp` the log of the target function (up to an additive constant) evaluated at this point.
  If pi is a posterior distribution, the double `llh` should contain the log-likelihood evaluated at this point (although this is only necessary for returning the log-likelihood to output file, and can contain any other value if preferred).

The examples provided, with comments, show typical examples of these user files for the problems under consideration. 
