.. _api:

Public API
==========

Functions:
----------

All public functions require as a first argument a pointer to :c:type:`amSampler` structure.

.. c:function:: int initAMSampler(amSampler *am, int nmodels, int *model_dims,
                targetFunc logpost, double *initRWM);

    To initialize :c:type:`amSampler` you need to provide:

     * the number of models in `nmodels`;
     * the dimensions for each model in `model_dims`;
     * the logarithm of the target distribution (or posterior distribution in a Bayesian analysis) in logpost (see :c:type:`targetFunc`);
     * and an array with initial conditions for each model in `initRWM`.

    `initRWM` is a flattened array, with the initial values in contiguous order for each model.
    If `initRWM` is passed as `NULL`, then the initial values are randomly selected from a uniform distribution in [0, 1).

    initAMSampler allocates necessary memory for most array in the rest of the structures.

.. c:function:: void freeAMSampler(amSampler *am);

    Once you finish working with :c:type:`amSampler`, you should free the memory with a call to `freeAMSampler`.

.. c:function:: void estimate_conditional_probs(amSampler *am, int nsweep2);

    Before running RJMCMC sweeps, you need to call :c:func:`estimate_conditional_probs` to
    estimate the conditional probabilities of your models, to adjust for an appropriate proposal distribution.

    The statistics collected during the :c:func:`estimate_conditional_probs` calls, are stored in the :c:member:`cpstats` member inside the :c:type:`amSampler` struct.
    For a full description of :c:member:`cpstats`, see :c:type:`condProbStats`.

.. c:function:: void burn_samples(amSampler *am, int nburn);

    It is advised to burn some samples at the begining of a RJMCMC run to let the Markov chain converge.

.. c:function:: void rjmcmc_samples(amSampler *am, int nsweep);

    This is the function that will provide the RJMCMC samples.

    The statistics collected during the :c:func:`rjmcmc_samples` calls, are stored in the :c:member:`st` member inside the :c:type:`amSampler` struct.
    For a full description of :c:member:`st`, see :c:type:`runStats`.

C struct's
----------

.. c:type:: runStats

  Contains statistics of the `rjmcmc_samples` call.

  .. c:member:: naccrwmb

    [unsigned long] Number of accepted Block-RWM.

  .. c:member:: ntryrwmb

    [unsigned long] Number of tried Block-RWM.

  .. c:member:: naccrwms

    [unsigned long] Number of accepted Single-RWM.

  .. c:member:: ntryrwms

    [unsigned long] Number of tried Single-RWM.

  .. c:member:: nacctd

    [unsigned long] Number of accepted Auto RJ.

  .. c:member:: ntrytd

    [unsigned long] Number of tried Auto RJ.

  .. c:member:: theta_summary

    [`double ***`] The theta parameter for each model, for each sweep.
    theta_summary[i][j][k] is the :math:`\theta_k` component for the j-th sweep for the i-th model.

  .. c:member:: theta_summary_len

    [`int *`] the number of sweeps for theta_summary in model k.
    theta_summary_len[1] = 20 means that the model 1 has 20 :math:`\theta` samples.

  .. c:member:: ksummary

    [`int *`] the number of times the model was visited.


.. c:type:: amSampler

  .. c:member:: ch

    A :c:type:`chainState` instance containing the chain state of the RJMCMC.

  .. c:member:: jd

    A :c:type:`proposalDist` instance containing the proposal distribution for the chain.

  .. c:member:: cpstats

    A :c:type:`condProbStats` instance that holds statistics for the estimation of conditional probabilities (if ran).

  .. c:member:: st

    A :c:type:`runStats` instance that holds statistics for the RJMCMC runs (if ran).

  .. c:member:: doAdapt

    A :c:type:`bool` value indicating if Adaptation is to be performed or not during the RJMCMC runs.
    Default value is `true`.

  .. c:member:: doPerm

    A :c:type:`bool` value indicating if Permutation is to be performed or not during the RJMCMC runs.
    Default value is `false`.

  .. c:member:: student_T_dof

    int type. The degree of freedom of the Student's T distribution to sample from.
    If 0, sample from a Normal distribution instead.
    Default value is 0.

  .. c:member:: am_mixfit

    A :c:type:`automix_mix_fit` value. Default is `FIGUEREIDO_MIX_FIT`.

  .. c:member:: seed

    An unsigned long value to initialize the random number generator. Defaults to clock time.

.. c:type:: condProbStats

  .. c:member:: rwm_summary_len

    TBD.

  .. c:member:: sig_k_rwm_summary

    TBD.

  .. c:member:: nacc_ntry_rwm

    TBD.

  .. c:member:: nfitmix

    TBD.

  .. c:member:: fitmix_annulations

    TBD.

  .. c:member:: fitmix_costfnnew

    TBD.

  .. c:member:: fitmix_lpn

    TBD.

  .. c:member:: fitmix_Lkk

    TBD.

.. c:type:: proposalDist

  The proposal distribution is of the form:

  .. math::

    p(\theta) = \sum_{k=1}^{k_{max}} L_k \exp \left(-\frac{1}{2} \left( \theta_k - \mu_k \right) B \cdot B^T \left( \theta_k - \mu_k \right) \right).

  .. c:member:: nmodels

    The number of models in M (referred as :math:`K_{max}` in thesis, p 143)

  .. c:member:: nMixComps

    The numberof mixture components for each model (refferred as :math:`L_k` in thesis, p144)

  .. c:member:: model_dims

    The dimension of the :math:`\theta` parameter space for each model.

  .. c:member:: lambda

    The relative weights for each mixture component for each model.
    lambda[i][j] is the weight of the j-th mixture component for model :math:`M_i`.

  .. c:member:: mu

    The mean vector for the mixture component for each model.
    mu[i][j][k] is the k-th vector index for the j-th mixture component for model :math:`M_i`.

  .. c:member:: B

    :math:`B \cdot B^T` is the covariance matrix for each mixture component for each model.
    B[i][j][k][m] is the k,m index of the B matrix for mixture component j of model :math:`M_i`.

  .. c:member:: sig

    Vector of adapted RWM scale parameters for each model.

.. c:type:: chainState

  Struct to hold the MCMC chain state

  .. c:member:: theta

    The current value of :math:`\theta_k` in the chain.

  .. c:member:: pk

    TBD.

  .. c:member:: log_posterior

    The current value of the log-posterior.

  .. c:member:: current_model_k

    The current model `k` in the chain.

  .. c:member:: mdim

    The dimension of the current model in the chain.

  .. c:member:: current_Lkk

    TBD

  .. c:member:: nreinit

    TBD.

  .. c:member:: reinit

    TBD.

  .. c:member:: pkllim

    TBD.

  .. c:member:: doBlockRWM

    A :c:type:`bool` value to indicate if the chain should do a Block-RWM.
    This is done every 10 sweeps.

  .. c:member:: isBurning

    A :c:type:`bool` that indicates wether the chain is burning samples.

  .. c:member:: sweep_i

    The index of the chain.


Typedef's and Enum's
--------------------

.. c:type:: targetFunc

    targetFunc is the prototype of the log-posterior function that must be passed to :c:type:`amSampler` upon initialization.

    .. code-block:: c

      typedef double (*targetFunc)(int model_k, double *x);

.. c:type:: bool

    bool is a typedef for int. C does not have a bool type but int is just as good.

    .. code-block:: c

      typedef int bool;

.. c:type:: automix_mix_fit

    Enum to specify whether using Figuereido or AutoRJ in conditional probabilities estimation.

    .. code-block:: c

      typedef enum { FIGUEREIDO_MIX_FIT = 0, AUTORJ_MIX_FIT } automix_mix_fit;
