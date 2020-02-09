/* The AutoMix program.

Last edited 25/11/04.
Developed by David Hastie, Department of Mathematics,
University of Bristol, UK as a part of a submission for the
degree of Ph.D. This Ph.D. was supervised by Prof. Peter Green (PJG),
University of Bristol. Special thanks also to Dr. Christophe Andrieu CA),
University of Bristol, for advice on adaptive schemes and mixture fitting.

The AutoMix sampler is free for personal and academic use, but must
reference the sampler as instructed below.  For commercial
use please permission must be sought from the author. To seek permission
for such use please send an e-mail to d_hastie@hotmail.com
outlining the desired usage.

Use of the AutoMix sampler is entirely at the user's own risk. It is the
responsibility of the user to ensure that any conclusions made through the
use of the AutoMix sampler are valid. The author accepts no responsibility
whatsoever for any loss, financial or otherwise, that may arise in
connection with the use of the sampler.

The AutoMix sampler is available from http://www.davidhastie.me.uk/AutoMix
Although the sampler may be modified and redistributed, the author
encourages users to register at the above site so that updates of the
software can be received.

Before use, please read the README file bundled with this software.

Users should reference the sampler as instructed on the AutoMix website
(see above). Initially this is likely to be the Ph.D. thesis that
introduces the AutoMix sampler. However, this will hopefully change to
be a published paper in the not too distant future.  */

#ifndef __AutoMix__AutoMix_h__
#define __AutoMix__AutoMix_h__

#define AUTOMIX_MAJOR_VERSION 2
#define AUTOMIX_MINOR_VERSION 0
#define AUTOMIX_REVISION 0
#define AUTOMIX_VERSION "2.0"

#include <time.h>

// The prototype of the log-posterior function that must be passed to amSampler
// upon initialization.
typedef double (*targetDist)(int model_k, double *x);

// C does not have a bool type but int is just as good.
typedef int bool;
// true and false are also defined for compatibility with C++.
#define true 1
#define false 0

// Enum to specify whether using Figuereido or AutoRJ in conditional
// probs estimation.
typedef enum { FIGUEREIDO_MIX_FIT = 0, AUTORJ_MIX_FIT } automix_mix_fit;

/* -------------------------------------------------------------------------

Struct Typedef's

------------------------------------------------------------------------- */

// Holds information for the AutoMix sampler
typedef struct amSampler amSampler;

// Holds the MCMC chain state
typedef struct chainState chainState;

// Proposal distribution parameters
typedef struct proposalDist proposalDist;

// Statistics gathered during the conditional probability estimation
typedef struct condProbStats condProbStats;

// Statistics gathered during the RJMCMC samples
typedef struct runStats runStats;

/* -------------------------------------------------------------------------

Public Functions

------------------------------------------------------------------------- */

// The initializer for AutoMix Sampler.
int initAMSampler(amSampler *am, int nmodels, int *model_dims,
                  targetDist logpost, double *initRWM);
// Memory free for AutoMix Sampler
void freeAMSampler(amSampler *am);

// Estimate the conditional probabilities of the models,
// to adjust for an appropriate proposal distribution.
void estimate_conditional_probs(amSampler *am, int nsweeps);

// Burns nsweeps samples at the begining of a RJMCMC run
// to let the Markov chain converge.
void burn_samples(amSampler *am, int nsweeps);

// The function that provides the RJMCMC samples.
void rjmcmc_samples(amSampler *am, int nsweeps);

/* -------------------------------------------------------------------------

Struct definitions

------------------------------------------------------------------------- */

/*
 * chainState
 *
 * Holds the MCMC chain state.
 */
struct chainState {
  double *theta;         /* The theta parameters of the model */
  double *pk;            /* TBD */
  double log_posterior;  /* The log_posterior for (k, theta) */
  int current_model_k;   /* The model index, k */
  int mdim;              /* The model dimension (dim theta) */
  int current_Lkk;       /* The */
  int nreinit;           /* Number of re-init's */
  int reinit;            /* Re-init */
  double pkllim;         /* TBD */
  bool doBlockRWM;       /* Flag to do a Block Random Walk Metropolis */
  bool isBurning;        /* Flag for wether the sweeps are burnt */
  unsigned long sweep_i; /* The sweep index */
  bool isInitialized;    /* Flag to know if the arrays have been allocated */
};

/*
 * proposalDist
 *
 * Proposal distribution parameters.
 */
struct proposalDist {
  int nmodels;     /* Number of models in M (referred a K_{max} in thesis,
                      p 143). */
  int *nMixComps;  /* Number of mixture components for each model (refferred as
                      L_k in thesis, p144). */
  int *model_dims; /* The dimension of the theta parameter space for each model.
                    */
  double **lambda; /* The relative weights for each mixture component for each
                      model. lambda[i][j] is the weight of the j-th mixture
                      component for model M_i. */
  double ***mu;    /* The mean vector for the mixture component for each model.
                      mu[i][j][k] is the k-th vector index for the j-th mixture
                      component for model M_i. */
  double ****B; /* B.B^T is the covariance matrix for each mixture component for
                   each model. B[i][j][k][m] is the k,m index of the B matrix
                   for mixture component j of model M_i. */
  double **sig; /* vector of adapted RWM scale parameters for each model */
  int NUM_MIX_COMPS_MAX; /* Maximum number of components in the mixture. */
  bool isInitialized;    /* Flag to know if the arrays have been allocated. */
};

/*
 * condProbStats
 *
 * Statistics gathered during the conditional probability estimation.
 */
struct condProbStats {
  int rwm_summary_len;         /* The length of sig_k_rwm_summary */
  double ***sig_k_rwm_summary; /* TBD */
  double ***nacc_ntry_rwm;     /* Ratio of accepted to tried RWM */
  int *nfitmix;                /* Fitmix parameter */
  int **fitmix_annulations;    /* Fitmix parameter */
  double **fitmix_costfnnew;   /* Fitmix parameter */
  double **fitmix_lpn;         /* Fitmix parameter */
  int **fitmix_Lkk;            /* Fitmix parameter */
  double timesecs_condprobs;   /* Time spent in estimate_conditional_probs */
  bool isInitialized; /* Flag to know if conditional probilities were estimated.
                       */
};

/*
 * runStats
 *
 * Statistics gathered during the RJMCMC samples.
 */
struct runStats {
  unsigned long naccrwmb;  /* Block RWM acceptances */
  unsigned long ntryrwmb;  /* Block RWM tries */
  unsigned long naccrwms;  /* Single RWM acceptances */
  unsigned long ntryrwms;  /* Single RWM tries */
  unsigned long nacctd;    /* Auto RJ acceptances */
  unsigned long ntrytd;    /* Auto RJ tries */
  double ***theta_summary; /* The theta values drawn. theta_summary[i][j][k]
                              corresponds to the kth component of the jth sweep
                              for model i. */
  int *theta_summary_len;  /* The length of theta_summary[i] (number of sweeps
                              for model i) */
  int *theta_summary_size; /* Size allocated for theta_summary[i]. Internal. */
  int nsokal;    /* nsokal is every how many samples are kept in array xr */
  int nkeep;     /* nkeep is how many samples are kept in array xr */
  int keep;      /* keep is the index from which to start keeping samples in xr,
                    (kept every nsokal samples) */
  int m;         /* TBD */
  double *xr;    /* TBD */
  double var;    /* TBD */
  double tau;    /* TBD */
  int *ksummary; /* The number of calls to model k */
  double **pk_summary;    /* TBD */
  int *k_which_summary;   /* The sequence of model index calls */
  double **logp_summary;  /* the logposterior values for each model. */
  double timesecs_rjmcmc; /* Time spent in rjmcmc_samples() */
  double timesecs_burn;   /* Time spent in burn_samples() */
  bool isInitialized;     /* Flag to determine if arrays were allocated */
};

/*
 * amSampler
 *
 * Holds information for the AutoMix sampler.
 */
struct amSampler {
  int NMODELS_MAX;       /* Maximum number of models */
  int NUM_MIX_COMPS_MAX; /* Max number of components in proposal distribution */
  int NUM_FITMIX_MAX;    /* TBD */
  chainState ch;         /* Internal instance of chainState struct */
  proposalDist jd;       /* Internal instance of proposalDist struct */
  condProbStats cpstats; /* Internal instance of condProbStats struct */
  runStats st;           /* Internal instance of rnStats struct */
  bool doAdapt;          /* Flag to perform adaptation */
  bool doPerm;           /* Flag to perform permutation */
  targetDist logposterior; /* Pointer to the target distribution function */
  double **initRWM;        /* Initial values to call model log-posteriors */
  int student_T_dof;       /* Degree of freedom for Student T's Distribution */
  automix_mix_fit am_mixfit; /* Enum to choose Cond Prob Estimation */
  unsigned long seed;        /* Seed for random number generator */
};

#endif /* defined(__AutoMix__AutoMix_h__) */
