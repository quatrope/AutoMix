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

#define AUTOMIX_MAJOR_VERSION 1
#define AUTOMIX_MINOR_VERSION 3
#define AUTOMIX_REVISION 0
#define AUTOMIX_VERSION "1.3"

#include <time.h>

// Global constants (please feel free to change as required)
// NMODELS_MAX = maximum number of models
#define NMODELS_MAX 15 // used in initProposalDist
// Lkmaxmax = initial number of mixture components fitted in stage 2 of
// AutoMix algorithm
#define NUM_MIX_COMPS_MAX                                                      \
  30 // used in initJD, freeJD and fit_mixture_from_samples
#define NUM_FITMIX_MAX                                                         \
  5000 // used in initCondProbStats and fit_mixture_from_samples

typedef double (*targetFunc)(int model_k, int mdim, double *x);
typedef void (*rwmInitFunc)(int model_k, int mdim, double *x);
// C does not have a bool type but int is just as good
typedef int bool;
typedef enum { FIGUEREIDO_MIX_FIT = 0, AUTORJ_MIX_FIT } AUTOMIX_MIX_FIT;

// Struct to hold the MCMC chain state
typedef struct {
  double *theta;
  double *pk;
  double log_posterior;
  int current_model_k;
  int mdim;
  int current_Lkk;
  int nreinit;
  int reinit;
  double pkllim;
  bool doBlockRWM;
  bool isBurning;
  unsigned long sweep_i;
  bool isInitialized;
} chainState;

typedef struct {
  // The number of models in M (referred a K_{max} in thesis, p 143)
  int nmodels;
  // The numberof mixture components for each model (refferred as L_k in thesis,
  // p144)
  int *nMixComps;
  // The dimension of the theta parameter space for each model.
  int *model_dims;
  // The relative weights for each mixture component for each model
  // lambda[i][j] is the weight of the j-th mixture component for model M_i.
  double **lambda;
  // The mean vector for the mixture component for each model
  // mu[i][j][k] is the k-th vector index for the j-th mixture component for
  // model M_i.
  double ***mu;
  // B.B^T is the covariance matrix for each mixture component for each model.
  // B[i][j][k][m] is the k,m index of the B matrix for mixture component j of
  // model M_i.
  double ****B;
  // vector of adapted RWM scale parameters for each model
  double **sig;
  bool isInitialized;
  bool isAllocated;
} proposalDist;

typedef struct {
  int rwm_summary_len;
  double ***sig_k_rwm_summary;
  double ***nacc_ntry_rwm;
  int *nfitmix;
  int **fitmix_annulations;
  double **fitmix_costfnnew;
  double **fitmix_lpn;
  int **fitmix_Lkk;
  // Booleans to specify whether using Figuereido or AutoRJ in conditional
  // probs estimation.
  bool useAutoMixFit;
  bool useAutoRJFit;
  double timesecs_condprobs;
  bool isInitialized;
} condProbStats;

typedef struct {
  // Block RWM acceptance and tries
  unsigned long naccrwmb;
  unsigned long ntryrwmb;
  // Single RWM acceptance and tries
  unsigned long naccrwms;
  unsigned long ntryrwms;
  // Auto RJ acceptance and tries
  unsigned long nacctd;
  unsigned long ntrytd;
  double ***theta_summary;
  int *theta_summary_len;
  int *theta_summary_size;
  // nsokal is every how many samples are kept in array xr
  int nsokal;
  // nkeep is how many samples are kept in array xr
  int nkeep;
  // keep is the index from which to start keeping samples in xr
  // (kept every nsokal samples)
  int keep;
  int m;
  double *xr;
  double var;
  double tau;
  int *ksummary;
  double **pk_summary;
  int *k_which_summary;
  double **logp_summary;
  double timesecs_rjmcmc;
  double timesecs_burn;
} runStats;

typedef struct {
  // int NMODELS_MAX;
  // int NUM_MIX_COMPS_MAX;
  // int NUM_FITMIX_MAX;
  chainState ch;
  proposalDist jd;
  condProbStats cpstats;
  runStats st;
  bool doAdapt;
  bool doPerm;
  targetFunc logposterior;
  rwmInitFunc initRWM;
  int student_T_dof;
  AUTOMIX_MIX_FIT am_mixfit;
  unsigned long seed;
} amSampler;

/*** Constructors and Destructors ***/
int initAMSampler(amSampler *am, int nmodels, int *model_dims,
                  targetFunc logpost, rwmInitFunc initRWM);
void freeAMSampler(amSampler *am);
void initRunStats(runStats *st, int nsweep, proposalDist jd);
void freeRunStats(runStats st, proposalDist jd);

/*** Public Functions ***/
int read_mixture_params(char *fname, amSampler *am);
void estimate_conditional_probs(amSampler *am, int nsweep2);
void burn_samples(amSampler *am, int nburn);
void rjmcmc_samples(amSampler *am, int nsweep);

#endif /* defined(__AutoMix__AutoMix_h__) */
