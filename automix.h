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

#define AUTOMIX_MAJOR_VERSION 1
#define AUTOMIX_MINOR_VERSION 3
#define AUTOMIX_REVISION 0
#define AUTOMIX_VERSION "1.3"
#define AUTOMIX_VERSION_CHECK(maj, min)                                        \
  ((maj == AUTOMIX_MAJOR_VERSION) && (min <= AUTOMIX_MINOR_VERSION))

// Global constants (please feel free to change as required)
// NMODELS_MAX = maximum number of models
#define NMODELS_MAX 15
// Lkmaxmax = initial number of mixture components fitted in stage 2 of
// AutoMix algorithm
#define NUM_MIX_COMPS_MAX 30
#define NUM_FITMIX_MAX 5000

#include <time.h>

// C does not have a bool type but int is just as good
typedef int bool;
#ifndef AUTOMIX_DATA_STRUCTS
#define AUTOMIX_DATA_STRUCTS
// Struct to hold the MCMC chain state
typedef struct {
  bool isInitialized;
  double *theta;
  double *pk;
  double log_posterior;
  double log_likelihood;
  int current_model_k;
  int mdim;
  int current_Lkk;
  int nreinit;
  int reinit;
  double pkllim;
  bool doBlockRWM;
  bool doAdapt;
  bool doPerm;
  bool isBurning;
  double gamma_sweep;
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
  bool isInitialized;
  bool isAllocated;
} proposalDist;

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
  int rwm_summary_len;
  double ***sig_k_rwm_summary;
  double ***nacc_ntry_rwm;
  int *nfitmix;
  int **fitmix_annulations;
  double **fitmix_costfnnew;
  double **fitmix_lpn;
  int **fitmix_Lkk;
} runStats;
#endif

void initChain(chainState *ch, proposalDist jd, int adapt);
void freeChain(chainState *aChain);
int initJD(proposalDist *jd);
int allocJD(proposalDist *jd);
void freeJD(proposalDist jd);
void initializeRunStats(runStats *st, int nsweep, int nsweep2, int nburn,
                        proposalDist jd);
void freeRunStats(runStats st, proposalDist jd);

int read_mixture_params(char *fname, proposalDist jd, double **sig);

void rwm_within_model(int k1, int *model_dims, int nsweep2, runStats st,
                      double *sig_k, int dof, double **samples);

void fit_mixture_from_samples(int model_k, proposalDist jd, double **samples,
                              int nsamples, runStats *st);

void fit_autorj(int model_k, proposalDist jd, double **samples, int nsamples);

void reversible_jump_move(chainState *ch, proposalDist jd, int dof,
                          runStats *st, double **sig);

void estimate_conditional_probs(proposalDist jd, int dof, int nsweep2,
                                runStats st, double **sig, int mode,
                                char *fname);

void burn_main_samples(chainState *ch, int nburn, proposalDist jd, int dof,
                       runStats *st, double **sig);

void rjmcmc_samples(chainState *ch, int nsweep, int nburn, proposalDist jd,
                    int dof, runStats *st, double **sig, char *fname);

void flush_final_stats(char *fname, chainState ch, double timesecs,
                       unsigned long seed, int mode, int nsweep, int nsweep2,
                       proposalDist jd, double **sig, runStats st);
