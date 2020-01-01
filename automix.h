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

#include "user.h"
#include "utils.h"
#include <stdio.h>

// Global constants (please feel free to change as required)
// NMODELS_MAX = maximum number of models
#define NMODELS_MAX 15
// Lkmaxmax = initial number of mixture components fitted in stage 2 of
// AutoMix algorithm
#define NUM_MIX_COMPS_MAX 30

// C does not have a bool type but int is just as good
typedef int bool;
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
  int nmodels;
  int *nMixComps;
  int *model_dims;
  double **lambda;
  double ***mu;
  double ****B;
  bool isInitialized;
  bool isAllocated;
} proposalDist;

void initChain(chainState *ch, proposalDist jd, int adapt);
void freeChain(chainState *aChain);
int initJD(proposalDist *jd);
int allocJD(proposalDist *jd);
void freeJD(proposalDist jd);

int read_mixture_params(char *fname, proposalDist jd, double **sig);

void rwm_within_model(int k1, int *model_dims, int nsweep2, FILE *fpl,
                      FILE *fpcf, FILE *fpad, double **sig, int dof,
                      double **data);

void fit_mixture_from_samples(int model_k, proposalDist jd, double **data,
                              int lendata, FILE *fpcf);

void fit_autorj(int model_k, proposalDist jd, double **data, int lendata);

void reversible_jump_move(chainState *ch, proposalDist jd, int dof,
                          int *naccrwmb, int *naccrwms, int *nacctd,
                          int *ntryrwmb, int *ntryrwms, int *ntrytd,
                          double **sig);
