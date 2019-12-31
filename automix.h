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

#include <stdio.h>
#include "user.h"
#include "utils.h"
// Global constants (please feel free to change as required)

// NMODELS_MAX = maximum number of models
#define NMODELS_MAX 15
// Lkmaxmax = initial number of mixture components fitted in stage 2 of
// AutoMix algorithm
#define Lkmaxmax 30

int read_mixture_params(char *fname, int kmax, int *model_dims, double **sig,
                        int *Lk, double **lambda, double ***mu, double ****B);

void rwn_within_model(int k1, int *model_dims, int nsweep2, FILE *fpl,
                      FILE *fpcf, FILE *fpad, double **sig, int dof,
                      double **data);

void fit_mixture_from_samples(int mdim, double **data, int lendata,
                              double **mu_k, double ***B_k, double *lambda_k,
                              FILE *fpcf, int *Lk_k);

void fit_autorj(double *lambda_k, int *Lk_k, int mdim, double **mu_k,
                double ***B_k, double **data, int lendata);

void reversible_jump_move(int is_burning, double ****B, int *Lk, int adapt,
                          int *current_Lkk, int *current_model_k, double **detB,
                          int do_block_RWM, int dof, int doperm,
                          double **lambda, double *llh, double *lp, int *mdim,
                          int *model_dims, double ***mu, int *naccrwmb,
                          int *naccrwms, int *nacctd, int nmodels, int *nreinit,
                          int *ntryrwmb, int *ntryrwms, int *ntrytd, double *pk,
                          int *pkllim, double *propk, int *reinit, double **sig,
                          double gamma_sweep, double *theta, int mdim_max,
                          int Lkmax);
