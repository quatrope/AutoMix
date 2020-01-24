// AutoMix - By David Hastie.
// See automix.h for full license and credits.

#include "automix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

/***********************************************************************
 *                                                                      *
 *                          UTILITY FUNCTIONS                           *
 *                                                                      *
 ************************************************************************/

void gauss(double *z, int n);
void rt(double *z, int n, int dof);
void chol(int n, double **B);
void perm(double *work, int n);
double ltprob(int dof, double z);
double lnormprob(int n, double *mu_k_l, double **B_k_l, double *datai);
double det(int n, double **B_k_l);
/* Combined congruential and Tauseworthe generators from SuperDuper
 * package. Should work on machines with unsigned long of at least 32
 * bits. JC and JT must be initialized to values with 0 < JC < 2^32 and
 * 0 < JT < 2^32. JC must be odd.
 * References: Marsaglia, Ananthanarayanan & Paul, 1973,
 * Learmonth & Lewis, 1973, and Dudewicz, 1976)
 */
void sdrni(unsigned long *seed);
double sdrand(void);
// Returns random sample from a Gamma distribution
double rgamma(double s);
// Calculates the logarithm of the Gamma function at x
double loggamma(double x);

/************************************************************************
 *                                                                      *
 *                   PRIVATE FUNCTION PROTOTYPES                        *
 *                                                                      *
 ************************************************************************/

// Constructors and Destructors
int initProposalDist(proposalDist *jd, int nmodels, int *model_dims,
                     int NUM_MIX_COMPS_MAX);
void freeProposalDist(proposalDist jd);
int initCondProbStats(condProbStats *cpstats, proposalDist jd, int nsweeps2,
                      int NUM_FITMIX_MAX);
void freeCondProbStats(condProbStats *cpstats, proposalDist jd);
void initChain(chainState *ch, proposalDist jd, double **initRWM,
               targetFunc logposterior);
void freeChain(chainState *ch);
void initRunStats(runStats *st, int nsweep, proposalDist jd);
void freeRunStats(runStats st, proposalDist jd);

// Helper functions
void rwm_within_model(int k1, int *model_dims, int nsweep2,
                      condProbStats *cpstats, double *sig_k, int dof,
                      double **samples, targetFunc logpost, double *initRWM);
void fit_mixture_from_samples(int model_k, proposalDist jd, double **samples,
                              int nsamples, condProbStats *cpstats,
                              int NUM_MIX_COMPS_MAX, int NUM_FITMIX_MAX);
void fit_autorj(int model_k, proposalDist jd, double **samples, int nsamples);
void reversible_jump_move(bool doPerm, bool doAdapt, chainState *ch,
                          proposalDist jd, int dof, runStats *st,
                          targetFunc logpost);

/***********************************************************************
 *                                                                      *
 *                   PRIVATE FUNCTION IMPLEMENTATION                    *
 *                                                                      *
 ************************************************************************/

void rjmcmc_samples(amSampler *am, int nsweep) {
  clock_t starttime = clock();
  initChain(&(am->ch), am->jd, am->initRWM, am->logposterior);
  initRunStats(&(am->st), nsweep, am->jd);
  chainState *ch = &(am->ch);
  proposalDist jd = am->jd;
  runStats *st = &(am->st);
  // Start here main sample
  int xr_i = 0;
  ch->isBurning = 0;
  for (unsigned long sweep_init = ch->sweep_i;
       ch->sweep_i < sweep_init + nsweep; ch->sweep_i++) {
    // for (int sweep = 1; sweep <= nsweep; sweep++) {
    unsigned long sweep = ch->sweep_i - sweep_init;
    // Every 10 sweeps to block RWM
    ch->doBlockRWM = (ch->sweep_i % 10 == 0);

    reversible_jump_move(am->doPerm, am->doAdapt, ch, am->jd, am->student_T_dof,
                         st, am->logposterior);

    (st->ksummary[ch->current_model_k])++;
    st->k_which_summary[sweep] = ch->current_model_k + 1;
    st->logp_summary[sweep][0] = ch->log_posterior;
    for (int k1 = 0; k1 < jd.nmodels; k1++) {
      st->pk_summary[sweep][k1] = ch->pk[k1];
    }
    int sample_i = st->theta_summary_len[ch->current_model_k];
    int size = st->theta_summary_size[ch->current_model_k];
    if (sample_i == size - 1) {
      int newsize = (size * 2) > nsweep ? nsweep : (2 * size);
      st->theta_summary[ch->current_model_k] = realloc(
          st->theta_summary[ch->current_model_k], newsize * sizeof(double **));
      st->theta_summary_size[ch->current_model_k] = newsize;
    }
    st->theta_summary[ch->current_model_k][sample_i] =
        (double *)malloc(ch->mdim * sizeof(double));
    double *theta_k_i = st->theta_summary[ch->current_model_k][sample_i];
    for (int i = 0; i < ch->mdim; i++) {
      theta_k_i[i] = ch->theta[i];
    }
    (st->theta_summary_len[ch->current_model_k])++;

    if (sweep > st->keep && ((sweep - st->keep) % st->nsokal == 0)) {
      st->xr[xr_i++] = ch->current_model_k;
    }
  }
  for (int model_k = 0; model_k < jd.nmodels; model_k++) {
    int finalsize = st->theta_summary_len[model_k];
    st->theta_summary[model_k] =
        realloc(st->theta_summary[model_k], finalsize * sizeof(double **));
  }
  clock_t endtime = clock();
  st->timesecs_rjmcmc = (endtime - starttime) / (double)CLOCKS_PER_SEC;
}

void burn_samples(amSampler *am, int nburn) {
  clock_t starttime = clock();
  initChain(&(am->ch), am->jd, am->initRWM, am->logposterior);
  chainState *ch = &(am->ch);
  proposalDist jd = am->jd;
  runStats *st = &(am->st);
  ch->isBurning = 1;
  for (unsigned long sweep_init = ch->sweep_i; ch->sweep_i < sweep_init + nburn;
       ch->sweep_i++) {
    // Every 10 sweeps to block RWM
    ch->doBlockRWM = (ch->sweep_i % 10 == 0);

    reversible_jump_move(am->doPerm, am->doAdapt, ch, jd, am->student_T_dof, st,
                         am->logposterior);
  }
  clock_t endtime = clock();
  st->timesecs_burn = (endtime - starttime) / (double)CLOCKS_PER_SEC;
}

void estimate_conditional_probs(amSampler *am, int nsweep2) {
  clock_t starttime = clock();
  proposalDist jd = am->jd;
  condProbStats *cpstats = &(am->cpstats);
  initCondProbStats(cpstats, jd, nsweep2, am->NUM_FITMIX_MAX);
  // Section 5.2 - Within-model runs if mixture parameters unavailable
  for (int model_k = 0; model_k < jd.nmodels; model_k++) {
    int mdim = jd.model_dims[model_k];
    int nsamples = 1000 * mdim;
    // samples holds the (sub-sampled) RWM output theta^1,...,theta^{1000 *
    // mdim} for each model k
    double **samples = (double **)malloc(nsamples * sizeof(*samples));
    samples[0] = (double *)malloc(nsamples * mdim * sizeof(**samples));
    for (int i = 1; i < nsamples; i++) {
      samples[i] = samples[i - 1] + mdim;
    }
    // Section 5.2.1 - Random Walk Metropolis (RWM) Within Model
    // Adapt within-model RWM samplers and to provide the next stage with
    // samples from pi(theta_k|k) for each value of k. (see thesis, p 144)
    rwm_within_model(model_k, jd.model_dims, nsweep2, cpstats, jd.sig[model_k],
                     am->student_T_dof, samples, am->logposterior,
                     (am->initRWM)[model_k]);
    if (am->am_mixfit == FIGUEREIDO_MIX_FIT) {
      // Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)
      // Fit a Normal mixture distribution to the conditional target
      // distributions pi(theta_k|k). See theis, p 144.
      fit_mixture_from_samples(model_k, jd, samples, nsamples, cpstats,
                               am->NUM_MIX_COMPS_MAX, am->NUM_FITMIX_MAX);
    }
    if (am->am_mixfit == AUTORJ_MIX_FIT) {
      //--- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --
      fit_autorj(model_k, jd, samples, nsamples);
    }
    free(samples[0]);
    free(samples);
  }
  clock_t endtime = clock();
  cpstats->timesecs_condprobs = (endtime - starttime) / (double)CLOCKS_PER_SEC;
}

int initAMSampler(amSampler *am, int nmodels, int *model_dims,
                  targetFunc logposterior, double **initRWM) {
  // Proposal Distribution must be initialized immediately
  if (nmodels < 0) {
    printf("Error: negative number of models.\n");
    return EXIT_FAILURE;
  }
  am->NMODELS_MAX = 15;
  am->NUM_MIX_COMPS_MAX = 30;
  am->NUM_FITMIX_MAX = 5000;
  am->seed = 0;
  sdrni(&(am->seed));
  initProposalDist(&(am->jd), nmodels, model_dims, am->NUM_MIX_COMPS_MAX);
  am->logposterior = logposterior;
  am->initRWM = (double **)malloc(nmodels * sizeof(double *));
  for (int i = 0; i < nmodels; i++) {
    (am->initRWM)[i] = (double *)malloc(model_dims[i] * sizeof(double));
  }
  if (initRWM == NULL) {
    for (int m_k = 0; m_k < nmodels; m_k++) {
      for (int i = 0; i < model_dims[m_k]; i++) {
        (am->initRWM)[m_k][i] = sdrand();
      }
    }
  } else {
    for (int m_k = 0; m_k < nmodels; m_k++) {
      for (int i = 0; i < model_dims[m_k]; i++) {
        (am->initRWM)[m_k][i] = initRWM[m_k][i];
      }
    }
  }
  // Set all structs as un-initialized
  (&(am->cpstats))->isInitialized = 0;
  (&(am->ch))->isInitialized = 0;
  (&(am->st))->isInitialized = 0;
  // Set default values
  am->doAdapt = 1;
  am->doPerm = 0;
  am->student_T_dof = 0;
  am->am_mixfit = FIGUEREIDO_MIX_FIT;
  return EXIT_SUCCESS;
}

void freeAMSampler(amSampler *am) {
  int nmodels = am->jd.nmodels;
  for (int i = 0; i < nmodels; i++) {
    free(am->initRWM[i]);
  }
  free(am->initRWM);
  freeCondProbStats(&(am->cpstats), am->jd);
  freeChain(&(am->ch));
  freeProposalDist(am->jd);
  return;
}

int initCondProbStats(condProbStats *cpstats, proposalDist jd, int nsweeps2,
                      int NUM_FITMIX_MAX) {
  if (cpstats->isInitialized) {
    return EXIT_SUCCESS;
  }
  cpstats->sig_k_rwm_summary =
      (double ***)malloc(jd.nmodels * sizeof(double **));
  cpstats->nacc_ntry_rwm = (double ***)malloc(jd.nmodels * sizeof(double **));
  for (int model_k = 0; model_k < jd.nmodels; model_k++) {
    int mdim = jd.model_dims[model_k];
    int nsweepr = max(nsweeps2, 10000 * mdim);
    int nburn = nsweepr / 10;
    cpstats->rwm_summary_len = max(1, (nsweepr + nburn) / 100);
    cpstats->sig_k_rwm_summary[model_k] =
        (double **)malloc(cpstats->rwm_summary_len * sizeof(double *));
    cpstats->sig_k_rwm_summary[model_k][0] =
        (double *)malloc(cpstats->rwm_summary_len * mdim * sizeof(double));
    for (int i = 1; i < cpstats->rwm_summary_len; i++) {
      cpstats->sig_k_rwm_summary[model_k][i] =
          cpstats->sig_k_rwm_summary[model_k][i - 1] + mdim;
    }
    cpstats->nacc_ntry_rwm[model_k] =
        (double **)malloc(cpstats->rwm_summary_len * sizeof(double *));
    cpstats->nacc_ntry_rwm[model_k][0] =
        (double *)malloc(cpstats->rwm_summary_len * mdim * sizeof(double));
    for (int i = 1; i < cpstats->rwm_summary_len; i++) {
      cpstats->nacc_ntry_rwm[model_k][i] =
          cpstats->nacc_ntry_rwm[model_k][i - 1] + mdim;
    }
  }
  cpstats->nfitmix = (int *)calloc(jd.nmodels, sizeof(int));
  cpstats->fitmix_annulations = (int **)malloc(jd.nmodels * sizeof(int *));
  cpstats->fitmix_costfnnew = (double **)malloc(jd.nmodels * sizeof(double *));
  cpstats->fitmix_lpn = (double **)malloc(jd.nmodels * sizeof(double *));
  cpstats->fitmix_Lkk = (int **)malloc(jd.nmodels * sizeof(int *));
  for (int i = 0; i < jd.nmodels; i++) {
    cpstats->fitmix_annulations[i] =
        (int *)malloc(NUM_FITMIX_MAX * sizeof(int));
    cpstats->fitmix_costfnnew[i] =
        (double *)malloc(NUM_FITMIX_MAX * sizeof(double));
    cpstats->fitmix_lpn[i] = (double *)malloc(NUM_FITMIX_MAX * sizeof(double));
    cpstats->fitmix_Lkk[i] = (int *)malloc(NUM_FITMIX_MAX * sizeof(int));
  }
  cpstats->isInitialized = 1;
  return EXIT_SUCCESS;
}

void freeCondProbStats(condProbStats *cpstats, proposalDist jd) {
  if (!cpstats->isInitialized) {
    return;
  }
  if (cpstats->sig_k_rwm_summary != NULL) {
    for (int model_k = 0; model_k < jd.nmodels; model_k++) {
      if (cpstats->sig_k_rwm_summary[model_k][0] != NULL) {
        free(cpstats->sig_k_rwm_summary[model_k][0]);
      }
      if (cpstats->sig_k_rwm_summary[model_k] != NULL) {
        free(cpstats->sig_k_rwm_summary[model_k]);
      }
      if (cpstats->nacc_ntry_rwm[model_k][0] != NULL) {
        free(cpstats->nacc_ntry_rwm[model_k][0]);
      }
      if (cpstats->nacc_ntry_rwm[model_k] != NULL) {
        free(cpstats->nacc_ntry_rwm[model_k]);
      }
    }
    free(cpstats->sig_k_rwm_summary);
  }
  if (cpstats->nacc_ntry_rwm != NULL) {
    free(cpstats->nacc_ntry_rwm);
  }
  if (cpstats->nfitmix != NULL) {
    free(cpstats->nfitmix);
  }
  if (cpstats->fitmix_annulations != NULL) {
    for (int i = 0; i < jd.nmodels; i++) {
      if (cpstats->fitmix_annulations[i] != NULL) {
        free(cpstats->fitmix_annulations[i]);
      }
      if (cpstats->fitmix_costfnnew[i] != NULL) {
        free(cpstats->fitmix_costfnnew[i]);
      }
      if (cpstats->fitmix_lpn[i] != NULL) {
        free(cpstats->fitmix_lpn[i]);
      }
      if (cpstats->fitmix_Lkk[i] != NULL) {
        free(cpstats->fitmix_Lkk[i]);
      }
    }
    free(cpstats->fitmix_annulations);
  }
  if (cpstats->fitmix_costfnnew != NULL) {
    free(cpstats->fitmix_costfnnew);
  }
  if (cpstats->fitmix_lpn != NULL) {
    free(cpstats->fitmix_lpn);
  }
  if (cpstats->fitmix_Lkk != NULL) {
    free(cpstats->fitmix_Lkk);
  }
  cpstats->isInitialized = 0;
}

void initRunStats(runStats *st, int nsweep, proposalDist jd) {
  if (st->isInitialized) {
    return;
  }
  st->naccrwmb = 0;
  st->ntryrwmb = 0;
  st->naccrwms = 0;
  st->ntryrwms = 0;
  st->nacctd = 0;
  st->ntrytd = 0;
  st->nsokal = 1;
  st->nkeep = (int)pow(
      2.0, min(15, (int)(log(nsweep / (2 * st->nsokal)) / log(2.0) + 0.001)));
  st->keep = nsweep - st->nkeep * st->nsokal;
  st->xr = (double *)malloc(st->nkeep * sizeof(double));
  st->ksummary = (int *)calloc(jd.nmodels, sizeof(int));
  st->pk_summary = (double **)malloc(nsweep * sizeof(double *));
  st->pk_summary[0] = (double *)malloc(nsweep * jd.nmodels * sizeof(double));
  for (int i = 1; i < nsweep; i++) {
    st->pk_summary[i] = st->pk_summary[i - 1] + jd.nmodels;
  }
  st->k_which_summary = (int *)malloc(nsweep * sizeof(int));
  st->logp_summary = (double **)malloc(nsweep * sizeof(double *));
  st->logp_summary[0] = (double *)malloc(nsweep * 2 * sizeof(double));
  for (int i = 1; i < nsweep; i++) {
    st->logp_summary[i] = st->logp_summary[i - 1] + 2;
  }
  st->theta_summary_len = (int *)calloc(jd.nmodels, sizeof(int));
  st->theta_summary_size = (int *)malloc(jd.nmodels * sizeof(int));
  int initial_size = 1000;
  for (int i = 0; i < jd.nmodels; i++) {
    st->theta_summary_size[i] = initial_size;
  }
  st->theta_summary = (double ***)malloc(jd.nmodels * sizeof(double **));
  for (int i = 0; i < jd.nmodels; i++) {
    st->theta_summary[i] = (double **)malloc(initial_size * sizeof(double *));
  }
}

void freeRunStats(runStats st, proposalDist jd) {
  if (!st.isInitialized) {
    return;
  }
  if (st.xr != NULL) {
    free(st.xr);
  }
  if (st.ksummary != NULL) {
    free(st.ksummary);
  }
  if (st.pk_summary != NULL) {
    free(st.pk_summary[0]);
  }
  if (st.pk_summary != NULL) {
    free(st.pk_summary);
  }
  if (st.k_which_summary != NULL) {
    free(st.k_which_summary);
  }
  if (st.logp_summary[0] != NULL) {
    free(st.logp_summary[0]);
  }
  if (st.logp_summary != NULL) {
    free(st.logp_summary);
  }
}

void initChain(chainState *ch, proposalDist jd, double **initRWM,
               targetFunc logposterior) {
  if (ch->isInitialized) {
    return;
  }
  ch->current_model_k = (int)floor(jd.nmodels * sdrand());
  ch->mdim = jd.model_dims[ch->current_model_k];
  int mdim_max = jd.model_dims[0];
  for (int i = 1; i < jd.nmodels; i++) {
    mdim_max = max(jd.model_dims[i], mdim_max);
  }
  ch->theta = (double *)malloc(mdim_max * sizeof(double));
  ch->pk = (double *)malloc(jd.nmodels * sizeof(double));
  for (int i = 0; i < ch->mdim; i++) {
    ch->theta[i] = initRWM[ch->current_model_k][i];
  }
  ch->current_Lkk = jd.nMixComps[ch->current_model_k];
  ch->log_posterior = logposterior(ch->current_model_k, ch->mdim, ch->theta);
  for (int i = 0; i < jd.nmodels; i++) {
    ch->pk[i] = 1.0 / jd.nmodels;
  }
  ch->nreinit = 1;
  ch->reinit = 0;
  ch->pkllim = 1.0 / 10.0;
  ch->sweep_i = 1;
  ch->isInitialized = 1;
}

void freeChain(chainState *ch) {
  if (!ch->isInitialized) {
    return;
  }
  if (ch->theta != NULL) {
    free(ch->theta);
  }
  if (ch->pk != NULL) {
    free(ch->pk);
  }
  ch->isInitialized = 0;
}

int initProposalDist(proposalDist *jd, int nmodels, int *model_dims,
                     int NUM_MIX_COMPS_MAX) {
  jd->nmodels = nmodels;
  jd->NUM_MIX_COMPS_MAX = NUM_MIX_COMPS_MAX;
  jd->model_dims = (int *)malloc(jd->nmodels * sizeof(int));
  if (jd->model_dims == NULL) {
    return EXIT_FAILURE;
  }
  for (int i = 0; i < nmodels; i++) {
    jd->model_dims[i] = model_dims[i];
  }
  jd->nMixComps = (int *)malloc(jd->nmodels * sizeof(int));
  if (jd->nMixComps == NULL) {
    return EXIT_FAILURE;
  }
  jd->lambda = (double **)malloc(jd->nmodels * sizeof(double *));
  if (jd->lambda == NULL) {
    return EXIT_FAILURE;
  }
  jd->mu = (double ***)malloc(jd->nmodels * sizeof(double **));
  if (jd->mu == NULL) {
    return EXIT_FAILURE;
  }
  jd->B = (double ****)malloc(jd->nmodels * sizeof(double ***));
  if (jd->B == NULL) {
    return EXIT_FAILURE;
  }
  for (int k = 0; k < jd->nmodels; k++) {
    int mdim = jd->model_dims[k];
    jd->lambda[k] = (double *)malloc(NUM_MIX_COMPS_MAX * sizeof(double));
    if (jd->lambda[k] == NULL) {
      return EXIT_FAILURE;
    }
    jd->mu[k] = (double **)malloc(NUM_MIX_COMPS_MAX * sizeof(double *));
    if (jd->mu[k] == NULL) {
      return EXIT_FAILURE;
    }
    jd->B[k] = (double ***)malloc(NUM_MIX_COMPS_MAX * sizeof(double **));
    if (jd->B[k] == NULL) {
      return EXIT_FAILURE;
    }
    for (int i = 0; i < NUM_MIX_COMPS_MAX; i++) {
      jd->mu[k][i] = (double *)malloc(mdim * sizeof(double));
      if (jd->mu[k][i] == NULL) {
        return EXIT_FAILURE;
      }
      jd->B[k][i] = (double **)malloc(mdim * sizeof(double *));
      if (jd->B[k][i] == NULL) {
        return EXIT_FAILURE;
      }
      for (int j = 0; j < mdim; j++) {
        jd->B[k][i][j] = (double *)malloc(mdim * sizeof(double));
        if (jd->B[k][i][j] == NULL) {
          return EXIT_FAILURE;
        }
      }
    }
  }
  jd->sig = (double **)malloc(jd->nmodels * sizeof(double *));
  for (int k = 0; k < jd->nmodels; k++) {
    int mdim = jd->model_dims[k];
    jd->sig[k] = (double *)malloc(mdim * sizeof(double));
  }

  jd->isInitialized = 1;
  return EXIT_SUCCESS;
}

void freeProposalDist(proposalDist jd) {
  for (int k = 0; k < jd.nmodels; k++) {
    int mdim = jd.model_dims[k];
    for (int i = 0; i < jd.NUM_MIX_COMPS_MAX; i++) {
      for (int j = 0; j < mdim; j++) {
        if (jd.B[k][i][j] != NULL) {
          free(jd.B[k][i][j]);
        }
      }
      if (jd.mu[k][i] != NULL) {
        free(jd.mu[k][i]);
      }
      if (jd.B[k][i] != NULL) {
        free(jd.B[k][i]);
      }
    }
    if (jd.lambda[k] != NULL) {
      free(jd.lambda[k]);
    }
    if (jd.mu[k] != NULL) {
      free(jd.mu[k]);
    }
    if (jd.B[k] != NULL) {
      free(jd.B[k]);
    }
  }
  if (jd.lambda != NULL) {
    free(jd.lambda);
  }
  if (jd.mu != NULL) {
    free(jd.mu);
  }
  if (jd.B != NULL) {
    free(jd.B);
  }
  if (jd.model_dims != NULL) {
    free(jd.model_dims);
  }
  if (jd.nMixComps != NULL) {
    free(jd.nMixComps);
  }
}

void rwm_within_model(int model_k, int *model_dims, int nsweep2,
                      condProbStats *cpstats, double *sig_k, int dof,
                      double **samples, targetFunc logpost, double *initRWM) {
  // --- Section 5.2.1 - RWM Within Model (Stage 1) -------
  int mdim = model_dims[model_k];
  int nsweepr = max(nsweep2, 10000 * mdim);
  int nburn = nsweepr / 10;
  double alphastar = 0.25;
  nsweepr += nburn;
  double *rwm = (double *)malloc(mdim * sizeof(double));
  double *rwmn = (double *)malloc(mdim * sizeof(double));
  int *nacc = (int *)malloc(mdim * sizeof(int));
  int *ntry = (int *)malloc(mdim * sizeof(int));
  double *Znkk = (double *)malloc(mdim * sizeof(double));

  int sig_k_rwm_n = 0; // This is used to load stats arrays
  for (int i = 0; i < mdim; i++) {
    rwm[i] = initRWM[i];
  }
  for (int j1 = 0; j1 < mdim; j1++) {
    rwmn[j1] = rwm[j1];
    sig_k[j1] = 10.0;
    nacc[j1] = 0;
    ntry[j1] = 0;
  }
  double lp = logpost(model_k, mdim, rwm);

  int i2 = 0;
  int remain = nsweepr;
  for (int sweep = 1; sweep <= nsweepr; sweep++) {
    remain--;
    double u = sdrand();
    if (sweep > nburn && u < 0.1) {
      rt(Znkk, mdim, dof);
      for (int i = 0; i < mdim; i++) {
        rwmn[i] = rwm[i] + sig_k[i] * Znkk[i];
      }
      double lpn = logpost(model_k, mdim, rwmn);
      if (sdrand() < exp(max(-30.0, min(0.0, lpn - lp)))) {
        for (int i = 0; i < mdim; i++) {
          rwm[i] = rwmn[i];
        }
        lp = lpn;
      }
    } else {
      double gamma = 10.0 * pow(1.0 / (sweep + 1), 2.0 / 3.0);
      for (int i = 0; i < mdim; i++) {
        rwmn[i] = rwm[i];
      }
      for (int i = 0; i < mdim; i++) {
        double Z;
        rt(&Z, 1, dof);
        rwmn[i] = rwm[i] + sig_k[i] * Z;
        double lpn = logpost(model_k, mdim, rwmn);
        double accept = min(1, exp(max(-30.0, min(0.0, lpn - lp))));
        if (sdrand() < accept) {
          (nacc[i])++;
          (ntry[i])++;
          rwm[i] = rwmn[i];
          lp = lpn;
          sig_k[i] = max(0, sig_k[i] - gamma * (alphastar - 1));
        } else {
          (ntry[i])++;
          rwmn[i] = rwm[i];
          sig_k[i] = max(0, sig_k[i] - gamma * (alphastar));
        }
      }
    }
    if (remain < (10000 * mdim) && (remain % 10 == 0)) {
      for (int i = 0; i < mdim; i++) {
        samples[i2][i] = rwm[i];
      }
      i2++;
    }
    if (sweep % 100 == 0) {
      for (int i = 0; i < mdim; i++) {
        cpstats->sig_k_rwm_summary[model_k][sig_k_rwm_n][i] = sig_k[i];
        cpstats->nacc_ntry_rwm[model_k][sig_k_rwm_n][i] =
            (double)nacc[i] / (double)ntry[i];
      }
      sig_k_rwm_n++;
    }
  }
  free(rwm);
  free(rwmn);
  free(nacc);
  free(ntry);
  free(Znkk);
}

void fit_mixture_from_samples(int model_k, proposalDist jd, double **samples,
                              int nsamples, condProbStats *cpstats,
                              int NUM_MIX_COMPS_MAX, int NUM_FITMIX_MAX) {
  // --- Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)-
  // Mixture fitting done component wise EM algorithm described in
  // Figueiredo and Jain, 2002 (see thesis for full reference)
  int mdim = jd.model_dims[model_k];
  double *lambda_k = jd.lambda[model_k];
  double **mu_k = jd.mu[model_k];
  double ***B_k = jd.B[model_k];
  int Lkmaxmax = NUM_MIX_COMPS_MAX;
  int Lkk = Lkmaxmax;
  int *init = (int *)malloc(Lkk * sizeof(int));
  int l1 = 0;
  double lpn = 0.0;
  double costfn = 0.0;
  double costfnmin = 0.0;

  while (l1 < Lkk) {
    int indic = 0;
    double u = sdrand();
    init[l1] = (int)floor(nsamples * u);
    if (l1 > 0) {
      for (int l2 = 0; l2 < l1; l2++) {
        if (init[l2] == init[l1]) {
          indic = 1;
          break;
        }
      }
    }
    if (indic == 0) {
      l1++;
    }
  }

  // sigma holds the trace of the covariance matrix plus some normalization.
  double sigma = 0.0;
  for (int j = 0; j < mdim; j++) {
    double datasum = 0.0;
    double datasqsum = 0.0;
    for (int i = 0; i < nsamples; i++) {
      datasum += samples[i][j];
      datasqsum += samples[i][j] * samples[i][j];
    }
    double len = (double)nsamples;
    sigma += (datasqsum - datasum * datasum / len) / len;
  }
  sigma /= (10.0 * mdim);

  for (int i = 0; i < Lkk; i++) {
    for (int j = 0; j < mdim; j++) {
      mu_k[i][j] = samples[init[i]][j];
      B_k[i][j][j] = sigma;
      for (int k = 0; k < j; k++) {
        B_k[i][j][k] = 0.0;
      }
    }
    chol(mdim, B_k[i]);
    lambda_k[i] = 1.0 / Lkk;
  }

  double **w = (double **)malloc(nsamples * sizeof(double *));
  double *logw = (double *)malloc(Lkk * sizeof(double));
  double **lpdatagivenl = (double **)malloc(nsamples * sizeof(double *));
  for (int i1 = 0; i1 < nsamples; i1++) {
    w[i1] = (double *)malloc(Lkk * sizeof(double));
    lpdatagivenl[i1] = (double *)malloc(Lkk * sizeof(double));
  }

  for (int i = 0; i < nsamples; i++) {
    double sum = 0.0;
    for (int j = 0; j < Lkk; j++) {
      lpdatagivenl[i][j] = lnormprob(mdim, mu_k[j], B_k[j], samples[i]);
      logw[j] = log(lambda_k[j]) + lpdatagivenl[i][j];
      w[i][j] = exp(logw[j]);
      sum += w[i][j];
    }
    for (int j = 0; j < Lkk; j++) {
      w[i][j] /= sum;
    }
  }

  double *sumw = (double *)malloc(Lkk * sizeof(double));

  int stop = 0;
  int count = 0;
  double wnewl1 = 0.0;
  int nparams = mdim + (mdim * (mdim + 1)) / 2;
  int Lkkmin = 0;

  double ***Bmin = (double ***)malloc(Lkmaxmax * sizeof(double **));
  for (int l1 = 0; l1 < Lkmaxmax; l1++) {
    Bmin[l1] = (double **)malloc(mdim * sizeof(double *));
    for (int j1 = 0; j1 < mdim; j1++) {
      Bmin[l1][j1] = (double *)malloc(mdim * sizeof(double));
    }
  }
  double **mumin = (double **)malloc(Lkmaxmax * sizeof(double *));
  for (int l1 = 0; l1 < Lkmaxmax; l1++) {
    mumin[l1] = (double *)malloc(mdim * sizeof(double));
  }
  double *lambdamin = (double *)malloc(Lkmaxmax * sizeof(double));

  while (!stop) {
    count++;
    l1 = 0;
    int natann = 0;
    int forceann = 0;
    while (l1 < Lkk) {
      double sumwnew = 0.0;
      for (int l2 = 0; l2 < Lkk; l2++) {
        sumw[l2] = 0.0;
        for (int i1 = 0; i1 < nsamples; i1++) {
          sumw[l2] += w[i1][l2];
        }
        double wnew = max(0.0, (sumw[l2] - nparams / 2.0));
        if (l2 == l1) {
          wnewl1 = wnew;
        }
        sumwnew += wnew;
      }
      lambda_k[l1] = wnewl1 / sumwnew;
      double sumlambda = 0.0;
      for (int l2 = 0; l2 < Lkk; l2++) {
        sumlambda += lambda_k[l2];
      }
      for (int l2 = 0; l2 < Lkk; l2++) {
        lambda_k[l2] /= sumlambda;
      }

      if (lambda_k[l1] > 0.005) {
        // changed to 0.005 from 0.0 -renormalise else
        for (int j1 = 0; j1 < mdim; j1++) {
          mu_k[l1][j1] = 0.0;
          for (int i1 = 0; i1 < nsamples; i1++) {
            mu_k[l1][j1] += samples[i1][j1] * w[i1][l1];
          }
          mu_k[l1][j1] /= sumw[l1];

          for (int j2 = 0; j2 <= j1; j2++) {
            double temp = 0.0;
            for (int i1 = 0; i1 < nsamples; i1++) {
              temp += (samples[i1][j1] - mu_k[l1][j1]) *
                      (samples[i1][j2] - mu_k[l1][j2]) * w[i1][l1];
            }
            B_k[l1][j1][j2] = temp / sumw[l1];
          }
        }

        chol(mdim, B_k[l1]);

        for (int i1 = 0; i1 < nsamples; i1++) {
          lpdatagivenl[i1][l1] =
              lnormprob(mdim, mu_k[l1], B_k[l1], samples[i1]);
        }
        l1++;

      } else {
        natann = 1;
        if (l1 < (Lkk - 1)) {
          for (int l2 = l1; l2 < (Lkk - 1); l2++) {
            lambda_k[l2] = lambda_k[l2 + 1];
            for (int j1 = 0; j1 < mdim; j1++) {
              mu_k[l2][j1] = mu_k[l2 + 1][j1];
              for (int j2 = 0; j2 <= j1; j2++) {
                B_k[l2][j1][j2] = B_k[l2 + 1][j1][j2];
              }
            }
            for (int i1 = 0; i1 < nsamples; i1++) {
              lpdatagivenl[i1][l2] = lpdatagivenl[i1][l2 + 1];
            }
          }
        }
        Lkk--;
        sumlambda = 0.0;
        for (int l2 = 0; l2 < Lkk; l2++) {
          sumlambda += lambda_k[l2];
        }
        for (int l2 = 0; l2 < Lkk; l2++) {
          lambda_k[l2] /= sumlambda;
        }
      }

      lpn = 0.0;
      for (int i1 = 0; i1 < nsamples; i1++) {
        double sum = 0.0;
        for (int l2 = 0; l2 < Lkk; l2++) {
          logw[l2] = log(lambda_k[l2]) + lpdatagivenl[i1][l2];
          w[i1][l2] = exp(logw[l2]);
          sum += w[i1][l2];
        }
        if (sum > 0) {
          for (int l2 = 0; l2 < Lkk; l2++) {
            w[i1][l2] /= sum;
          }
          lpn += log(sum);
        } else {
          // if no component fits point well make equally likely
          for (int l2 = 0; l2 < Lkk; l2++) {
            w[i1][l2] = 1.0 / Lkk;
          }
          lpn -= 500.0;
        }
      }
    }

    double sum = 0.0;
    for (int l1 = 0; l1 < Lkk; l1++) {
      sum += log(nsamples * lambda_k[l1] / 12.0);
    }
    double costfnnew = (nparams / 2.0) * sum +
                       (Lkk / 2.0) * log(nsamples / 12.0) +
                       Lkk * (nparams + 1) / 2.0 - lpn;

    if (count == 1) {
      costfn = costfnnew;
    }
    if (count == 1 || costfnnew < costfnmin) {
      Lkkmin = Lkk;
      costfnmin = costfnnew;
      for (int l1 = 0; l1 < Lkk; l1++) {
        lambdamin[l1] = lambda_k[l1];
        for (int j1 = 0; j1 < mdim; j1++) {
          mumin[l1][j1] = mu_k[l1][j1];
          for (int j2 = 0; j2 <= j1; j2++) {
            Bmin[l1][j1][j2] = B_k[l1][j1][j2];
          }
        }
      }
    }
    if ((fabs(costfn - costfnnew) < min(1E-5 * fabs(costfn), 0.01)) &&
        (count > 1)) {
      if (Lkk == 1) {
        stop = 1;
      } else {
        forceann = 2;
        double minlambda = lambda_k[0];
        int ldel = 0;
        for (int l1 = 1; l1 < Lkk; l1++) {
          if (minlambda > lambda_k[l1]) {
            minlambda = lambda_k[l1];
            ldel = l1;
          }
        }
        if (ldel < (Lkk - 1)) {
          for (int l1 = ldel; l1 < (Lkk - 1); l1++) {
            lambda_k[l1] = lambda_k[l1 + 1];
            for (int j1 = 0; j1 < mdim; j1++) {
              mu_k[l1][j1] = mu_k[l1 + 1][j1];
              for (int j2 = 0; j2 <= j1; j2++) {
                B_k[l1][j1][j2] = B_k[l1 + 1][j1][j2];
              }
            }
            for (int i1 = 0; i1 < nsamples; i1++) {
              lpdatagivenl[i1][l1] = lpdatagivenl[i1][l1 + 1];
            }
          }
        }
        Lkk--;
        double sumlambda = 0.0;
        for (int l1 = 0; l1 < Lkk; l1++) {
          sumlambda += lambda_k[l1];
        }
        for (int l1 = 0; l1 < Lkk; l1++) {
          lambda_k[l1] /= sumlambda;
        }

        lpn = 0.0;
        for (int i1 = 0; i1 < nsamples; i1++) {
          sum = 0.0;
          for (int l2 = 0; l2 < Lkk; l2++) {
            logw[l2] = log(lambda_k[l2]) + lpdatagivenl[i1][l2];
            w[i1][l2] = exp(logw[l2]);
            sum += w[i1][l2];
          }
          if (sum > 0) {
            for (int l2 = 0; l2 < Lkk; l2++) {
              w[i1][l2] /= sum;
            }
            lpn += log(sum);
          } else {
            // if no component fits point well make equally likely
            for (int l2 = 0; l2 < Lkk; l2++) {
              w[i1][l2] = 1.0 / Lkk;
            }
            lpn -= 500.0;
          }
        }

        sum = 0.0;
        for (int l1 = 0; l1 < Lkk; l1++) {
          sum += log(nsamples * lambda_k[l1] / 12.0);
        }
        costfnnew = (nparams / 2.0) * sum + (Lkk / 2.0) * log(nsamples / 12.0) +
                    Lkk * (nparams + 1) / 2.0 - lpn;
      }
    }
    if (count > NUM_FITMIX_MAX) {
      stop = 1;
    }
    costfn = costfnnew;
    int nfitmix = cpstats->nfitmix[model_k];
    cpstats->fitmix_annulations[model_k][nfitmix] = natann + forceann;
    cpstats->fitmix_costfnnew[model_k][nfitmix] = costfnnew;
    cpstats->fitmix_lpn[model_k][nfitmix] = lpn;
    cpstats->fitmix_Lkk[model_k][nfitmix] = Lkk;
    (cpstats->nfitmix[model_k])++;
  }

  for (int i = 0; i < nsamples; i++) {
    free(w[i]);
    free(lpdatagivenl[i]);
  }
  free(w);
  free(lpdatagivenl);
  free(logw);
  free(sumw);
  free(init);
  jd.nMixComps[model_k] = Lkkmin;
  for (int i = 0; i < Lkkmin; i++) {
    lambda_k[i] = lambdamin[i];
    for (int j = 0; j < mdim; j++) {
      mu_k[i][j] = mumin[i][j];
    }
    for (int j = 0; j < mdim; j++) {
      for (int l = 0; l <= j; l++) {
        B_k[i][j][l] = Bmin[i][j][l];
      }
    }
  }
  for (int i = 0; i < Lkmaxmax; i++) {
    for (int j = 0; j < mdim; j++) {
      free(Bmin[i][j]);
    }
    free(Bmin[i]);
  }
  free(Bmin);
  for (int i = 0; i < Lkmaxmax; i++) {
    free(mumin[i]);
  }
  free(mumin);
  free(lambdamin);
}

void fit_autorj(int model_k, proposalDist jd, double **samples, int nsamples) {
  // --- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --
  jd.nMixComps[model_k] = 1;
  jd.lambda[model_k][0] = 1.0;
  int mdim = jd.model_dims[model_k];
  double **mu_k = jd.mu[model_k];
  double ***B_k = jd.B[model_k];
  for (int j = 0; j < mdim; j++) {
    mu_k[0][j] = 0.0;
    for (int i = 0; i < nsamples; i++) {
      mu_k[0][j] += samples[i][j];
    }
    mu_k[0][j] /= ((double)nsamples);
  }
  for (int l = 0; l < mdim; l++) {
    for (int j = 0; j <= l; j++) {
      B_k[0][l][j] = 0.0;
      for (int i = 0; i < nsamples; i++) {
        B_k[0][l][j] +=
            (samples[i][l] - mu_k[0][l]) * (samples[i][j] - mu_k[0][j]);
      }
      B_k[0][l][j] /= ((double)(nsamples - 1));
    }
  }
  chol(mdim, B_k[0]);
}

void reversible_jump_move(bool doPerm, bool doAdapt, chainState *ch,
                          proposalDist jd, int dof, runStats *st,
                          targetFunc logpost) {
  int Lkmax = jd.nMixComps[0];
  for (int k1 = 1; k1 < jd.nmodels; k1++) {
    Lkmax = max(Lkmax, jd.nMixComps[k1]);
  }
  int mdim_max = jd.model_dims[0];
  for (int i = 1; i < jd.nmodels; i++) {
    mdim_max = max(jd.model_dims[i], mdim_max);
  }
  double *theta = ch->theta;
  double *thetan = (double *)malloc(mdim_max * sizeof(double));
  double *work = (double *)malloc(mdim_max * sizeof(double));
  double *palloc = (double *)malloc(Lkmax * sizeof(double));
  double *pallocn = (double *)malloc(Lkmax * sizeof(double));
  double *Znkk = (double *)malloc(mdim_max * sizeof(double));
  const double logrtpi = 0.9189385332046727; // 0.5 * log(2.0 * pi)

  // --Section 8 - RWM within-model moves ---
  // Every 10 sweeps to block RWM */
  if (ch->doBlockRWM) {
    (st->ntryrwmb)++;
    rt(Znkk, ch->mdim, dof);
    for (int i = 0; i < ch->mdim; i++) {
      thetan[i] = theta[i] + jd.sig[ch->current_model_k][i] * Znkk[i];
    }
    double lpn = logpost(ch->current_model_k, ch->mdim, thetan);
    if (sdrand() < exp(max(-30.0, min(0.0, lpn - ch->log_posterior)))) {
      (st->naccrwmb)++;
      memcpy(theta, thetan, ch->mdim * sizeof(*thetan));
      ch->log_posterior = lpn;
    }
  } else {
    // else do component-wise RWM
    memcpy(thetan, theta, (ch->mdim) * sizeof(*thetan));
    for (int j1 = 0; j1 < ch->mdim; j1++) {
      (st->ntryrwms)++;
      double Z;
      rt(&Z, 1, dof);
      thetan[j1] = theta[j1] + jd.sig[ch->current_model_k][j1] * Z;
      double lpn = logpost(ch->current_model_k, ch->mdim, thetan);
      if (sdrand() < exp(max(-30.0, min(0.0, lpn - ch->log_posterior)))) {
        (st->naccrwms)++;
        theta[j1] = thetan[j1];
        ch->log_posterior = lpn;
      } else {
        thetan[j1] = theta[j1];
      }
    }
  }

  // --- Section 9 Reversible Jump Moves -------------------

  // --Section 9.1 - Allocate current position to a component --
  int l = 0;
  (st->ntrytd)++;
  int ln = 0;
  if (ch->current_Lkk > 1) {
    double sum = 0.0;
    for (int i = 0; i < ch->current_Lkk; i++) {
      palloc[i] = log(jd.lambda[ch->current_model_k][i]) +
                  lnormprob(ch->mdim, jd.mu[ch->current_model_k][i],
                            jd.B[ch->current_model_k][i], theta);
      palloc[i] = exp(palloc[i]);
      sum += palloc[i];
    }
    if (sum > 0) {
      for (int i = 0; i < ch->current_Lkk; i++) {
        palloc[i] /= sum;
      }
    } else {
      for (int i = 0; i < ch->current_Lkk; i++) {
        palloc[i] = 1.0 / ch->current_Lkk;
      }
    }
    double u = sdrand();
    double thresh = 0.0;
    for (int i = 0; i < ch->current_Lkk; i++) {
      thresh += palloc[i];
      if (u < thresh) {
        l = i;
        break;
      }
    }
  } else {
    l = 0;
    palloc[l] = 1.0;
  }

  // --Section 9.2 - Standardise state variable ---------------

  for (int i = 0; i < ch->mdim; i++) {
    work[i] = theta[i] - jd.mu[ch->current_model_k][l][i];
  }
  for (int i = 0; i < ch->mdim; i++) {
    for (int j = 0; j < i; j++) {
      work[i] = work[i] - jd.B[ch->current_model_k][l][i][j] * work[j];
    }
    work[i] = work[i] / jd.B[ch->current_model_k][l][i][i];
  }

  // --Section 9.3 - Choose proposed new model and component ----
  int kn = 0;
  double gamma = 0.0;
  double logratio;
  if (jd.nmodels == 1) {
    kn = ch->current_model_k;
    logratio = 0.0;
  } else {
    gamma = pow(1.0 / (ch->sweep_i + 1), (2.0 / 3.0));
    double u = sdrand();
    double thresh = 0.0;
    for (int k1 = 0; k1 < jd.nmodels; k1++) {
      thresh += ch->pk[k1];
      if (u < thresh) {
        kn = k1;
        break;
      }
    }
    logratio = log(ch->pk[ch->current_model_k]) - log(ch->pk[kn]);
  }

  int mdim_kn = jd.model_dims[kn];
  int Lkkn = jd.nMixComps[kn];

  double u = sdrand();
  double thresh = 0.0;
  for (int l1 = 0; l1 < Lkkn; l1++) {
    thresh += jd.lambda[kn][l1];
    if (u < thresh) {
      ln = l1;
      break;
    }
  }

  // --Section 9.4 Propose new state ----------------

  if (ch->mdim < mdim_kn) {
    rt(&(work[ch->mdim]), mdim_kn - ch->mdim, dof);
    if (dof > 0) {
      for (int i = ch->mdim; i < mdim_kn; i++) {
        logratio -= ltprob(dof, work[i]);
      }
    } else {
      for (int j1 = ch->mdim; j1 < mdim_kn; j1++) {
        logratio += 0.5 * pow(work[j1], 2.0) + logrtpi;
      }
    }
    if (doPerm) {
      perm(work, mdim_kn);
    }
  } else if (ch->mdim == mdim_kn) {
    if (doPerm) {
      perm(work, ch->mdim);
    }
  } else {
    if (doPerm) {
      perm(work, ch->mdim);
    }
    if (dof > 0) {
      for (int j1 = mdim_kn; j1 < ch->mdim; j1++) {
        logratio += ltprob(dof, work[j1]);
      }
    } else {
      for (int j1 = mdim_kn; j1 < ch->mdim; j1++) {
        logratio -= (0.5 * pow(work[j1], 2.0) + logrtpi);
      }
    }
  }

  for (int j1 = 0; j1 < mdim_kn; j1++) {
    thetan[j1] = jd.mu[kn][ln][j1];
    for (int j2 = 0; j2 <= j1; j2++) {
      thetan[j1] += jd.B[kn][ln][j1][j2] * work[j2];
    }
  }

  // --Section 9.5 - Work out probability of allocating to component
  // for acceptance ratio (reverse move) ---------

  if (Lkkn > 1) {
    double sum = 0.0;
    for (int l1 = 0; l1 < Lkkn; l1++) {
      pallocn[l1] = log(jd.lambda[kn][l1]) +
                    lnormprob(mdim_kn, jd.mu[kn][l1], jd.B[kn][l1], thetan);
      pallocn[l1] = exp(pallocn[l1]);
      sum += pallocn[l1];
    }
    if (sum > 0) {
      for (int l1 = 0; l1 < Lkkn; l1++) {
        pallocn[l1] /= sum;
      }
    } else {
      for (int l1 = 0; l1 < Lkkn; l1++) {
        pallocn[l1] = 1.0 / Lkkn;
      }
    }
  } else {
    pallocn[ln] = 1.0;
  }

  // --Section 9.6 - Work out acceptance probability  and new state --
  double lpn = logpost(kn, mdim_kn, thetan);

  int mdim = jd.model_dims[ch->current_model_k];
  logratio += (lpn - ch->log_posterior);
  logratio += (log(pallocn[ln]) - log(palloc[l]));
  logratio += (log(jd.lambda[ch->current_model_k][l]) - log(jd.lambda[kn][ln]));
  logratio += (log(det(mdim_kn, jd.B[kn][ln])) -
               log(det(mdim, jd.B[ch->current_model_k][l])));

  if (sdrand() < exp(max(-30.0, min(0.0, logratio)))) {
    for (int j1 = 0; j1 < mdim_kn; j1++) {
      theta[j1] = thetan[j1];
    }
    ch->log_posterior = lpn;
    ch->current_model_k = kn;
    ch->mdim = mdim_kn;
    ch->current_Lkk = Lkkn;
    (st->nacctd)++;
  }

  if (doAdapt && !ch->isBurning) {
    for (int k1 = 0; k1 < jd.nmodels; k1++) {
      double propk;
      if (k1 == ch->current_model_k) {
        propk = 1.0;
      } else {
        propk = 0.0;
      }
      ch->pk[k1] += (gamma * (propk - ch->pk[k1]));
    }
    for (int k1 = 0; k1 < jd.nmodels; k1++) {
      if (ch->pk[k1] < ch->pkllim) {
        ch->reinit = 1;
        break;
      }
    }
    if (ch->reinit == 1) {
      ch->reinit = 0;
      (ch->nreinit)++;
      ch->pkllim = 1.0 / (10.0 * ch->nreinit);
      for (int k1 = 0; k1 < jd.nmodels; k1++) {
        ch->pk[k1] = 1.0 / jd.nmodels;
      }
    }
  }
  free(thetan);
  free(work);
  free(palloc);
  free(pallocn);
  free(Znkk);
}

/* Combined congruential and Tauseworthe generators from SuperDuper
 * package. Should work on machines with unsigned long of at least 32
 * bits. JC and JT must be initialized to values with 0 < JC < 2^32 and
 * 0 < JT < 2^32. JC must be odd.
 * References: Marsaglia, Ananthanarayanan & Paul, 1973,
 * Learmonth & Lewis, 1973, and Dudewicz, 1976)
 */
static unsigned long JC, JT;
static double Norm = 4.656612873E-10;

double sdrand(void) {
  JC = (JC * 69069) & 037777777777; /* congruential part */
  JT ^= JT >> 15;                   /* tausworthe part */
  JT ^= (JT << 17) & 037777777777;
  return (((JT ^ JC) >> 1) * Norm);
}

void sdrni(unsigned long *i) {
  unsigned long k = *i;
  if (k == 0)
    k = time(0);
  JT = k / 65536;
  JC = k - 65536 * JT;
  JT = 65536 * JT + 1;
  JC = 32768 * JC + 1;
  *i = k;
}

/* Routines calculating quantities related to the gamma function */

/* Taken from algama.f (PJG) - converted to C using appropriate machine
   constants by DIH 04/11/03 */

double loggamma(double X) {

  /* This routine calculates the LOG GAMMA function for a positive real
     argument X.  Computation is based on an algorithm outlined in
     references 1 and 2.  The program uses rational functions that
     theoretically approximate LOG GAMMA to at least 18 significant
     decimal digits.  The approximation for X > 12 is from reference
     3, while approximations for X < 12.0 are similar to those in
     reference 1, but are unpublished.  The accuracy achieved depends
     on the arithmetic system, the compiler, the intrinsic functions,
     and proper selection of the machine-dependent constants.

     Values taken from float.h and algama.f (for XBIG)

     ---------------------------------------------------------------

     Explanation of machine-dependent constants

  beta   - radix for the floating-point representation
  maxexp - the smallest positive power of beta that overflows
  XBIG   - largest argument for which LN(GAMMA(X)) is representable
           in the machine, i.e., the solution to the equation
                   LN(GAMMA(XBIG)) = beta**maxexp
  XINF   - largest machine representable floating-point number;
           approximately beta**maxexp.
  EPS    - The smallest positive floating-point number such that
           1.0+EPS > 1.0
  FRTBIG - Rough estimate of the fourth root of XBIG

     ---------------------------------------------------------------

     Error returns

     The program returns the value XINF for X <= 0.0 or when
     overflow would occur.  The computation is believed to
     be free of underflow and overflow.

     ---------------------------------------------------------------
     References:

     1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
     1967, pp. 198-203.

     2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
     1969.

     3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
     York, 1968.

  -----------------------------------------------------------------
  Start of code
  -----------------------------------------------------------------*/

  int I;
  double CORR, D1, D2, D4, EPS, FRTBIG, FOUR, HALF, ONE, PNT68;
  double RES, SQRTPI, THRHAL, TWELVE, TWO, XBIG, XDEN, XINF;
  double XM1, XM2, XM4, XNUM, Y, YSQ, ZERO;
  double C[7], P1[8], P2[8], P4[8], Q1[8], Q2[8], Q4[8];

  /*---------------------------------------------------------------------
  Mathematical constants
  ---------------------------------------------------------------------*/

  ONE = 1.0E0;
  HALF = 0.5E0;
  TWELVE = 12.0E0;
  ZERO = 0.0E0;
  FOUR = 4.0E0;
  THRHAL = 1.5E0;
  TWO = 2.0E0;
  PNT68 = 0.6796875E0;
  SQRTPI = 0.9189385332046727417803297E0; /* eh? */

  /*---------------------------------------------------------------------
    Machine dependent parameters
    -------------------------------------------------------------------*/

  XBIG = 2.55E305;
  XINF = 1.79E308;
  EPS = 2.22E-16;
  FRTBIG = 2.25E76;

  /*--------------------------------------------------------------------
    Numerator and denominator coefficients for rational minimax
    approximation over (EPS,1.5).
    -------------------------------------------------------------------*/
  D1 = -5.772156649015328605195174E-1;
  P1[0] = 4.945235359296727046734888E0;
  P1[1] = 2.018112620856775083915565E2;
  P1[2] = 2.290838373831346393026739E3;
  P1[3] = 1.131967205903380828685045E4;
  P1[4] = 2.855724635671635335736389E4;
  P1[5] = 3.848496228443793359990269E4;
  P1[6] = 2.637748787624195437963534E4;
  P1[7] = 7.225813979700288197698961E3;
  Q1[0] = 6.748212550303777196073036E1;
  Q1[1] = 1.113332393857199323513008E3;
  Q1[2] = 7.738757056935398733233834E3;
  Q1[3] = 2.763987074403340708898585E4;
  Q1[4] = 5.499310206226157329794414E4;
  Q1[5] = 6.161122180066002127833352E4;
  Q1[6] = 3.635127591501940507276287E4;
  Q1[7] = 8.785536302431013170870835E3;

  /*---------------------------------------------------------------------
     Numerator and denominator coefficients for rational minimax
     Approximation over (1.5,4.0).
     ------------------------------------------------------------------*/

  D2 = 4.227843350984671393993777E-1;
  P2[0] = 4.974607845568932035012064E0;
  P2[1] = 5.424138599891070494101986E2;
  P2[2] = 1.550693864978364947665077E4;
  P2[3] = 1.847932904445632425417223E5;
  P2[4] = 1.088204769468828767498470E6;
  P2[5] = 3.338152967987029735917223E6;
  P2[6] = 5.106661678927352456275255E6;
  P2[7] = 3.074109054850539556250927E6;
  Q2[0] = 1.830328399370592604055942E2;
  Q2[1] = 7.765049321445005871323047E3;
  Q2[2] = 1.331903827966074194402448E5;
  Q2[3] = 1.136705821321969608938755E6;
  Q2[4] = 5.267964117437946917577538E6;
  Q2[5] = 1.346701454311101692290052E7;
  Q2[6] = 1.782736530353274213975932E7;
  Q2[7] = 9.533095591844353613395747E6;

  /*--------------------------------------------------------------------
    Numerator and denominator coefficients for rational minimax
    Approximation over (4.0,12.0).
    -------------------------------------------------------------------*/

  D4 = 1.791759469228055000094023E0;
  P4[0] = 1.474502166059939948905062E4;
  P4[1] = 2.426813369486704502836312E6;
  P4[2] = 1.214755574045093227939592E8;
  P4[3] = 2.663432449630976949898078E9;
  P4[4] = 2.940378956634553899906876E10;
  P4[5] = 1.702665737765398868392998E11;
  P4[6] = 4.926125793377430887588120E11;
  P4[7] = 5.606251856223951465078242E11;
  Q4[0] = 2.690530175870899333379843E3;
  Q4[1] = 6.393885654300092398984238E5;
  Q4[2] = 4.135599930241388052042842E7;
  Q4[3] = 1.120872109616147941376570E9;
  Q4[4] = 1.488613728678813811542398E10;
  Q4[5] = 1.016803586272438228077304E11;
  Q4[6] = 3.417476345507377132798597E11;
  Q4[7] = 4.463158187419713286462081E11;

  /*---------------------------------------------------------------------
    Coefficients for minimax approximation over (12, INF).
    -------------------------------------------------------------------*/
  C[0] = -1.910444077728E-03;
  C[1] = 8.4171387781295E-04;
  C[2] = -5.952379913043012E-04;
  C[3] = 7.93650793500350248E-04;
  C[4] = -2.777777777777681622553E-03;
  C[5] = 8.333333333333333331554247E-02;
  C[6] = 5.7083835261E-03;

  /*----------------------------------------------------------------------
    0 < X <= EPS
    --------------------------------------------------------------------*/
  Y = X;
  if ((Y > 0) && (Y <= XBIG)) {
    if (Y <= EPS) {
      RES = -log(Y);
    } else if (Y <= THRHAL) {
      /*-----------------------------------------------------------------------
        EPS < X <= 1.5
        ---------------------------------------------------------------------*/
      if (Y < PNT68) {
        CORR = -log(Y);
        XM1 = Y;
      } else {
        CORR = ZERO;
        XM1 = (Y - HALF) - HALF;
      }

      if ((Y <= HALF) || (Y >= PNT68)) {
        XDEN = ONE;
        XNUM = ZERO;
        for (I = 0; I < 8; I++) {
          XNUM = XNUM * XM1 + P1[I];
          XDEN = XDEN * XM1 + Q1[I];
        }
        RES = CORR + (XM1 * (D1 + XM1 * (XNUM / XDEN)));
      } else {
        XM2 = (Y - HALF) -
              HALF; /*Is XM2 symbol used to agree with other 2 symbols*/
        XDEN = ONE;
        XNUM = ZERO;
        for (I = 0; I < 8; I++) {
          XNUM = XNUM * XM2 + P2[I];
          XDEN = XDEN * XM2 + Q2[I];
        }
        RES = CORR + XM2 * (D2 + XM2 * (XNUM / XDEN));
      }
    } else if (Y <= FOUR) {

      /*---------------------------------------------------------------------
        1.5 < X <= 4.0
        -------------------------------------------------------------------*/
      XM2 = Y - TWO;
      XDEN = ONE;
      XNUM = ZERO;
      for (I = 0; I < 8; I++) {
        XNUM = XNUM * XM2 + P2[I];
        XDEN = XDEN * XM2 + Q2[I];
      }
      RES = XM2 * (D2 + XM2 * (XNUM / XDEN));
    } else if (Y <= TWELVE) {

      /*----------------------------------------------------------------------
        4.0 < X <= 12.0
        --------------------------------------------------------------------*/
      XM4 = Y - FOUR;
      XDEN = -ONE;
      XNUM = ZERO;
      for (I = 0; I < 8; I++) {
        XNUM = XNUM * XM4 + P4[I];
        XDEN = XDEN * XM4 + Q4[I];
      }
      RES = D4 + XM4 * (XNUM / XDEN);
    } else {

      /*---------------------------------------------------------------------
        X > 12.0,
        -------------------------------------------------------------------*/
      RES = ZERO;
      if (Y <= FRTBIG) {
        RES = C[6];
        YSQ = Y * Y;
        for (I = 0; I < 6; I++) {
          RES = RES / YSQ + C[I];
        }
      }
      RES = RES / Y;
      CORR = log(Y);
      RES = RES + SQRTPI - HALF * CORR;
      RES = RES + Y * (CORR - ONE);
    }
  } else {

    /*----------------------------------------------------------------------
      Return for bad arguments
      --------------------------------------------------------------------*/
    RES = XINF;
  }

  /*----------------------------------------------------------------------
    Final return
    --------------------------------------------------------------------*/
  return (RES);
}

/* Function for generating a Gamma(s,1) random variable
   Converted to C from PJG rgamma FORTRAN function
   Note: (1/t)*GAMMA(s,1) is GAMMA(s,t) */

double rgamma(double s) {

  double exp1, b, c1, c2, c3, c4, c5, u1, u2, w, bu, out;

  exp1 = exp(1.0);

  if (s < 1) {
    b = (s + exp1) / exp1;
    c1 = 1.0 / s;
  LAB_1:
    bu = b * sdrand();
    if (bu <= 1.0) {
      double max = (-30.0) > (c1 * log(bu)) ? (-30.0) : (c1 * log(bu));
      out = exp(max);
      if (sdrand() >= exp(-out)) {
        goto LAB_1;
      }
    } else {
      out = -log((b - bu) / s);
      if (sdrand() >= pow(out, (s - 1.0))) {
        goto LAB_1;
      }
    }
  } else if (s == 1.0) {
    out = -log(sdrand());
  } else {
    c1 = s - 1.0;
    c2 = (s - 1.0 / (6.0 * s)) / c1;
    c3 = 2.0 / c1;
    c4 = c3 + 2.0;
    c5 = 1.0 / sqrt(s);
  LAB_2:
    u1 = sdrand();
    u2 = sdrand();
    if (s > 2.5) {
      u1 = u2 + c5 * (1.0 - 1.86 * u1);
    }
    if (u1 <= 0.0 || u1 >= 1.0) {
      goto LAB_2;
    }
    w = c2 * u2 / u1;
    if ((c3 * u1 + w + 1.0 / w) <= c4) {
      goto LAB_3;
    }
    if ((c3 * log(u1) - log(w) + w) >= 1.0) {
      goto LAB_2;
    }
  LAB_3:
    out = c1 * w;
  }

  return out;
}

void gauss(double *z, int n) {
  /* Uses Box mueller method to simulate n N(0,1) variables and stores them
     in z */

  int n1;
  double u, v;

  n1 = n - 1;
  for (int i = 0; i < n1; i += 2) {
    u = sqrt(-2.0 * log(sdrand()));
    v = 2.0 * M_PI * sdrand();
    z[i] = u * sin(v);
    z[i + 1] = u * cos(v);
  }
  if (fmod(n, 2) < 0.5) {
    return;
  } else {
    u = sqrt(-2.0 * log(sdrand()));
    v = 2.0 * M_PI * sdrand();
    z[n - 1] = u * sin(v);
  }
  return;
}

void rt(double *z, int n, int dof) {
  /* Simulates n random t variable with dof degrees of freedom
     by simulating standard normals and chi-squared random variables.
     Chi-squared rvs simulated by rgamma function that simulates random gamma
     variables (see gammafns file for details).
     If dof is 0 return n gaussian random variables instead.
   */

  gauss(z, n);
  if (dof > 0) {
    double s = 0.5 * dof;
    double denom = sqrt(rgamma(s) / s);
    for (int j1 = 0; j1 < n; j1++) {
      z[j1] /= denom;
    }
  }
  return;
}

void chol(int n, double **A) {
  /* Performs cholesky decompositon of A and returns result in the
     same matrix - adapted from PJG Fortran function*/

  for (int j1 = 0; j1 < n; j1++) {
    double sum = A[j1][j1];
    for (int j2 = 0; j2 < j1; j2++) {
      sum -= pow(A[j1][j2], 2);
    }
    A[j1][j1] = sqrt(sum);

    for (int j2 = j1 + 1; j2 < n; j2++) {
      sum = A[j2][j1];
      for (int j3 = 0; j3 < j1; j3++) {
        sum -= A[j2][j3] * A[j1][j3];
      }
      A[j2][j1] = sum / A[j1][j1];
    }
  }
}

void perm(double *work, int n) {
  /* Randomly permutes the n-dim work vector */

  for (int i = 0; i < (n - 1); i++) {
    int j = i + (int)((n - i) * sdrand());
    if (j != i) {
      double temp = work[j];
      work[j] = work[i];
      work[i] = temp;
    }
  }
  return;
}

double ltprob(int dof, double z) {
  /* Evaluates the log of p.d.f. of a t variable with dof degrees of freedom
     at point z */

  double constt =
      loggamma(0.5 * (dof + 1)) - loggamma(0.5 * dof) - 0.5 * log(dof * M_PI);
  double out = constt - 0.5 * (dof + 1) * log(1.0 + pow(z, 2.0) / dof);
  return out;
}

double lnormprob(int n, double *mu_k_l, double **B_k_l, double *datai) {
  /* Evaluates log of p.d.f. for a multivariate normal for model
     k, of dimension n, component l. The summary of means and
     sqrt of cov matrices (for all models and all component)
     are supplied in mu and B */

  double work[n];

  for (int i = 0; i < n; i++) {
    work[i] = datai[i] - mu_k_l[i];
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      (work[i]) -= B_k_l[i][j] * work[j];
    }
    (work[i]) /= B_k_l[i][i];
  }
  double out = 0.0;
  for (int i = 0; i < n; i++) {
    out += (work[i] * work[i]);
  }
  out = -0.5 * out - (n / 2.0) * log(2.0 * M_PI) - log(det(n, B_k_l));
  return out;
}

double det(int n, double **B_k_l) {

  /* Evaluates the determinant of a matrix in B corresponding to model k,
     component l. */
  double out = 1.0;
  for (int i = 0; i < n; i++) {
    out *= B_k_l[i][i];
  }
  return out;
}
