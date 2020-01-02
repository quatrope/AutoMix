// AutoMix - By David Hastie.
// See automix.h for full license and credits.

#define VERSION "1.2"

#include "automix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

void usage(char *invocation);
void write_stats_to_file(char *fname, unsigned long seed, int mode, int adapt,
                         int doperm, int nsweep2, int nsweep, proposalDist jd,
                         double **sig, int nkeep, int nsokal, runStats st,
                         double timesecs);
void write_mix_to_file(char *fname, proposalDist jd, double **sig);
void write_log_to_file(char *fname, unsigned long seed, int mode, int adapt,
                       int doperm, int nsweep2, int nsweep, proposalDist jd,
                       double **sig, int nkeep, int nsokal, runStats st,
                       double timesecs);
void write_ac_to_file(char *fname, int m, double *xr);
void write_theta_to_file(char *fname, int current_model_k, int mdim,
                         double *theta);
void write_pk_to_file(char *fname, int nsweep, int nmodels,
                      double **pk_summary);
void write_k_to_file(char *fname, int nsweep, int *k_which_summary);
void write_lp_to_file(char *fname, int nsweep, double **logp_summary);

int main(int argc, char *argv[]) {

  clock_t starttime = clock();

  // Default values for command line arguments
  // mode ~ m ~ 0 if mixture fitting, 1 if user supplied mixture params, 2 if
  // AutoRJ
  int mode = 0;
  // nsweep ~ N ~ no. of reversible jump sweeps in stage 3
  int nsweep = 1E5;
  // nsweep2 ~ n ~ max(n,10000*mdim,100000) sweeps in within-model RWM in stage
  // 1
  int nsweep2 = 1E5;
  // adapt ~ a ~ 1 if RJ adaptation done in stage 3, 0 otherwise
  int adapt = 1;
  // doperm ~ p ~ 1 if random permutation done in stage 3 RJ move, 0 otherwise
  int doperm = 1;
  // seed ~ s ~ random no. seed, 0 uses clock seed
  unsigned long seed = 0;
  // dof ~ t ~ 0 if Normal random variables used in RWM and RJ moves,
  // otherwise specify integer degrees of freedom of student t variables
  int dof = 0;
  // fname ~ f ~ filename base (default = "output")
  char *const fname_default = "output";
  char *fname = fname_default;

  // Override defaults if user supplies command line options
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-f")) {
      fname = argv[++i];
      continue;
    } else if (!strcmp(argv[i], "-N")) {
      nsweep = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-n")) {
      nsweep2 = atoi(argv[++i]);
      nsweep2 = max(nsweep2, 1E5);
      continue;
    } else if (!strcmp(argv[i], "-s")) {
      seed = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-p")) {
      doperm = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-m")) {
      mode = atoi(argv[++i]);
      if (mode > 2 || mode < 0) {
        printf("Error: Invalid mode entered. Mode must be {0, 1, 2}.\n");
        usage(argv[0]);
        return EXIT_FAILURE;
      }
      continue;
    } else if (!strcmp(argv[i], "-a")) {
      adapt = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-t")) {
      dof = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      usage(argv[0]);
      return EXIT_SUCCESS;
    } else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--version")) {
      printf("Version %s\n", VERSION);
      return EXIT_SUCCESS;
    } else {
      printf("Unrecognized argument: %s\n", argv[i]);
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }
  sdrni(&seed);
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));

  // Set running sweep variables
  int nburn = max(10000, (int)(nsweep / 10));
  int nsokal = 1;
  int nkeep = nsweep / (2 * nsokal);
  nkeep = (int)pow(2.0, min(15, (int)(log(nkeep) / log(2.0) + 0.001)));
  int keep = nburn + (nsweep - nkeep * nsokal);

  // Initialize the Proposal (jumping) Distribution
  proposalDist jd;
  initJD(&jd);
  // Struct to hold run statistic variables
  runStats st;
  initializeRunStats(&st, nkeep, nsweep, jd);

  double **sig = (double **)malloc(jd.nmodels * sizeof(double *));
  for (int k = 0; k < jd.nmodels; k++) {
    int mdim = jd.model_dims[k];
    sig[k] = (double *)malloc(mdim * sizeof(double));
  }

  // --- Section 5.1 - Read in mixture parameters if mode 1 (m=1) ---
  sprintf(datafname, "%s_cf.data", fname);
  FILE *fp_cf = fopen(datafname, "w");
  sprintf(datafname, "%s_adapt.data", fname);
  FILE *fp_adapt = fopen(datafname, "w");
  if (mode == 1) {
    // Read AutoMix parameters from file if mode = 1
    int ok = read_mixture_params(fname, jd, sig);
    if (ok == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
  } else {
    // --- Section 5.2 - Within-model runs if mixture parameters unavailable -
    for (int model_k = 0; model_k < jd.nmodels; model_k++) {
      int mdim = jd.model_dims[model_k];
      int lendata = 1000 * mdim;
      double **data = (double **)malloc(lendata * sizeof(*data));
      data[0] = (double *)malloc(lendata * mdim * sizeof(**data));
      for (int i = 1; i < lendata; i++) {
        data[i] = data[i - 1] + mdim;
      }
      int nsweepr = max(nsweep2, 10000 * mdim);

      int rwm_summary_len = max(1, (nsweepr + nburn) / 100);
      st.sig_k_rwm_summary =
          (double **)malloc(rwm_summary_len * sizeof(double *));
      st.sig_k_rwm_summary[0] =
          (double *)malloc(rwm_summary_len * mdim * sizeof(double));
      for (int i = 1; i < rwm_summary_len; i++) {
        st.sig_k_rwm_summary[i] = st.sig_k_rwm_summary[i - 1] + mdim;
      }
      st.nacc_ntry_rwm = (double **)malloc(rwm_summary_len * sizeof(double *));
      st.nacc_ntry_rwm[0] =
          (double *)malloc(rwm_summary_len * mdim * sizeof(double));
      for (int i = 1; i < rwm_summary_len; i++) {
        st.nacc_ntry_rwm[i] = st.nacc_ntry_rwm[i - 1] + mdim;
      }

      // --- Section 5.2.1 - RWM Within Model (Stage 1) -------
      fprintf(fp_cf, "RWM for Model %d\n", model_k + 1);
      fprintf(fp_adapt, "RWM for Model %d\n", model_k + 1);
      rwm_within_model(model_k, jd.model_dims, nsweepr, st, sig[model_k], dof,
                       data);
      for (int j = 0; j < rwm_summary_len; j++) {
        for (int i = 0; i < mdim; i++) {
          fprintf(fp_adapt, "%lf %lf ", st.sig_k_rwm_summary[j][i],
                  st.nacc_ntry_rwm[j][i]);
        }
        fprintf(fp_adapt, "\n");
      }
      free(st.sig_k_rwm_summary[0]);
      free(st.sig_k_rwm_summary);
      free(st.nacc_ntry_rwm[0]);
      free(st.nacc_ntry_rwm);

      printf("\nMixture Fitting: Model %d", model_k + 1);
      if (mode == 0) {
        // --- Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)-
        // Mixture fitting done component wise EM algorithm described in
        // Figueiredo and Jain, 2002 (see thesis for full reference)
        fit_mixture_from_samples(model_k, jd, data, lendata, fp_cf);
      }
      if (mode == 2) {
        //--- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --
        fit_autorj(model_k, jd, data, lendata);
      }
      free(data[0]);
      free(data);
    }
  }
  fclose(fp_adapt);
  fclose(fp_cf);
  free(datafname);

  // Print mixture parameters to file
  write_mix_to_file(fname, jd, sig);

  // Initialization of the MC Markov Chain parameters
  chainState ch;
  initChain(&ch, jd, adapt);

  // -----Start of main loop ----------------
  // Burn some samples first
  printf("\nBurning in");
  ch.isBurning = 1;
  for (int sweep = 1; sweep <= nburn; sweep++) {
    // Every 10 sweeps to block RWM
    ch.doBlockRWM = (sweep % 10 == 0);
    ch.gamma_sweep = pow(1.0 / (sweep + 1), (2.0 / 3.0));

    reversible_jump_move(&ch, jd, dof, &st, sig);
    if ((10 * sweep) % nburn == 0) {
      printf(" .");
      fflush(NULL);
    }
  }
  printf("\n");

  // Start here main sample
  printf("Start of main sample:");
  ch.isBurning = 0;
  for (int sweep = nburn + 1; sweep <= (nburn + nsweep); sweep++) {
    // Every 10 sweeps to block RWM
    ch.doBlockRWM = (sweep % 10 == 0);
    ch.gamma_sweep = pow(1.0 / (sweep + 1), (2.0 / 3.0));

    reversible_jump_move(&ch, jd, dof, &st, sig);

    (st.ksummary[ch.current_model_k])++;
    st.k_which_summary[sweep - nburn - 1] = ch.current_model_k + 1;
    st.logp_summary[sweep - nburn - 1][0] = ch.log_posterior;
    st.logp_summary[sweep - nburn - 1][1] = ch.log_likelihood;
    for (int k1 = 0; k1 < jd.nmodels; k1++) {
      st.pk_summary[sweep - nburn - 1][k1] = ch.pk[k1];
    }
    write_theta_to_file(fname, ch.current_model_k, ch.mdim, ch.theta);

    if (sweep > keep && ((sweep - keep) % nsokal == 0)) {
      st.xr[((sweep - keep) / nsokal) - 1] = ch.current_model_k;
    }
    if ((10 * (sweep - nburn)) % nsweep == 0) {
      printf("\nNo. of iterations remaining: %d", nsweep + nburn - sweep);
    }
    fflush(NULL);
  }
  printf("\n");
  freeChain(&ch);

  clock_t endtime = clock();
  double timesecs = (endtime - starttime) / ((double)CLOCKS_PER_SEC);

  // --- Section 10 - Write statistics to files ---------
  write_stats_to_file(fname, seed, mode, adapt, doperm, nsweep2, nsweep, jd,
                      sig, nkeep, nsokal, st, timesecs);
  freeJD(jd);

  return EXIT_SUCCESS;
}

void write_lp_to_file(char *fname, int nsweep, double **logp_summary) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_lp.data", fname);
  FILE *fp_lp = fopen(datafname, "w");
  free(datafname);
  for (int i = 0; i < nsweep; i++) {
    fprintf(fp_lp, "%lf %lf\n", logp_summary[i][0], logp_summary[i][1]);
  }
  fclose(fp_lp);
}

void write_k_to_file(char *fname, int nsweep, int *k_which_summary) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_k.data", fname);
  FILE *fp_k = fopen(datafname, "w");
  free(datafname);
  for (int i = 0; i < nsweep; i++) {
    fprintf(fp_k, "%d\n", k_which_summary[i]);
  }
  fclose(fp_k);
}

void write_pk_to_file(char *fname, int nsweep, int nmodels,
                      double **pk_summary) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_pk.data", fname);
  FILE *fp_pk = fopen(datafname, "w");
  free(datafname);
  for (int i = 0; i < nsweep; i++) {
    for (int j = 0; j < nmodels; j++) {
      fprintf(fp_pk, "%lf ", pk_summary[i][j]);
    }
    fprintf(fp_pk, "\n");
  }
  fclose(fp_pk);
}

void write_theta_to_file(char *fname, int current_model_k, int mdim,
                         double *theta) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_theta%d.data", fname, current_model_k + 1);
  FILE *fp_theta = fopen(datafname, "w");
  free(datafname);
  for (int j1 = 0; j1 < mdim; j1++) {
    fprintf(fp_theta, "%lf ", theta[j1]);
  }
  fprintf(fp_theta, "\n");
  fclose(fp_theta);
}

void write_stats_to_file(char *fname, unsigned long seed, int mode, int adapt,
                         int doperm, int nsweep2, int nsweep, proposalDist jd,
                         double **sig, int nkeep, int nsokal, runStats st,
                         double timesecs) {
  write_pk_to_file(fname, nsweep, jd.nmodels, st.pk_summary);
  write_k_to_file(fname, nsweep, st.k_which_summary);
  write_lp_to_file(fname, nsweep, st.logp_summary);
  sokal(nkeep, st.xr, &(st.var), &(st.tau), &(st.m));
  write_log_to_file(fname, seed, mode, adapt, doperm, nsweep2, nsweep, jd, sig,
                    nkeep, nsokal, st, timesecs);
  write_ac_to_file(fname, st.m, st.xr);
}

void write_ac_to_file(char *fname, int m, double *xr) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_ac.data", fname);
  FILE *fp_ac = fopen(datafname, "w");
  free(datafname);
  for (int i1 = 0; i1 < m; i1++) {
    fprintf(fp_ac, "%lf\n", xr[i1]);
  }
  fclose(fp_ac);
}

void write_mix_to_file(char *fname, proposalDist jd, double **sig) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_mix.data", fname);
  FILE *fp_mix = fopen(datafname, "w");
  free(datafname);
  fprintf(fp_mix, "%d\n", jd.nmodels);
  for (int k1 = 0; k1 < jd.nmodels; k1++) {
    fprintf(fp_mix, "%d\n", jd.model_dims[k1]);
  }
  for (int k1 = 0; k1 < jd.nmodels; k1++) {
    int Lkk = jd.nMixComps[k1];
    int mdim = jd.model_dims[k1];
    for (int j1 = 0; j1 < mdim; j1++) {
      fprintf(fp_mix, "%lf\n", sig[k1][j1]);
    }
    fprintf(fp_mix, "%d\n", Lkk);
    for (int l1 = 0; l1 < Lkk; l1++) {
      fprintf(fp_mix, "%lf\n", jd.lambda[k1][l1]);
      for (int j1 = 0; j1 < mdim; j1++) {
        fprintf(fp_mix, "%lf\n", jd.mu[k1][l1][j1]);
      }
      for (int j1 = 0; j1 < mdim; j1++) {
        for (int j2 = 0; j2 <= j1; j2++) {
          fprintf(fp_mix, "%lf\n", jd.B[k1][l1][j1][j2]);
        }
      }
    }
  }
  fclose(fp_mix);
}

void write_log_to_file(char *fname, unsigned long seed, int mode, int adapt,
                       int doperm, int nsweep2, int nsweep, proposalDist jd,
                       double **sig, int nkeep, int nsokal, runStats st,
                       double timesecs) {

  // Print user options to log file
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_log.data", fname);
  FILE *fp_log = fopen(datafname, "w");
  free(datafname);
  fprintf(fp_log, "seed: %ld\n", seed);
  fprintf(fp_log, "m: %d\n", mode);
  fprintf(fp_log, "a: %d\n", adapt);
  fprintf(fp_log, "p: %d\n", doperm);
  fprintf(fp_log, "n: %d\n", nsweep2);
  fprintf(fp_log, "N: %d\n", nsweep);
  if (mode != 1) {
    for (int model_k = 0; model_k < jd.nmodels; model_k++) {
      fprintf(fp_log, "\nRWM for Model %d", model_k + 1);
    }
  }
  for (int k1 = 0; k1 < jd.nmodels; k1++) {
    fprintf(fp_log, "\nModel:%d\n", k1 + 1);
    int Lkk = jd.nMixComps[k1];
    int mdim = jd.model_dims[k1];
    fprintf(fp_log, "\nARW params:\n");
    for (int j1 = 0; j1 < mdim; j1++) {
      fprintf(fp_log, "%lf ", sig[k1][j1]);
    }
    fprintf(fp_log, "\n");
    fprintf(fp_log, "\nLkk:%d\n", Lkk);
    for (int l1 = 0; l1 < Lkk; l1++) {
      fprintf(fp_log, "\nComponent:%d\n", l1 + 1);
      fprintf(fp_log, "lambda:%lf\n", jd.lambda[k1][l1]);
      fprintf(fp_log, "mu:\n");
      for (int j1 = 0; j1 < mdim; j1++) {
        fprintf(fp_log, "%lf ", jd.mu[k1][l1][j1]);
      }
      fprintf(fp_log, "\nB:\n");
      for (int j1 = 0; j1 < mdim; j1++) {
        for (int j2 = 0; j2 <= j1; j2++) {
          fprintf(fp_log, "%lf ", jd.B[k1][l1][j1][j2]);
        }
        fprintf(fp_log, "\n");
      }
    }
  }
  fprintf(fp_log, "\nAutocorrelation Time:\n");
  fprintf(fp_log, "nkeep:%d, nsokal:%d, var:%lf, tau:%lf\n", nkeep, nsokal,
          st.var, st.tau);
  fprintf(fp_log, "\nPosterior Model Probabilities:\n");
  for (int i = 0; i < jd.nmodels; i++) {
    fprintf(fp_log, "Model %d: %lf\n", i + 1,
            (double)st.ksummary[i] / (double)nsweep);
  }
  fprintf(fp_log, "\nAcceptance Rates:\n");
  fprintf(fp_log, "Block RWM: %lf\n",
          (double)st.naccrwmb / (double)st.ntryrwmb);
  fprintf(fp_log, "Single RWM: %lf\n",
          (double)st.naccrwms / (double)st.ntryrwms);
  fprintf(fp_log, "Auto RJ: %lf\n", (double)st.nacctd / (double)st.ntrytd);
  fprintf(fp_log, "\nRun time:\n");
  fprintf(fp_log, "Time: %lf\n", timesecs);
  fclose(fp_log);
}

void initializeRunStats(runStats *st, int nkeep, int nsweep, proposalDist jd) {
  st->naccrwmb = 0;
  st->ntryrwmb = 0;
  st->naccrwms = 0;
  st->ntryrwms = 0;
  st->nacctd = 0;
  st->ntrytd = 0;
  st->xr = (double *)malloc(nkeep * sizeof(double));
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
}

void freeRunStats(runStats st) {
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

void initChain(chainState *ch, proposalDist jd, int adapt) {
  ch->current_model_k = (int)floor(jd.nmodels * sdrand());
  ch->mdim = jd.model_dims[ch->current_model_k];
  int mdim_max = jd.model_dims[0];
  for (int i = 1; i < jd.nmodels; i++) {
    mdim_max = max(jd.model_dims[i], mdim_max);
  }
  ch->theta = (double *)malloc(mdim_max * sizeof(double));
  ch->pk = (double *)malloc(jd.nmodels * sizeof(double));
  get_rwm_init(ch->current_model_k, ch->mdim, ch->theta);
  ch->current_Lkk = jd.nMixComps[ch->current_model_k];
  logpost(ch->current_model_k, ch->mdim, ch->theta, &(ch->log_posterior),
          &(ch->log_likelihood));
  for (int i = 0; i < jd.nmodels; i++) {
    ch->pk[i] = 1.0 / jd.nmodels;
  }
  ch->nreinit = 1;
  ch->reinit = 0;
  ch->pkllim = 1.0 / 10.0;
  ch->doAdapt = adapt;
  ch->isInitialized = 1;
}

void freeChain(chainState *ch) {
  if (ch->theta != NULL) {
    free(ch->theta);
  }
  if (ch->pk != NULL) {
    free(ch->pk);
  }
  ch->isInitialized = 0;
}

int initJD(proposalDist *jd) {
  jd->nmodels = get_nmodels();
  if (jd->nmodels > NMODELS_MAX) {
    printf("\nError:kmax too large \n");
    return EXIT_FAILURE;
  } else if (jd->nmodels < 0) {
    printf("\nError:negative kmax \n");
    return EXIT_FAILURE;
  }
  allocJD(jd);
  jd->isInitialized = 1;
  return EXIT_SUCCESS;
}

int allocJD(proposalDist *jd) {
  // nmodels is the number of models
  int Lkmaxmax = NUM_MIX_COMPS_MAX;
  jd->model_dims = (int *)malloc(jd->nmodels * sizeof(int));
  if (jd->model_dims == NULL) {
    return EXIT_FAILURE;
  }
  load_model_dims(jd->nmodels, jd->model_dims);
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
    jd->lambda[k] = (double *)malloc(Lkmaxmax * sizeof(double));
    if (jd->lambda[k] == NULL) {
      return EXIT_FAILURE;
    }
    jd->mu[k] = (double **)malloc(Lkmaxmax * sizeof(double *));
    if (jd->mu[k] == NULL) {
      return EXIT_FAILURE;
    }
    jd->B[k] = (double ***)malloc(Lkmaxmax * sizeof(double **));
    if (jd->B[k] == NULL) {
      return EXIT_FAILURE;
    }
    for (int i = 0; i < Lkmaxmax; i++) {
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
  return EXIT_SUCCESS;
}

void freeJD(proposalDist jd) {
  int Lkmaxmax = NUM_MIX_COMPS_MAX;
  for (int k = 0; k < jd.nmodels; k++) {
    int mdim = jd.model_dims[k];
    for (int i = 0; i < Lkmaxmax; i++) {
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

int read_mixture_params(char *fname, proposalDist jd, double **sig) {
  // Check user has supplied mixture parameters if trying to use mode 1.
  // If not, abort and exit.
  char *datafname = (char *)malloc((strlen(fname) + 20) * sizeof(*datafname));
  sprintf(datafname, "%s_mix.data", fname);
  FILE *fpmix = fopen(datafname, "r");
  free(datafname);
  if (fpmix == NULL) {
    printf("\nProblem opening mixture file:");
    return EXIT_FAILURE;
  }
  int k1;
  if (fscanf(fpmix, "%d", &k1) == EOF) {
    printf("\nEnd of file encountered before parameters read:");
    return EXIT_FAILURE;
  }
  if (k1 != jd.nmodels) {
    printf("\nFile nmodels contradicts getkmax function:");
    return EXIT_FAILURE;
  }
  for (int k = 0; k < jd.nmodels; k++) {
    int mdim;
    if (fscanf(fpmix, "%d", &mdim) == EOF) {
      printf("\nEnd of file encountered before parameters read:");
      return EXIT_FAILURE;
    }
    if (mdim != jd.model_dims[k]) {
      printf("\nFile kmax contradicts getnk function:");
      return EXIT_FAILURE;
    }
  }
  for (int k = 0; k < jd.nmodels; k++) {
    int mdim = jd.model_dims[k];
    for (int l = 0; l < mdim; l++) {
      if (fscanf(fpmix, "%lf", &(sig[k][l])) == EOF) {
        printf("\nEnd of file encountered before parameters read:");
        return EXIT_FAILURE;
      }
    }
    if (fscanf(fpmix, "%d", &(jd.nMixComps[k])) == EOF) {
      printf("\nEnd of file encountered before parameters read:");
      return EXIT_FAILURE;
    }
    int Lkk = jd.nMixComps[k];
    for (int l = 0; l < Lkk; l++) {
      if (fscanf(fpmix, "%lf", &(jd.lambda[k][l])) == EOF) {
        printf("\nEnd of file encountered before parameters read:");
        return EXIT_FAILURE;
      }
      for (int i = 0; i < mdim; i++) {
        if (fscanf(fpmix, "%lf", &(jd.mu[k][l][i])) == EOF) {
          printf("\nEnd of file encountered before parameters read:");
          return EXIT_FAILURE;
        }
      }
      for (int i = 0; i < mdim; i++) {
        for (int j = 0; j <= i; j++) {
          if (fscanf(fpmix, "%lf", &(jd.B[k][l][i][j])) == EOF) {
            printf("\nEnd of file encountered before parameters read:");
            return EXIT_FAILURE;
          }
        }
      }
    }
    double sumlambda = 0.0;
    for (int l = 0; l < Lkk; l++) {
      sumlambda += jd.lambda[k][l];
    }
    double sumlambda_tol = 1E-5;
    if (fabs(sumlambda - 1.0) > sumlambda_tol) {
      printf("\nComponents weights read do not sum to one for k=%d:", k);
      return EXIT_FAILURE;
    }
    if (sumlambda != 1.0) {
      for (int l = 0; l < Lkk; l++) {
        jd.lambda[k][l] /= sumlambda;
      }
    }
  }
  fclose(fpmix);
  return EXIT_SUCCESS;
}

void rwm_within_model(int model_k, int *model_dims, int nsweepr, runStats st,
                      double *sig_k, int dof, double **data) {
  // --- Section 5.2.1 - RWM Within Model (Stage 1) -------

  int mdim = model_dims[model_k];
  int nburn = nsweepr / 10;
  double alphastar = 0.25;
  nsweepr += nburn;
  double *rwm = (double *)malloc(mdim * sizeof(double));
  double *rwmn = (double *)malloc(mdim * sizeof(double));
  int *nacc = (int *)malloc(mdim * sizeof(int));
  int *ntry = (int *)malloc(mdim * sizeof(int));
  double *Znkk = (double *)malloc(mdim * sizeof(double));

  int sig_k_rwm_n = 0; // This is used to load stats arrays
  printf("\nRWM for Model %d", model_k + 1);
  fflush(NULL);
  get_rwm_init(model_k, mdim, rwm);
  for (int j1 = 0; j1 < mdim; j1++) {
    rwmn[j1] = rwm[j1];
    sig_k[j1] = 10.0;
    nacc[j1] = 0;
    ntry[j1] = 0;
  }
  double lp, llh;
  logpost(model_k, mdim, rwm, &lp, &llh);

  int i2 = 0;
  int remain = nsweepr;
  for (int sweep = 1; sweep <= nsweepr; sweep++) {
    remain--;
    if ((sweep >= nburn) && ((10 * (sweep - nburn)) % (nsweepr - nburn) == 0)) {
      printf("\nNo. of iterations remaining: %d", remain);
      fflush(NULL);
    }
    double u = sdrand();
    if (sweep > nburn && u < 0.1) {
      rt(Znkk, mdim, dof);
      for (int i = 0; i < mdim; i++) {
        rwmn[i] = rwm[i] + sig_k[i] * Znkk[i];
      }
      double lpn;
      logpost(model_k, mdim, rwmn, &lpn, &llh);
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
        double lpn;
        logpost(model_k, mdim, rwmn, &lpn, &llh);
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
    if (remain < (10000 * mdim) && fmod(remain, 10.0) < 0.05) {
      for (int i = 0; i < mdim; i++) {
        data[i2][i] = rwm[i];
      }
      i2++;
    }
    if (sweep % 100 == 0) {
      for (int i = 0; i < mdim; i++) {
        st.sig_k_rwm_summary[sig_k_rwm_n][i] = sig_k[i];
        st.nacc_ntry_rwm[sig_k_rwm_n][i] = (double)nacc[i] / (double)ntry[i];
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

void fit_mixture_from_samples(int model_k, proposalDist jd, double **data,
                              int lendata, FILE *fpcf) {
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
    init[l1] = (int)floor(lendata * u);
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
    for (int i = 0; i < lendata; i++) {
      datasum += data[i][j];
      datasqsum += data[i][j] * data[i][j];
    }
    double len = (double)lendata;
    sigma += (datasqsum - datasum * datasum / len) / len;
  }
  sigma /= (10.0 * mdim);

  for (int i = 0; i < Lkk; i++) {
    for (int j = 0; j < mdim; j++) {
      mu_k[i][j] = data[init[i]][j];
      B_k[i][j][j] = sigma;
      for (int k = 0; k < j; k++) {
        B_k[i][j][k] = 0.0;
      }
    }
    chol(mdim, B_k[i]);
    lambda_k[i] = 1.0 / Lkk;
  }

  double **w = (double **)malloc(lendata * sizeof(double *));
  double *logw = (double *)malloc(Lkk * sizeof(double));
  double **lpdatagivenl = (double **)malloc(lendata * sizeof(double *));
  for (int i1 = 0; i1 < lendata; i1++) {
    w[i1] = (double *)malloc(Lkk * sizeof(double));
    lpdatagivenl[i1] = (double *)malloc(Lkk * sizeof(double));
  }

  for (int i = 0; i < lendata; i++) {
    double sum = 0.0;
    for (int j = 0; j < Lkk; j++) {
      lpdatagivenl[i][j] = lnormprob(mdim, mu_k[j], B_k[j], data[i]);
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
        for (int i1 = 0; i1 < lendata; i1++) {
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
          for (int i1 = 0; i1 < lendata; i1++) {
            mu_k[l1][j1] += data[i1][j1] * w[i1][l1];
          }
          mu_k[l1][j1] /= sumw[l1];

          for (int j2 = 0; j2 <= j1; j2++) {
            double temp = 0.0;
            for (int i1 = 0; i1 < lendata; i1++) {
              temp += (data[i1][j1] - mu_k[l1][j1]) *
                      (data[i1][j2] - mu_k[l1][j2]) * w[i1][l1];
            }
            B_k[l1][j1][j2] = temp / sumw[l1];
          }
        }

        chol(mdim, B_k[l1]);

        for (int i1 = 0; i1 < lendata; i1++) {
          lpdatagivenl[i1][l1] = lnormprob(mdim, mu_k[l1], B_k[l1], data[i1]);
        }
        l1++;

      } else {
        if (Lkk % 5 == 0) {
          printf("\n");
        }
        printf("%d(%d-n) ", Lkk, count);
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
            for (int i1 = 0; i1 < lendata; i1++) {
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
      for (int i1 = 0; i1 < lendata; i1++) {
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
      sum += log(lendata * lambda_k[l1] / 12.0);
    }
    double costfnnew = (nparams / 2.0) * sum +
                       (Lkk / 2.0) * log(lendata / 12.0) +
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
        if (Lkk % 5 == 0) {
          printf("\n");
        }
        printf("%d(%d-f) ", Lkk, count);
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
            for (int i1 = 0; i1 < lendata; i1++) {
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
        for (int i1 = 0; i1 < lendata; i1++) {
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
          sum += log(lendata * lambda_k[l1] / 12.0);
        }
        costfnnew = (nparams / 2.0) * sum + (Lkk / 2.0) * log(lendata / 12.0) +
                    Lkk * (nparams + 1) / 2.0 - lpn;
      }
    }
    if (count > 5000) {
      stop = 1;
    }
    costfn = costfnnew;
    fprintf(fpcf, "%d %lf %lf %d\n", Lkk, lpn, costfnnew, (natann + forceann));
    fflush(NULL);
  }

  for (int i = 0; i < lendata; i++) {
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

void fit_autorj(int model_k, proposalDist jd, double **data, int lendata) {
  // --- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --
  jd.nMixComps[model_k] = 1;
  jd.lambda[model_k][0] = 1.0;
  int mdim = jd.model_dims[model_k];
  double **mu_k = jd.mu[model_k];
  double ***B_k = jd.B[model_k];
  for (int j = 0; j < mdim; j++) {
    mu_k[0][j] = 0.0;
    for (int i = 0; i < lendata; i++) {
      mu_k[0][j] += data[i][j];
    }
    mu_k[0][j] /= ((double)lendata);
  }
  for (int l = 0; l < mdim; l++) {
    for (int j = 0; j <= l; j++) {
      B_k[0][l][j] = 0.0;
      for (int i = 0; i < lendata; i++) {
        B_k[0][l][j] += (data[i][l] - mu_k[0][l]) * (data[i][j] - mu_k[0][j]);
      }
      B_k[0][l][j] /= ((double)(lendata - 1));
    }
  }
  chol(mdim, B_k[0]);
}

void reversible_jump_move(chainState *ch, proposalDist jd, int dof,
                          runStats *st, double **sig) {
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
      thetan[i] = theta[i] + sig[ch->current_model_k][i] * Znkk[i];
    }
    double lpn, llhn;
    logpost(ch->current_model_k, ch->mdim, thetan, &lpn, &llhn);
    if (sdrand() < exp(max(-30.0, min(0.0, lpn - ch->log_posterior)))) {
      (st->naccrwmb)++;
      memcpy(theta, thetan, ch->mdim * sizeof(*thetan));
      ch->log_posterior = lpn;
      ch->log_likelihood = llhn;
    }
  } else {
    // else do component-wise RWM
    memcpy(thetan, theta, (ch->mdim) * sizeof(*thetan));
    for (int j1 = 0; j1 < ch->mdim; j1++) {
      (st->ntryrwms)++;
      double Z;
      rt(&Z, 1, dof);
      thetan[j1] = theta[j1] + sig[ch->current_model_k][j1] * Z;
      double lpn, llhn;
      logpost(ch->current_model_k, ch->mdim, thetan, &lpn, &llhn);
      if (sdrand() < exp(max(-30.0, min(0.0, lpn - ch->log_posterior)))) {
        (st->naccrwms)++;
        theta[j1] = thetan[j1];
        ch->log_posterior = lpn;
        ch->log_likelihood = llhn;
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
    gamma = ch->gamma_sweep;
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
    if (ch->doPerm) {
      perm(work, mdim_kn);
    }
  } else if (ch->mdim == mdim_kn) {
    if (ch->doPerm) {
      perm(work, ch->mdim);
    }
  } else {
    if (ch->doPerm) {
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
  double lpn, llhn;
  logpost(kn, mdim_kn, thetan, &lpn, &llhn);

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
    ch->log_likelihood = llhn;
    ch->current_model_k = kn;
    ch->mdim = mdim_kn;
    ch->current_Lkk = Lkkn;
    (st->nacctd)++;
  }

  if (ch->doAdapt && !ch->isBurning) {
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

void usage(char *invocation) {
  char *name = strrchr(invocation, '/');
  if (name == NULL) {
    name = invocation;
  } else {
    name += 1;
  }
  printf("%s version %s\n", name, VERSION);
  printf("Usage: %s [-m int] [-N int] [-n int] [-a bool] \
           [-p bool] [-s int] [-t int] [-f string] [-h, --help]\n",
         name);
  printf(
      "-m int: Specify mode. 0 if mixture fitting, 1 if user supplied mixture params, \
           2 if AutoRJ.\n");
  printf("-N int: Number of reversible jump sweeps in stage 3.\n");
  printf("-n int: max(n, 10000 * mdim, 100000) sweeps in within-model RWM in "
         "stage 1.\n");
  printf("-a bool: 1 if RJ adaptation done in stage 3, 0 otherwise.\n");
  printf("-p bool: 1 if random permutation done in stage 3 RJ move, 0 "
         "otherwise.\n");
  printf("-s int: random no. seed. 0 uses clock seed.\n");
  printf(
      "-t int: 0 if Normal random variables used in RWM and RJ moves, otherwise \
           specify integer degrees of freedom of student t variables.\n");
  printf("-f string: filename base.\n");
  printf("-h, --help: Print this help and exit.");
  printf(
      "\n(c) Original code by David Hastie. Modifications by Martin Beroiz.\n");
}
