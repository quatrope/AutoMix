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

  // --- Section 3 - Initial File handling ---------------------
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_log.data", fname);
  FILE *fpl = fopen(datafname, "w");
  sprintf(datafname, "%s_pk.data", fname);
  FILE *fpp = fopen(datafname, "w");
  sprintf(datafname, "%s_ac.data", fname);
  FILE *fpac = fopen(datafname, "w");
  sprintf(datafname, "%s_adapt.data", fname);
  FILE *fpad = fopen(datafname, "w");
  sprintf(datafname, "%s_cf.data", fname);
  FILE *fpcf = fopen(datafname, "w");

  // Print user options to log file

  fprintf(fpl, "seed: %ld\n", seed);
  fprintf(fpl, "m: %d\n", mode);
  fprintf(fpl, "a: %d\n", adapt);
  fprintf(fpl, "p: %d\n", doperm);
  fprintf(fpl, "n: %d\n", nsweep2);
  fprintf(fpl, "N: %d\n", nsweep);

  // --- Section 4.0 - Read in key variables from user functions -

  // nmodels is the number of models
  int nmodels = get_nmodels();
  if (nmodels > NMODELS_MAX) {
    printf("\nError:kmax too large \n");
    return 0;
  } else if (nmodels < 0) {
    printf("\nError:negative kmax \n");
    return 0;
  }

  int *model_dims = (int *)malloc(nmodels * sizeof(int));
  int *Lk = (int *)malloc(nmodels * sizeof(int));
  load_model_dims(nmodels, model_dims);

  double **lambda = (double **)malloc(nmodels * sizeof(double *));
  double ***mu = (double ***)malloc(nmodels * sizeof(double **));
  double ****B = (double ****)malloc(nmodels * sizeof(double ***));
  double **sig = (double **)malloc(nmodels * sizeof(double *));
  for (int k = 0; k < nmodels; k++) {
    int mdim = model_dims[k];
    lambda[k] = (double *)malloc(Lkmaxmax * sizeof(double));
    mu[k] = (double **)malloc(Lkmaxmax * sizeof(double *));
    B[k] = (double ***)malloc(Lkmaxmax * sizeof(double **));
    sig[k] = (double *)malloc(mdim * sizeof(double));
    for (int i = 0; i < Lkmaxmax; i++) {
      mu[k][i] = (double *)malloc(mdim * sizeof(double));
      B[k][i] = (double **)malloc(mdim * sizeof(double *));
      for (int j = 0; j < mdim; j++) {
        B[k][i][j] = (double *)malloc(mdim * sizeof(double));
      }
    }
  }

  // --- Section 5.1 - Read in mixture parameters if mode 1 (m=1) ---

  if (mode == 1) {
    // Read AutoMix parameters from file if mode = 1
    int ok =
        read_mixture_params(fname, nmodels, model_dims, sig, Lk, lambda, mu, B);
    if (ok == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
  } else {
    // --- Section 5.2 - Within-model runs if mixture parameters unavailable -
    for (int model_k = 0; model_k < nmodels; model_k++) {
      int mdim = model_dims[model_k];
      int lendata = 1000 * mdim;
      double **data = (double **)malloc(lendata * sizeof(*data));
      data[0] = (double *)malloc(lendata * mdim * sizeof(**data));
      for (int i = 1; i < lendata; i++) {
        data[i] = data[i - 1] + mdim;
      }
      // --- Section 5.2.1 - RWM Within Model (Stage 1) -------
      rwn_within_model(model_k, model_dims, nsweep2, fpl, fpcf, fpad, sig, dof,
                       data);

      printf("\nMixture Fitting: Model %d", model_k + 1);
      if (mode == 0) {
        // --- Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)-
        // Mixture fitting done component wise EM algorithm described in
        // Figueiredo and Jain, 2002 (see thesis for full reference)
        fit_mixture_from_samples(model_dims[model_k], data, lendata,
                                 mu[model_k], B[model_k], lambda[model_k], fpcf,
                                 Lk + model_k);
      }
      if (mode == 2) {
        //--- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --
        fit_autorj(lambda[model_k], Lk + model_k, model_dims[model_k],
                   mu[model_k], B[model_k], data, lendata);
      }
      free(data[0]);
      free(data);
    }
  }

  // Print mixture parameters to file (log and mix files)
  sprintf(datafname, "%s_mix.data", fname);
  FILE *fpmix = fopen(datafname, "w");
  fprintf(fpmix, "%d\n", nmodels);
  for (int k1 = 0; k1 < nmodels; k1++) {
    fprintf(fpmix, "%d\n", model_dims[k1]);
  }

  for (int k1 = 0; k1 < nmodels; k1++) {
    fprintf(fpl, "\nModel:%d\n", k1 + 1);
    int Lkk = Lk[k1];
    int mdim = model_dims[k1];
    fprintf(fpl, "\nARW params:\n");
    for (int j1 = 0; j1 < mdim; j1++) {
      fprintf(fpl, "%lf ", sig[k1][j1]);
    }
    fprintf(fpl, "\n");
    fprintf(fpl, "\nLkk:%d\n", Lkk);
    for (int j1 = 0; j1 < mdim; j1++) {
      fprintf(fpmix, "%lf\n", sig[k1][j1]);
    }
    fprintf(fpmix, "%d\n", Lkk);
    for (int l1 = 0; l1 < Lkk; l1++) {
      fprintf(fpl, "\nComponent:%d\n", l1 + 1);
      fprintf(fpl, "lambda:%lf\n", lambda[k1][l1]);
      fprintf(fpmix, "%lf\n", lambda[k1][l1]);
      fprintf(fpl, "mu:\n");
      for (int j1 = 0; j1 < mdim; j1++) {
        fprintf(fpl, "%lf ", mu[k1][l1][j1]);
        fprintf(fpmix, "%lf\n", mu[k1][l1][j1]);
      }
      fprintf(fpl, "\nB:\n");
      for (int j1 = 0; j1 < mdim; j1++) {
        for (int j2 = 0; j2 <= j1; j2++) {
          fprintf(fpl, "%lf ", B[k1][l1][j1][j2]);
          fprintf(fpmix, "%lf\n", B[k1][l1][j1][j2]);
        }
        fprintf(fpl, "\n");
      }
    }
  }
  fflush(NULL);

  // --Section 6 - Secondary file handling -------------

  sprintf(datafname, "%s_k.data", fname);
  FILE *fpk = fopen(datafname, "w");
  sprintf(datafname, "%s_lp.data", fname);
  FILE *fplp = fopen(datafname, "w");

  FILE *fpt[NMODELS_MAX];
  for (int k1 = 0; k1 < nmodels; k1++) {
    sprintf(datafname, "%s_theta%d.data", fname, k1 + 1);
    fpt[k1] = fopen(datafname, "w");
  }
  free(datafname);

  // --Section 7 - Final initialisation of variables ----

  int naccrwmb = 0;
  int ntryrwmb = 0;
  int naccrwms = 0;
  int ntryrwms = 0;
  int nacctd = 0;
  int ntrytd = 0;

  // Initialization of the MC Markov Chain parameters
  chainState ch;
  initChain(&ch, nmodels, model_dims, Lk, adapt);
  double llh;

  // propk and detB are auxiliary arrays
  double *propk = (double *)malloc(nmodels * sizeof(double));
  double **detB = (double **)malloc(nmodels * sizeof(double *));
  for (int k = 0; k < nmodels; k++) {
    detB[k] = (double *)malloc(Lkmaxmax * sizeof(double));
  }
  for (int k1 = 0; k1 < nmodels; k1++) {
    int Lkk = Lk[k1];
    int mdim = model_dims[k1];
    for (int l1 = 0; l1 < Lkk; l1++) {
      detB[k1][l1] = det(mdim, l1, B[k1]);
    }
  }
  for (int k1 = 0; k1 < nmodels; k1++) {
    if (k1 == ch.current_model_k) {
      propk[k1] = 1.0;
    } else {
      propk[k1] = 0.0;
    }
  }

  int nburn = max(10000, (int)(nsweep / 10));

  int nsokal = 1;
  int nkeep = nsweep / (2 * nsokal);
  nkeep = (int)pow(2.0, min(15, (int)(log(nkeep) / log(2.0) + 0.001)));
  int keep = nburn + (nsweep - nkeep * nsokal);
  double *xr = (double *)malloc(nkeep * sizeof(double));
  int *ksummary = (int *)calloc(nmodels, sizeof(int));

  // -----Start of main loop ----------------
  // Burn some samples first
  printf("\nBurning in");
  ch.isBurning = 1;
  for (int sweep = 1; sweep <= nburn; sweep++) {
    // Every 10 sweeps to block RWM
    ch.doBlockRWM = (sweep % 10 == 0);
    ch.gamma_sweep = pow(1.0 / (sweep + 1), (2.0 / 3.0));

    reversible_jump_move(&ch, B, Lk, detB, dof, lambda, &llh, nmodels,
                         model_dims, mu, &naccrwmb, &naccrwms, &nacctd,
                         &ntryrwmb, &ntryrwms, &ntrytd, propk, sig);
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

    reversible_jump_move(&ch, B, Lk, detB, dof, lambda, &llh, nmodels,
                         model_dims, mu, &naccrwmb, &naccrwms, &nacctd,
                         &ntryrwmb, &ntryrwms, &ntrytd, propk, sig);

    (ksummary[ch.current_model_k])++;
    // --- Section 10 - Write variables to files ---------

    fprintf(fpk, "%d\n", ch.current_model_k + 1);
    fprintf(fplp, "%lf %lf\n", ch.lp, llh);
    for (int k1 = 0; k1 < nmodels; k1++) {
      fprintf(fpp, "%lf ", ch.pk[k1]);
    }
    fprintf(fpp, "\n");
    for (int j1 = 0; j1 < ch.mdim; j1++) {
      fprintf(fpt[ch.current_model_k], "%lf ", ch.theta[j1]);
    }
    fprintf(fpt[ch.current_model_k], "\n");
    if (sweep > keep && ((sweep - keep) % nsokal == 0)) {
      xr[((sweep - keep) / nsokal) - 1] = ch.current_model_k;
    }
    if ((10 * (sweep - nburn)) % nsweep == 0) {
      printf("\nNo. of iterations remaining: %d", nsweep + nburn - sweep);
    }
    fflush(NULL);
  }
  printf("\n");
  freeChain(&ch);

  // --- Section 11 - Write log file ----------------------

  double var, tau;
  int m;
  sokal(nkeep, xr, &var, &tau, &m);
  fprintf(fpl, "\nAutocorrelation Time:\n");
  fprintf(fpl, "nkeep:%d, nsokal:%d, var:%lf, tau:%lf\n", nkeep, nsokal, var,
          tau);
  for (int i1 = 0; i1 < m; i1++) {
    fprintf(fpac, "%lf\n", xr[i1]);
  }

  fprintf(fpl, "\nPosterior Model Probabilities:\n");
  for (int k1 = 0; k1 < nmodels; k1++) {
    fprintf(fpl, "Model %d: %lf\n", k1 + 1,
            (double)ksummary[k1] / (double)nsweep);
  }

  fprintf(fpl, "\nAcceptance Rates:\n");
  fprintf(fpl, "Block RWM: %lf\n", (double)naccrwmb / (double)ntryrwmb);
  fprintf(fpl, "Single RWM: %lf\n", (double)naccrwms / (double)ntryrwms);
  fprintf(fpl, "Auto RJ: %lf\n", (double)nacctd / (double)ntrytd);

  clock_t endtime = clock();
  double timesecs = (endtime - starttime) / ((double)CLOCKS_PER_SEC);
  fprintf(fpl, "\nRun time:\n");
  fprintf(fpl, "Time: %lf\n", timesecs);

  return 0;
}

void initChain(chainState *ch, int nmodels, int *model_dims, int *Lk,
               int adapt) {
  ch->current_model_k = (int)floor(nmodels * sdrand());
  ch->mdim = model_dims[ch->current_model_k];
  int mdim_max = model_dims[0];
  for (int i = 1; i < nmodels; i++) {
    mdim_max = max(model_dims[i], mdim_max);
  }
  ch->theta = (double *)malloc(mdim_max * sizeof(double));
  ch->pk = (double *)malloc(nmodels * sizeof(double));
  get_rwm_init(ch->current_model_k, ch->mdim, ch->theta);
  ch->current_Lkk = Lk[ch->current_model_k];
  double llh;
  ch->lp = logpost(ch->current_model_k, ch->mdim, ch->theta, &llh);
  for (int i = 0; i < nmodels; i++) {
    ch->pk[i] = 1.0 / nmodels;
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

int read_mixture_params(char *fname, int nmodels, int *model_dims, double **sig,
                        int *Lk, double **lambda, double ***mu, double ****B) {
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
  if (k1 != nmodels) {
    printf("\nFile nmodels contradicts getkmax function:");
    return EXIT_FAILURE;
  }
  for (int k = 0; k < nmodels; k++) {
    int mdim;
    if (fscanf(fpmix, "%d", &mdim) == EOF) {
      printf("\nEnd of file encountered before parameters read:");
      return EXIT_FAILURE;
    }
    if (mdim != model_dims[k]) {
      printf("\nFile kmax contradicts getnk function:");
      return EXIT_FAILURE;
    }
  }
  for (int k = 0; k < nmodels; k++) {
    int mdim = model_dims[k];
    for (int l = 0; l < mdim; l++) {
      if (fscanf(fpmix, "%lf", &(sig[k][l])) == EOF) {
        printf("\nEnd of file encountered before parameters read:");
        return EXIT_FAILURE;
      }
    }
    if (fscanf(fpmix, "%d", &(Lk[k])) == EOF) {
      printf("\nEnd of file encountered before parameters read:");
      return EXIT_FAILURE;
    }
    int Lkk = Lk[k];
    for (int l = 0; l < Lkk; l++) {
      if (fscanf(fpmix, "%lf", &(lambda[k][l])) == EOF) {
        printf("\nEnd of file encountered before parameters read:");
        return EXIT_FAILURE;
      }
      for (int i = 0; i < mdim; i++) {
        if (fscanf(fpmix, "%lf", &(mu[k][l][i])) == EOF) {
          printf("\nEnd of file encountered before parameters read:");
          return EXIT_FAILURE;
        }
      }
      for (int i = 0; i < mdim; i++) {
        for (int j = 0; j <= i; j++) {
          if (fscanf(fpmix, "%lf", &(B[k][l][i][j])) == EOF) {
            printf("\nEnd of file encountered before parameters read:");
            return EXIT_FAILURE;
          }
        }
      }
    }
    double sumlambda = 0.0;
    for (int l = 0; l < Lkk; l++) {
      sumlambda += lambda[k][l];
    }
    double sumlambda_tol = 1E-5;
    if (fabs(sumlambda - 1.0) > sumlambda_tol) {
      printf("\nComponents weights read do not sum to one for k=%d:", k);
      return EXIT_FAILURE;
    }
    if (sumlambda != 1.0) {
      for (int l = 0; l < Lkk; l++) {
        lambda[k][l] /= sumlambda;
      }
    }
  }
  fclose(fpmix);
  return EXIT_SUCCESS;
}

void rwn_within_model(int model_k, int *model_dims, int nsweep2, FILE *fpl,
                      FILE *fpcf, FILE *fpad, double **sig, int dof,
                      double **data) {
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

  printf("\nRWM for Model %d", model_k + 1);
  fprintf(fpl, "\nRWM for Model %d", model_k + 1);
  fprintf(fpcf, "RWM for Model %d\n", model_k + 1);
  fprintf(fpad, "RWM for Model %d\n", model_k + 1);
  fflush(NULL);
  get_rwm_init(model_k, mdim, rwm);
  for (int j1 = 0; j1 < mdim; j1++) {
    rwmn[j1] = rwm[j1];
    sig[model_k][j1] = 10.0;
    nacc[j1] = 0;
    ntry[j1] = 0;
  }
  double llh, llhn, lpn;
  double lp = logpost(model_k, mdim, rwm, &llh);

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
        rwmn[i] = rwm[i] + sig[model_k][i] * Znkk[i];
      }
      lpn = logpost(model_k, mdim, rwmn, &llhn);
      if (sdrand() < exp(max(-30.0, min(0.0, lpn - lp)))) {
        for (int i = 0; i < mdim; i++) {
          rwm[i] = rwmn[i];
        }
        lp = lpn;
        llh = llhn; // This may not ever be used!
      }
    } else {
      double gamma = 10.0 * pow(1.0 / (sweep + 1), 2.0 / 3.0);
      for (int i = 0; i < mdim; i++) {
        rwmn[i] = rwm[i];
      }
      for (int i = 0; i < mdim; i++) {
        double Z;
        rt(&Z, 1, dof);
        rwmn[i] = rwm[i] + sig[model_k][i] * Z;
        lpn = logpost(model_k, mdim, rwmn, &llhn);
        double accept = min(1, exp(max(-30.0, min(0.0, lpn - lp))));
        if (sdrand() < accept) {
          (nacc[i])++;
          (ntry[i])++;
          rwm[i] = rwmn[i];
          lp = lpn;
          llh = llhn;
          sig[model_k][i] = max(0, sig[model_k][i] - gamma * (alphastar - 1));
        } else {
          (ntry[i])++;
          rwmn[i] = rwm[i];
          sig[model_k][i] = max(0, sig[model_k][i] - gamma * (alphastar));
        }
      }
    }
    if (remain < (10000 * mdim) && fmod(remain, 10.0) < 0.05) {
      for (int i = 0; i < mdim; i++) {
        data[i2][i] = rwm[i];
      }
      i2++;
    }
    if (fmod(sweep, 100.0) < 0.05) {
      for (int i = 0; i < mdim; i++) {
        fprintf(fpad, "%lf %lf ", sig[model_k][i],
                (double)nacc[i] / (double)ntry[i]);
      }
      fprintf(fpad, "\n");
    }
  }
  free(rwm);
  free(rwmn);
  free(nacc);
  free(ntry);
  free(Znkk);
}

void fit_mixture_from_samples(int mdim, double **data, int lendata,
                              double **mu_k, double ***B_k, double *lambda_k,
                              FILE *fpcf, int *Lk_k) {
  // --- Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)-
  // Mixture fitting done component wise EM algorithm described in
  // Figueiredo and Jain, 2002 (see thesis for full reference)
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
      lpdatagivenl[i][j] = lnormprob(mdim, j, mu_k, B_k, data[i]);
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
          lpdatagivenl[i1][l1] = lnormprob(mdim, l1, mu_k, B_k, data[i1]);
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
  *Lk_k = Lkkmin;
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

void fit_autorj(double *lambda_k, int *Lk_k, int mdim, double **mu_k,
                double ***B_k, double **data, int lendata) {
  // --- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --
  *Lk_k = 1;
  lambda_k[0] = 1.0;
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

void reversible_jump_move(chainState *ch, double ****B, int *Lk, double **detB,
                          int dof, double **lambda, double *llh, int nmodels,
                          int *model_dims, double ***mu, int *naccrwmb,
                          int *naccrwms, int *nacctd, int *ntryrwmb,
                          int *ntryrwms, int *ntrytd, double *propk,
                          double **sig) {
  int Lkmax = Lk[0];
  for (int k1 = 1; k1 < nmodels; k1++) {
    Lkmax = max(Lkmax, Lk[k1]);
  }
  int mdim_max = model_dims[0];
  for (int i = 1; i < nmodels; i++) {
    mdim_max = max(model_dims[i], mdim_max);
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
    (*ntryrwmb)++;
    rt(Znkk, ch->mdim, dof);
    for (int i = 0; i < ch->mdim; i++) {
      thetan[i] = theta[i] + sig[ch->current_model_k][i] * Znkk[i];
    }
    double llhn;
    double lpn = logpost(ch->current_model_k, ch->mdim, thetan, &llhn);
    if (sdrand() < exp(max(-30.0, min(0.0, lpn - ch->lp)))) {
      (*naccrwmb)++;
      memcpy(theta, thetan, ch->mdim * sizeof(*thetan));
      ch->lp = lpn;
      *llh = llhn;
    }
  } else {
    // else do component-wise RWM
    memcpy(thetan, theta, (ch->mdim) * sizeof(*thetan));
    for (int j1 = 0; j1 < ch->mdim; j1++) {
      (*ntryrwms)++;
      double Z;
      rt(&Z, 1, dof);
      thetan[j1] = theta[j1] + sig[ch->current_model_k][j1] * Z;
      double llhn;
      double lpn = logpost(ch->current_model_k, ch->mdim, thetan, &llhn);
      if (sdrand() < exp(max(-30.0, min(0.0, lpn - ch->lp)))) {
        (*naccrwms)++;
        theta[j1] = thetan[j1];
        ch->lp = lpn;
        *llh = llhn;
      } else {
        thetan[j1] = theta[j1];
      }
    }
  }

  // --- Section 9 Reversible Jump Moves -------------------

  // --Section 9.1 - Allocate current position to a component --
  int l = 0;
  (*ntrytd)++;
  int ln = 0;
  if (ch->current_Lkk > 1) {
    double sum = 0.0;
    for (int i = 0; i < ch->current_Lkk; i++) {
      palloc[i] = log(lambda[ch->current_model_k][i]) +
                  lnormprob(ch->mdim, i, mu[ch->current_model_k],
                            B[ch->current_model_k], theta);
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
    work[i] = theta[i] - mu[ch->current_model_k][l][i];
  }
  for (int i = 0; i < ch->mdim; i++) {
    for (int j = 0; j < i; j++) {
      work[i] = work[i] - B[ch->current_model_k][l][i][j] * work[j];
    }
    work[i] = work[i] / B[ch->current_model_k][l][i][i];
  }

  // --Section 9.3 - Choose proposed new model and component ----
  int kn = 0;
  double gamma = 0.0;
  double logratio;
  if (nmodels == 1) {
    kn = ch->current_model_k;
    logratio = 0.0;
  } else {
    gamma = ch->gamma_sweep;
    double u = sdrand();
    double thresh = 0.0;
    for (int k1 = 0; k1 < nmodels; k1++) {
      thresh += ch->pk[k1];
      if (u < thresh) {
        kn = k1;
        break;
      }
    }
    logratio = log(ch->pk[ch->current_model_k]) - log(ch->pk[kn]);
  }

  int mdim_kn = model_dims[kn];
  int Lkkn = Lk[kn];

  double u = sdrand();
  double thresh = 0.0;
  for (int l1 = 0; l1 < Lkkn; l1++) {
    thresh += lambda[kn][l1];
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
    thetan[j1] = mu[kn][ln][j1];
    for (int j2 = 0; j2 <= j1; j2++) {
      thetan[j1] += B[kn][ln][j1][j2] * work[j2];
    }
  }

  // --Section 9.5 - Work out probability of allocating to component
  // for acceptance ratio (reverse move) ---------

  if (Lkkn > 1) {
    double sum = 0.0;
    for (int l1 = 0; l1 < Lkkn; l1++) {
      pallocn[l1] =
          log(lambda[kn][l1]) + lnormprob(mdim_kn, l1, mu[kn], B[kn], thetan);
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
  double llhn;
  double lpn = logpost(kn, mdim_kn, thetan, &llhn);

  logratio += (lpn - ch->lp);
  logratio += (log(pallocn[ln]) - log(palloc[l]));
  logratio += (log(lambda[ch->current_model_k][l]) - log(lambda[kn][ln]));
  logratio += (log(detB[kn][ln]) - log(detB[ch->current_model_k][l]));

  if (sdrand() < exp(max(-30.0, min(0.0, logratio)))) {
    for (int j1 = 0; j1 < mdim_kn; j1++) {
      theta[j1] = thetan[j1];
    }
    ch->lp = lpn;
    *llh = llhn;
    ch->current_model_k = kn;
    ch->mdim = mdim_kn;
    ch->current_Lkk = Lkkn;
    (*nacctd)++;
  }

  if (ch->doAdapt && !ch->isBurning) {
    for (int k1 = 0; k1 < nmodels; k1++) {
      if (k1 == ch->current_model_k) {
        propk[k1] = 1.0;
      } else {
        propk[k1] = 0.0;
      }
      ch->pk[k1] += (gamma * (propk[k1] - ch->pk[k1]));
    }
    for (int k1 = 0; k1 < nmodels; k1++) {
      if (ch->pk[k1] < ch->pkllim) {
        ch->reinit = 1;
        break;
      }
    }
    if (ch->reinit == 1) {
      ch->reinit = 0;
      (ch->nreinit)++;
      ch->pkllim = 1.0 / (10.0 * ch->nreinit);
      for (int k1 = 0; k1 < nmodels; k1++) {
        ch->pk[k1] = 1.0 / nmodels;
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
