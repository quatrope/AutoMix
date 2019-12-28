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

/* Standard library files */

#include "user.h"
#include "utils.h"
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

/* Global constants (please feel free to change as required)

   nkmaxmax = maximum dimension of any one model under consideration
   kmaxmax = maximum number of models
   Lkmaxmax = initial number of mixture components fitted in stage 2 of
              AutoMix algorithm */

#define nkmaxmax 20
#define kmaxmax 15
#define Lkmaxmax 30
#define tpi 6.283185307179586477
#define pi 3.141592653589793238
#define logrtpi 0.5 * log(tpi)

/* --- Internal functions (described below) ----------------- */

void gauss(double *z, int nkk);

void rt(double *z, int nkk, int dof);

void chol(int nkk, double **B);

void perm(double *work, int nkk);

double ltprob(int dof, double z, double *constt);

double lnormprob(int k, int nkk, int l, double ***mu, double ****B,
                 double *datai);

double det(int k, int nkk, int l, double ****B);

int read_mixture_params(char *fname, int kmax, int *nk, double **sig, int *Lk,
                        double **lambda, double ***mu, double ****B);

void rwn_within_model(int k1, int *nk, int nsweep2, FILE *fpl, FILE *fpcf,
                      FILE *fpad, double **sig, int dof, double **data);

void fit_mixture_from_samples(int k1, int *nk, double **data, double **sig,
                              double ***mu, double ****B, double **lambda,
                              FILE *fpcf, int *Lk);

void usage(char *invocation);

/* ---main program-------------------------- */

int main(int argc, char *argv[]) {

  /*---Section 1 Declare Variables -------------------------*/

  /* ---counting variables ------------------- */
  int nburn, nsokal, nkeep, keep;

  /* ---random no. variables ----------------- */
  double u, constt;

  /* ---State parameters and variables ------- */
  int k, kmax, nkk, nkkn, nkmax;
  int kn = 0;

  /* ---Mixture parameters --------------------*/
  int Lkk, Lkkn, Lkmax;
  int l = 0;
  int ln = 0;
  double tol = 0.00001;

  /* ---RWM parameters ------------------------*/
  double Z[1];
  double gamma = 0.0;

  /* ---Probabilities ------------------------ */
  double lp, logratio, llh, llhn;
  double lpn = 0.0;

  /* ---working arrays and variables --------- */
  double sum, thresh;

  /* ---autocorrelation variables ------------ */
  double var, tau;
  int m;

  /* ---adaptation parameters ---------------- */
  int reinit, nreinit;
  double pkllim;

  /* --- Section 2 - Read in Comand Line Variables ----------------- */

  clock_t starttime = clock();

  /* Definition of command line variables and explanation

     Prog variable ~ Command line variable ~ Explanation

     mode ~ m ~ 0 if mixture fitting, 1 if user supplied mixture params,
                2 if AutoRJ
     nsweep ~ N ~ no. of reversible jump sweeps in stage 3
     nsweep2 ~ n ~ max(n,10000*nk,100000) sweeps in within-model RWM in stage 1
     adapt ~ a ~ 1 if RJ adaptation done in stage 3, 0 otherwise
     doperm ~ p ~ 1 if random permutation done in stage 3 RJ move, 0 otherwise
     seed ~ s ~ random no. seed, 0 uses clock seed
     dof ~ t ~ 0 if Normal random variables used in RWM and RJ moves, otherwise
               specify integer degrees of freedom of student t variables
     fname ~ f ~ filename base */

  // Default values
  int nsweep = 1E5;
  int nsweep2 = 1E5;
  char fname_default[] = "output";
  char *fname = fname_default;
  int doperm = 1;
  unsigned long seed = 0;
  int mode = 0;
  int adapt = 1;
  int dof = 0;

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
    } else {
      printf("Unrecognized argument: %s\n", argv[i]);
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  sdrni(&seed);

  /* --- Section 3 - Initial File handling ---------------------  */
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

  /* Print user options to log file */

  fprintf(fpl, "seed: %ld\n", seed);
  fprintf(fpl, "m: %d\n", mode);
  fprintf(fpl, "a: %d\n", adapt);
  fprintf(fpl, "p: %d\n", doperm);
  fprintf(fpl, "n: %d\n", nsweep2);
  fprintf(fpl, "N: %d\n", nsweep);

  /* --- Section 4.0 - Read in key variables from user functions -*/

  getkmax(&kmax);
  if (kmax > kmaxmax) {
    printf("\nError:kmax too large \n");
    return 0;
  } else if (kmax < 0) {
    printf("\nError:negative kmax \n");
    return 0;
  }

  int *nk = (int *)malloc(kmax * sizeof(int));
  int *Lk = (int *)malloc(kmax * sizeof(int));
  int *ksummary = (int *)malloc(kmax * sizeof(int));

  getnk(kmax, nk);
  nkmax = nk[0];
  ksummary[0] = 0;
  for (int k1 = 1; k1 < kmax; k1++) {
    nkmax = max(nk[k1], nkmax);
    ksummary[k1] = 0;
  }

  double **lambda = (double **)malloc(kmax * sizeof(double *));
  double ***mu = (double ***)malloc(kmax * sizeof(double **));
  double ****B = (double ****)malloc(kmax * sizeof(double ***));
  double **detB = (double **)malloc(kmax * sizeof(double *));
  double **sig = (double **)malloc(kmax * sizeof(double *));
  for (int k1 = 0; k1 < kmax; k1++) {
    nkk = nk[k1];
    lambda[k1] = (double *)malloc(Lkmaxmax * sizeof(double));
    mu[k1] = (double **)malloc(Lkmaxmax * sizeof(double *));
    B[k1] = (double ***)malloc(Lkmaxmax * sizeof(double **));
    detB[k1] = (double *)malloc(Lkmaxmax * sizeof(double));
    sig[k1] = (double *)malloc(nkk * sizeof(double));
    for (int l1 = 0; l1 < Lkmaxmax; l1++) {
      mu[k1][l1] = (double *)malloc(nkk * sizeof(double));
      B[k1][l1] = (double **)malloc(nkk * sizeof(double *));
      for (int j1 = 0; j1 < nkk; j1++) {
        B[k1][l1][j1] = (double *)malloc(nkk * sizeof(double));
      }
    }
  }

  /* --- Section 5.1 - Read in mixture parameters if mode 1 (m=1) --- */

  /* These parameters are used if mode 1 (m=1) of the AutoMix sampler is
     used. Note that if the parameters are unavailable or inconsistent with
     the user supplied functions or unavailable go straight to section 5.2
     where initial within-model RWM are performed */

  if (mode == 1) {
    int ok = read_mixture_params(fname, kmax, nk, sig, Lk, lambda, mu, B);
    if (ok == EXIT_FAILURE) {
      mode = 0;
    }
  } else if (mode == 0 || mode == 2) {

    /* --- Section 5.2 - Within-model runs if mixture parameters unavailable -*/
    for (int k1 = 0; k1 < kmax; k1++) {
      int nkk = nk[k1];
      int lendata = 1000 * nkk;
      double **data = (double **)malloc(lendata * sizeof(double *));
      for (int i = 0; i < lendata; i++) {
        data[i] = (double *)malloc(nkk * sizeof(double));
      }
      /* --- Section 5.2.1 - RWM Within Model (Stage 1) -------*/
      rwn_within_model(k1, nk, nsweep2, fpl, fpcf, fpad, sig, dof, data);

      /* --- Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)- */
      /* Note only done if mode 0 (m=0) if mode m=2, go to section 5.2.3*/
      /* Mixture fitting done component wise EM algorithm described in
         Figueiredo and Jain, 2002 (see thesis for full reference) */

      printf("\nMixture Fitting: Model %d", k1 + 1);
      if (mode == 0) {
        fit_mixture_from_samples(k1, nk, data, sig, mu, B, lambda, fpcf, Lk);
      } else if (mode == 2) {
        /* --- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --*/
        /* Note only done if mode 2 (m=2).*/
        Lk[k1] = 1;
        lambda[k1][0] = 1.0;
        for (int j1 = 0; j1 < nkk; j1++) {
          mu[k1][0][j1] = 0.0;
          for (int i1 = 0; i1 < lendata; i1++) {
            mu[k1][0][j1] += data[i1][j1];
          }
          mu[k1][0][j1] /= ((double)lendata);
        }
        for (int j1 = 0; j1 < nkk; j1++) {
          for (int j2 = 0; j2 <= j1; j2++) {
            B[k1][0][j1][j2] = 0.0;
            for (int i1 = 0; i1 < lendata; i1++) {
              B[k1][0][j1][j2] += (data[i1][j1] - mu[k1][0][j1]) *
                                  (data[i1][j2] - mu[k1][0][j2]);
            }
            B[k1][0][j1][j2] /= ((double)(lendata - 1));
          }
        }
        chol(nkk, B[k1][0]);

        for (int i1 = 0; i1 < lendata; i1++) {
          free(data[i1]);
        }
        free(data);
      }
    }
  }

  /* Print mixture parameters to file (log and mix files) for reference
     and use in future runs. */

  sprintf(datafname, "%s_mix.data", fname);
  FILE *fpmix = fopen(datafname, "w");
  fprintf(fpmix, "%d\n", kmax);
  for (int k1 = 0; k1 < kmax; k1++) {
    fprintf(fpmix, "%d\n", nk[k1]);
  }

  for (int k1 = 0; k1 < kmax; k1++) {
    fprintf(fpl, "\nModel:%d\n", k1 + 1);
    Lkk = Lk[k1];
    nkk = nk[k1];
    fprintf(fpl, "\nARW params:\n");
    for (int j1 = 0; j1 < nkk; j1++) {
      fprintf(fpl, "%lf ", sig[k1][j1]);
    }
    fprintf(fpl, "\n");
    fprintf(fpl, "\nLkk:%d\n", Lkk);
    for (int j1 = 0; j1 < nkk; j1++) {
      fprintf(fpmix, "%lf\n", sig[k1][j1]);
    }
    fprintf(fpmix, "%d\n", Lkk);
    for (int l1 = 0; l1 < Lkk; l1++) {
      fprintf(fpl, "\nComponent:%d\n", l1 + 1);
      fprintf(fpl, "lambda:%lf\n", lambda[k1][l1]);
      fprintf(fpmix, "%lf\n", lambda[k1][l1]);
      fprintf(fpl, "mu:\n");
      for (int j1 = 0; j1 < nkk; j1++) {
        fprintf(fpl, "%lf ", mu[k1][l1][j1]);
        fprintf(fpmix, "%lf\n", mu[k1][l1][j1]);
      }
      fprintf(fpl, "\nB:\n");
      for (int j1 = 0; j1 < nkk; j1++) {
        for (int j2 = 0; j2 <= j1; j2++) {
          fprintf(fpl, "%lf ", B[k1][l1][j1][j2]);
          fprintf(fpmix, "%lf\n", B[k1][l1][j1][j2]);
        }
        fprintf(fpl, "\n");
      }
      detB[k1][l1] = det(k1, nkk, l1, B);
    }
  }
  fflush(NULL);

  /* --Section 6 - Secondary file handling -------------*/

  sprintf(datafname, "%s_k.data", fname);
  FILE *fpk = fopen(datafname, "w");
  sprintf(datafname, "%s_lp.data", fname);
  FILE *fplp = fopen(datafname, "w");

  FILE *fpt[kmaxmax];
  for (int k1 = 0; k1 < kmax; k1++) {
    sprintf(datafname, "%s_theta%d.data", fname, k1 + 1);
    fpt[k1] = fopen(datafname, "w");
  }
  free(datafname);

  /* --Section 7 - Final initialisation of variables ----*/

  int naccrwmb = 0;
  int ntryrwmb = 0;
  int naccrwms = 0;
  int ntryrwms = 0;
  int nacctd = 0;
  int ntrytd = 0;

  constt = 100000.0;
  Lkmax = Lk[0];
  for (int k1 = 1; k1 < kmax; k1++) {
    Lkmax = max(Lkmax, Lk[k1]);
  }
  k = (int)floor(kmax * sdrand());
  nkk = nk[k];
  Lkk = Lk[k];

  double *theta = (double *)malloc(nkmax * sizeof(double));
  double *thetan = (double *)malloc(nkmax * sizeof(double));
  double *work = (double *)malloc(nkmax * sizeof(double));
  double *palloc = (double *)malloc(Lkmax * sizeof(double));
  double *pallocn = (double *)malloc(Lkmax * sizeof(double));
  double *Znkk = (double *)malloc(nkmax * sizeof(double));
  double *propk = (double *)malloc(kmax * sizeof(double));
  double *pk = (double *)malloc(kmax * sizeof(double));

  getic(k, nkk, theta);

  lp = lpost(k, nkk, theta, &llh);

  for (int k1 = 0; k1 < kmax; k1++) {
    pk[k1] = 1.0 / kmax;
    if (k1 == k) {
      propk[k1] = 1.0;
    } else {
      propk[k1] = 0.0;
    }
  }
  nreinit = 1;
  reinit = 0;
  pkllim = 1.0 / 10.0;

  tol = 0.5 / nsweep;
  nburn = max(10000, (int)(nsweep / 10));

  nsokal = 1;
  nkeep = nsweep / (2 * nsokal);
  nkeep = (int)pow(2.0, min(15, (int)(log(nkeep) / log(2.0) + 0.001)));
  keep = nburn + (nsweep - nkeep * nsokal);
  double *xr = (double *)malloc(nkeep * sizeof(double));

  /* -----Start of main loop ----------------*/
  for (int sweep = 1; sweep <= (nburn + nsweep); sweep++) {

    /* --Section 8 - RWM within-model moves ---*/

    /* Every 10 sweeps to block RWM */
    if (fmod(sweep, 10) < 0.05) {
      ntryrwmb++;
      rt(Znkk, nkk, dof);
      for (int j1 = 0; j1 < nkk; j1++) {
        thetan[j1] = theta[j1] + sig[k][j1] * Znkk[j1];
      }
      lpn = lpost(k, nkk, thetan, &llhn);
      if (sdrand() < exp(max(-30.0, min(0.0, lpn - lp)))) {
        naccrwmb++;
        for (int j1 = 0; j1 < nkk; j1++) {
          theta[j1] = thetan[j1];
        }
        lp = lpn;
        llh = llhn;
      }
    } else {
      /* else do component-wise RWM */
      for (int j1 = 0; j1 < nkk; j1++) {
        thetan[j1] = theta[j1];
      }
      for (int j1 = 0; j1 < nkk; j1++) {
        ntryrwms++;
        rt(Z, 1, dof);
        thetan[j1] = theta[j1] + sig[k][j1] * Z[0];
        lpn = lpost(k, nkk, thetan, &llhn);
        if (sdrand() < exp(max(-30.0, min(0.0, lpn - lp)))) {
          naccrwms++;
          theta[j1] = thetan[j1];
          lp = lpn;
          llh = llhn;
        } else {
          thetan[j1] = theta[j1];
        }
      }
    }

    /* --- Section 9 Reversible Jump Moves -------------------*/

    /* --Section 9.1 - Allocate current position to a component --*/

    ntrytd++;
    if (Lkk > 1) {
      sum = 0.0;
      for (int l1 = 0; l1 < Lkk; l1++) {
        palloc[l1] = log(lambda[k][l1]) + lnormprob(k, nkk, l1, mu, B, theta);
        palloc[l1] = exp(palloc[l1]);
        sum += palloc[l1];
      }
      if (sum > 0) {
        for (int l1 = 0; l1 < Lkk; l1++) {
          palloc[l1] /= sum;
        }
      } else {
        for (int l1 = 0; l1 < Lkk; l1++) {
          palloc[l1] = 1.0 / Lkk;
        }
      }
      u = sdrand();
      thresh = 0.0;
      for (int l1 = 0; l1 < Lkk; l1++) {
        thresh += palloc[l1];
        if (u < thresh) {
          l = l1;
          break;
        }
      }
    } else {
      l = 0;
      palloc[l] = 1.0;
    }

    /* --Section 9.2 - Standardise state variable --------------- */

    for (int j1 = 0; j1 < nkk; j1++) {
      work[j1] = theta[j1] - mu[k][l][j1];
    }
    for (int j1 = 0; j1 < nkk; j1++) {
      for (int j2 = 0; j2 < j1; j2++) {
        work[j1] = work[j1] - B[k][l][j1][j2] * work[j2];
      }
      work[j1] = work[j1] / B[k][l][j1][j1];
    }

    /* --Section 9.3 - Choose proposed new model and component ----*/

    if (kmax == 1) {
      kn = k;
      logratio = 0.0;
    } else {
      gamma = pow(1.0 / (sweep + 1), (2.0 / 3.0));
      u = sdrand();
      thresh = 0.0;
      for (int k1 = 0; k1 < kmax; k1++) {
        thresh += pk[k1];
        if (u < thresh) {
          kn = k1;
          break;
        }
      }
      logratio = log(pk[k]) - log(pk[kn]);
    }

    nkkn = nk[kn];
    Lkkn = Lk[kn];

    u = sdrand();
    thresh = 0.0;
    for (int l1 = 0; l1 < Lkkn; l1++) {
      thresh += lambda[kn][l1];
      if (u < thresh) {
        ln = l1;
        break;
      }
    }

    /* --Section 9.4 Propose new state ----------------*/

    if (nkk < nkkn) {
      rt(&(work[nkk]), nkkn - nkk, dof);
      if (dof > 0) {
        for (int j1 = nkk; j1 < nkkn; j1++) {
          logratio -= ltprob(dof, work[j1], &constt);
        }
      } else {
        for (int j1 = nkk; j1 < nkkn; j1++) {
          logratio += 0.5 * pow(work[j1], 2.0) + logrtpi;
        }
      }
      if (doperm) {
        perm(work, nkkn);
      }
    } else if (nkk == nkkn) {
      if (doperm) {
        perm(work, nkk);
      }
    } else {
      if (doperm) {
        perm(work, nkk);
      }
      if (dof > 0) {
        for (int j1 = nkkn; j1 < nkk; j1++) {
          logratio += ltprob(dof, work[j1], &constt);
        }
      } else {
        for (int j1 = nkkn; j1 < nkk; j1++) {
          logratio -= (0.5 * pow(work[j1], 2.0) + logrtpi);
        }
      }
    }

    for (int j1 = 0; j1 < nkkn; j1++) {
      thetan[j1] = mu[kn][ln][j1];
      for (int j2 = 0; j2 <= j1; j2++) {
        thetan[j1] += B[kn][ln][j1][j2] * work[j2];
      }
    }

    /* --Section 9.5 - Work out probability of allocating to component
       for acceptance ratio (reverse move) ---------*/

    if (Lkkn > 1) {
      sum = 0.0;
      for (int l1 = 0; l1 < Lkkn; l1++) {
        pallocn[l1] =
            log(lambda[kn][l1]) + lnormprob(kn, nkkn, l1, mu, B, thetan);
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

    /* --Section 9.6 - Work out acceptance probability  and new state --*/

    lpn = lpost(kn, nkkn, thetan, &llhn);

    logratio += (lpn - lp);
    logratio += (log(pallocn[ln]) - log(palloc[l]));
    logratio += (log(lambda[k][l]) - log(lambda[kn][ln]));
    logratio += (log(detB[kn][ln]) - log(detB[k][l]));

    if (sdrand() < exp(max(-30.0, min(0.0, logratio)))) {
      for (int j1 = 0; j1 < nkkn; j1++) {
        theta[j1] = thetan[j1];
      }
      lp = lpn;
      llh = llhn;
      k = kn;
      nkk = nkkn;
      Lkk = Lkkn;
      nacctd++;
    }

    if (adapt == 1) {
      if (sweep > nburn) {
        for (int k1 = 0; k1 < kmax; k1++) {
          if (k1 == k) {
            propk[k1] = 1.0;
          } else {
            propk[k1] = 0.0;
          }
          pk[k1] += (gamma * (propk[k1] - pk[k1]));
        }
        for (int k1 = 0; k1 < kmax; k1++) {
          if (pk[k1] < pkllim) {
            reinit = 1;
          }
        }
        if (reinit == 1) {
          reinit = 0;
          nreinit++;
          pkllim = 1.0 / (10.0 * nreinit);
          for (int k1 = 0; k1 < kmax; k1++) {
            pk[k1] = 1.0 / kmax;
          }
        }
      }
    }

    if (sweep > nburn) {
      (ksummary[k])++;
    }
    /* --- Section 10 - Write variables to files --------- */

    if (sweep > nburn) {

      fprintf(fpk, "%d\n", k + 1);
      fprintf(fplp, "%lf %lf\n", lp, llh);
      for (int k1 = 0; k1 < kmax; k1++) {
        fprintf(fpp, "%lf ", pk[k1]);
      }
      fprintf(fpp, "\n");
      for (int j1 = 0; j1 < nkk; j1++) {
        fprintf(fpt[k], "%lf ", theta[j1]);
      }
      fprintf(fpt[k], "\n");
    }

    if (sweep > keep && fmod(sweep - keep, nsokal) < 0.005) {
      xr[((sweep - keep) / nsokal) - 1] = k;
    }

    if (sweep == 1) {
      printf("\nBurning in");
    }
    if ((sweep <= nburn) && (fmod(sweep, (nburn / 10)) < tol)) {
      printf(" .");
      fflush(NULL);
    }

    if (sweep == (nburn + 1)) {
      printf("\nStart of main sample:");
    }

    if ((sweep > nburn) && (fmod(sweep - nburn, (nsweep / 10)) < tol)) {
      printf("\nNo. of iterations remaining: %d", nsweep + nburn - sweep);
    }
    fflush(NULL);
  }
  printf("\n");

  /* --- Section 11 - Write log file ----------------------*/

  sokal(nkeep, xr, &var, &tau, &m);
  fprintf(fpl, "\nAutocorrelation Time:\n");
  fprintf(fpl, "nkeep:%d, nsokal:%d, var:%lf, tau:%lf\n", nkeep, nsokal, var,
          tau);
  for (int i1 = 0; i1 < m; i1++) {
    fprintf(fpac, "%lf\n", xr[i1]);
  }

  fprintf(fpl, "\nPosterior Model Probabilities:\n");
  for (int k1 = 0; k1 < kmax; k1++) {
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

void gauss(double *z, int n) {

  /* Uses Box mueller method to simulate n N(0,1) variables and stores them
     in z */

  int n1;
  double u, v;

  n1 = n - 1;
  for (int i = 0; i < n1; i += 2) {
    u = sqrt(-2.0 * log(sdrand()));
    v = tpi * sdrand();
    z[i] = u * sin(v);
    z[i + 1] = u * cos(v);
  }
  if (fmod(n, 2) < 0.5) {
    return;
  } else {
    u = sqrt(-2.0 * log(sdrand()));
    v = tpi * sdrand();
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

void chol(int nkk, double **A) {

  /* Performs cholesky decompositon of A and returns result in the
     same matrix - adapted from PJG Fortran function*/

  for (int j1 = 0; j1 < nkk; j1++) {
    double sum = A[j1][j1];
    for (int j2 = 0; j2 < j1; j2++) {
      sum -= pow(A[j1][j2], 2);
    }
    A[j1][j1] = sqrt(sum);

    for (int j2 = j1 + 1; j2 < nkk; j2++) {
      sum = A[j2][j1];
      for (int j3 = 0; j3 < j1; j3++) {
        sum -= A[j2][j3] * A[j1][j3];
      }
      A[j2][j1] = sum / A[j1][j1];
    }
  }
}

void perm(double *work, int nkk) {

  /* Randomly permutes the nkk-vector work */

  for (int j1 = 0; j1 < (nkk - 1); j1++) {
    int j2 = j1 + (int)((nkk - j1) * sdrand());
    if (j2 != j1) {
      double temp = work[j2];
      work[j2] = work[j1];
      work[j1] = temp;
    }
  }
  return;
}

double ltprob(int dof, double z, double *constt) {

  /* Evaluates the log of p.d.f. of a t variable with dof degrees of freedom
     at point z */

  double out;

  /* only calculate const of proportionality once */
  if ((*constt) > 10000.0) {
    *constt =
        loggamma(0.5 * (dof + 1)) - loggamma(0.5 * dof) - 0.5 * log(dof * pi);
  }
  out = (*constt) - 0.5 * (dof + 1) * log(1.0 + pow(z, 2.0) / dof);
  return out;
}

double lnormprob(int k, int nkk, int l, double ***mu, double ****B,
                 double *datai) {

  /* Evaluates log of p.d.f. for a multivariate normal for model
     k, of dimension nkk, component l. The summary of means and
     sqrt of cov matrices (for all models and all component)
     are supplied in mu and B */

  double work[nkk];
  double out;

  for (int j1 = 0; j1 < nkk; j1++) {
    work[j1] = datai[j1] - mu[k][l][j1];
  }
  for (int j1 = 0; j1 < nkk; j1++) {
    for (int j2 = 0; j2 < j1; j2++) {
      (work[j1]) -= B[k][l][j1][j2] * work[j2];
    }
    (work[j1]) /= B[k][l][j1][j1];
  }
  out = 0.0;
  for (int j1 = 0; j1 < nkk; j1++) {
    out += (work[j1] * work[j1]);
  }
  out = -0.5 * out - (nkk / 2.0) * log(tpi) - log(det(k, nkk, l, B));
  return out;
}

double det(int k, int nkk, int l, double ****B) {

  /* Evaluates the determinant of a matrix in B corresponding to model k,
     component l. */
  double out = 1.0;
  for (int j1 = 0; j1 < nkk; j1++) {
    out *= B[k][l][j1][j1];
  }
  return out;
}

int read_mixture_params(char *fname, int kmax, int *nk, double **sig, int *Lk,
                        double **lambda, double ***mu, double ****B) {
  /* Check user has supplied mixture parameters if trying to use mode 1.
     If not default back to mode 0 */
  char *datafname = (char *)malloc((strlen(fname) + 20) * sizeof(*datafname));
  sprintf(datafname, "%s_mix.data", fname);
  FILE *fpmix = fopen(datafname, "r");
  free(datafname);
  if (fpmix == NULL) {
    printf("\nProblem opening mixture file:");
    printf("\nContinuing using RWM to estimate parameters");
    return EXIT_FAILURE;
  }
  int k1;
  if (fscanf(fpmix, "%d", &k1) == EOF) {
    printf("\nEnd of file encountered before parameters read:");
    printf("\nContinuing using RWM to estimate parameters");
    return EXIT_FAILURE;
  }
  if (k1 != kmax) {
    printf("\nFile kmax contradicts getkmax function:");
    printf("\nContinuing using RWM to estimate parameters");
    return EXIT_FAILURE;
  }
  for (int k = 0; k < kmax; k++) {
    int nkk;
    if (fscanf(fpmix, "%d", &nkk) == EOF) {
      printf("\nEnd of file encountered before parameters read:");
      printf("\nContinuing using RWM to estimate parameters");
      return EXIT_FAILURE;
    }
    if (nkk != nk[k]) {
      printf("\nFile kmax contradicts getnk function:");
      printf("\nContinuing using RWM to estimate parameters");
      return EXIT_FAILURE;
    }
  }
  for (int k = 0; k < kmax; k++) {
    int nkk = nk[k];
    for (int l = 0; l < nkk; l++) {
      if (fscanf(fpmix, "%lf", &(sig[k][l])) == EOF) {
        printf("\nEnd of file encountered before parameters read:");
        printf("\nContinuing using RWM to estimate parameters");
        return EXIT_FAILURE;
      }
    }
    if (fscanf(fpmix, "%d", &(Lk[k])) == EOF) {
      printf("\nEnd of file encountered before parameters read:");
      printf("\nContinuing using RWM to estimate parameters");
      return EXIT_FAILURE;
    }
    int Lkk = Lk[k];
    for (int l = 0; l < Lkk; l++) {
      if (fscanf(fpmix, "%lf", &(lambda[k][l])) == EOF) {
        printf("\nEnd of file encountered before parameters read:");
        printf("\nContinuing using RWM to estimate parameters");
        return EXIT_FAILURE;
      }
      for (int i = 0; i < nkk; i++) {
        if (fscanf(fpmix, "%lf", &(mu[k][l][i])) == EOF) {
          printf("\nEnd of file encountered before parameters read:");
          printf("\nContinuing using RWM to estimate parameters");
          return EXIT_FAILURE;
        }
      }
      for (int i = 0; i < nkk; i++) {
        for (int j = 0; j <= i; j++) {
          if (fscanf(fpmix, "%lf", &(B[k][l][i][j])) == EOF) {
            printf("\nEnd of file encountered before parameters read:");
            printf("\nContinuing using RWM to estimate parameters");
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
      printf("\nContinuing using RWM to estimate parameters");
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

void rwn_within_model(int k1, int *nk, int nsweep2, FILE *fpl, FILE *fpcf,
                      FILE *fpad, double **sig, int dof, double **data) {
  int nkk = nk[k1];
  int nsweepr = max(nsweep2, 10000 * nkk);
  int nburn = nsweepr / 10;
  double alphastar = 0.25;
  nsweepr += nburn;
  double *rwm = (double *)malloc(nkk * sizeof(double));
  double *rwmn = (double *)malloc(nkk * sizeof(double));
  int *nacc = (int *)malloc(nkk * sizeof(int));
  int *ntry = (int *)malloc(nkk * sizeof(int));
  double *Znkk = (double *)malloc(nkk * sizeof(double));

  printf("\nRWM for Model %d", k1 + 1);
  fprintf(fpl, "\nRWM for Model %d", k1 + 1);
  fprintf(fpcf, "RWM for Model %d\n", k1 + 1);
  fprintf(fpad, "RWM for Model %d\n", k1 + 1);
  fflush(NULL);
  getic(k1, nkk, rwm);
  for (int j1 = 0; j1 < nkk; j1++) {
    rwmn[j1] = rwm[j1];
    sig[k1][j1] = 10.0;
    nacc[j1] = 0;
    ntry[j1] = 0;
  }
  double llh, llhn, lpn;
  double lp = lpost(k1, nkk, rwm, &llh);

  int i2 = 0;
  int remain = nsweepr;
  double tol = 1E-5;
  for (int sweep = 1; sweep <= nsweepr; sweep++) {
    remain--;
    if ((sweep >= nburn) &&
        (fmod((sweep - nburn), ((nsweepr - nburn) / 10)) < tol)) {
      printf("\nNo. of iterations remaining: %d", remain);
      fflush(NULL);
    }
    double u = sdrand();
    if (sweep > nburn && u < 0.1) {
      rt(Znkk, nkk, dof);
      for (int i = 0; i < nkk; i++) {
        rwmn[i] = rwm[i] + sig[k1][i] * Znkk[i];
      }
      lpn = lpost(k1, nkk, rwmn, &llhn);
      if (sdrand() < exp(max(-30.0, min(0.0, lpn - lp)))) {
        for (int i = 0; i < nkk; i++) {
          rwm[i] = rwmn[i];
        }
        lp = lpn;
        llh = llhn;
      }
    } else {
      double gamma = 10.0 * pow(1.0 / (sweep + 1), 2.0 / 3.0);
      for (int i = 0; i < nkk; i++) {
        rwmn[i] = rwm[i];
      }
      for (int i = 0; i < nkk; i++) {
        double Z;
        rt(&Z, 1, dof);
        rwmn[i] = rwm[i] + sig[k1][i] * Z;
        lpn = lpost(k1, nkk, rwmn, &llhn);
        double accept = min(1, exp(max(-30.0, min(0.0, lpn - lp))));
        if (sdrand() < accept) {
          (nacc[i])++;
          (ntry[i])++;
          rwm[i] = rwmn[i];
          lp = lpn;
          llh = llhn;
          sig[k1][i] = max(0, sig[k1][i] - gamma * (alphastar - 1));
        } else {
          (ntry[i])++;
          rwmn[i] = rwm[i];
          sig[k1][i] = max(0, sig[k1][i] - gamma * (alphastar));
        }
      }
    }
    if (remain < (10000 * nkk) && fmod(remain, 10.0) < 0.05) {
      for (int i = 0; i < nkk; i++) {
        data[i2][i] = rwm[i];
      }
      i2++;
    }
    if (fmod(sweep, 100.0) < 0.05) {
      for (int i = 0; i < nkk; i++) {
        fprintf(fpad, "%lf %lf ", sig[k1][i],
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

void fit_mixture_from_samples(int k1, int *nk, double **data, double **sig,
                              double ***mu, double ****B, double **lambda,
                              FILE *fpcf, int *Lk) {
  int Lkk = Lkmaxmax;
  int *init = (int *)malloc(Lkk * sizeof(int));
  int nkk = nk[k1];
  int lendata = 1000 * nkk;
  int l1 = 0;
  double lpn = 0.0;
  double costfn = 0.0;
  double costfnmin = 0.0;
  double tol = 1E-5;

  double ***BBT = (double ***)malloc(Lkmaxmax * sizeof(double **));
  for (int l1 = 0; l1 < Lkmaxmax; l1++) {
    BBT[l1] = (double **)malloc(nkk * sizeof(double *));
    for (int j1 = 0; j1 < nkk; j1++) {
      BBT[l1][j1] = (double *)malloc(nkk * sizeof(double));
    }
  }
  double ***Bmin = (double ***)malloc(Lkmaxmax * sizeof(double **));
  for (int l1 = 0; l1 < Lkmaxmax; l1++) {
    Bmin[l1] = (double **)malloc(nkk * sizeof(double *));
    for (int j1 = 0; j1 < nkk; j1++) {
      Bmin[l1][j1] = (double *)malloc(nkk * sizeof(double));
    }
  }
  double *lambdamin = (double *)malloc(Lkmaxmax * sizeof(double));
  double **mumin = (double **)malloc(Lkmaxmax * sizeof(double *));
  for (int l1 = 0; l1 < Lkmaxmax; l1++) {
    mumin[l1] = (double *)malloc(nkk * sizeof(double));
  }

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

  double *datamean = (double *)malloc(nkk * sizeof(double));
  double **M1 = (double **)malloc(nkk * sizeof(double *));
  for (int j1 = 0; j1 < nkk; j1++) {
    M1[j1] = (double *)malloc(nkk * sizeof(double));
  }
  for (int j1 = 0; j1 < nkk; j1++) {
    datamean[j1] = 0.0;
    for (int i1 = 0; i1 < lendata; i1++) {
      datamean[j1] += data[i1][j1];
    }
    datamean[j1] /= ((double)lendata);
  }
  for (int j1 = 0; j1 < nkk; j1++) {
    for (int j2 = 0; j2 < nkk; j2++) {
      M1[j1][j2] = 0;
      for (int i1 = 0; i1 < lendata; i1++) {
        M1[j1][j2] +=
            (data[i1][j1] - datamean[j1]) * (data[i1][j2] - datamean[j2]);
      }
      M1[j1][j2] /= ((double)lendata);
    }
  }
  double sigma = 0.0;
  for (int j1 = 0; j1 < nkk; j1++) {
    sigma += M1[j1][j1];
  }
  sigma /= (10.0 * nkk);

  for (int l1 = 0; l1 < Lkk; l1++) {
    for (int j1 = 0; j1 < nkk; j1++) {
      mu[k1][l1][j1] = data[init[l1]][j1];
      BBT[l1][j1][j1] = sigma;
      B[k1][l1][j1][j1] = BBT[l1][j1][j1];
      for (int j2 = 0; j2 < j1; j2++) {
        BBT[l1][j1][j2] = 0.0;
        B[k1][l1][j1][j2] = BBT[l1][j1][j2];
      }
    }
    chol(nkk, B[k1][l1]);
    lambda[k1][l1] = 1.0 / Lkk;
  }

  double **w = (double **)malloc(lendata * sizeof(double *));
  double *logw = (double *)malloc(Lkk * sizeof(double));
  double **lpdatagivenl = (double **)malloc(lendata * sizeof(double *));
  for (int i1 = 0; i1 < lendata; i1++) {
    w[i1] = (double *)malloc(Lkk * sizeof(double));
    lpdatagivenl[i1] = (double *)malloc(Lkk * sizeof(double));
  }

  for (int i1 = 0; i1 < lendata; i1++) {
    double sum = 0.0;
    for (int l1 = 0; l1 < Lkk; l1++) {
      lpdatagivenl[i1][l1] = lnormprob(k1, nkk, l1, mu, B, data[i1]);
      logw[l1] = log(lambda[k1][l1]) + lpdatagivenl[i1][l1];
      w[i1][l1] = exp(logw[l1]);
      sum += w[i1][l1];
    }
    for (int l1 = 0; l1 < Lkk; l1++) {
      w[i1][l1] /= sum;
    }
  }

  double *sumw = (double *)malloc(Lkk * sizeof(double));

  int stop = 0;
  int count = 0;
  double wnewl1 = 0.0;
  int nparams = nkk + (nkk * (nkk + 1)) / 2;
  int Lkkmin = 0;

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
      lambda[k1][l1] = wnewl1 / sumwnew;
      double sumlambda = 0.0;
      for (int l2 = 0; l2 < Lkk; l2++) {
        sumlambda += lambda[k1][l2];
      }
      for (int l2 = 0; l2 < Lkk; l2++) {
        lambda[k1][l2] /= sumlambda;
      }

      if (lambda[k1][l1] > 0.005) {
        /*changed to 0.005 from 0.0 -renormalise else */
        for (int j1 = 0; j1 < nkk; j1++) {
          mu[k1][l1][j1] = 0.0;
          for (int i1 = 0; i1 < lendata; i1++) {
            mu[k1][l1][j1] += data[i1][j1] * w[i1][l1];
          }
          mu[k1][l1][j1] /= sumw[l1];

          for (int j2 = 0; j2 <= j1; j2++) {
            BBT[l1][j1][j2] = 0.0;
            for (int i1 = 0; i1 < lendata; i1++) {
              BBT[l1][j1][j2] += (data[i1][j1] - mu[k1][l1][j1]) *
                                 (data[i1][j2] - mu[k1][l1][j2]) * w[i1][l1];
            }
            BBT[l1][j1][j2] /= sumw[l1];
            B[k1][l1][j1][j2] = BBT[l1][j1][j2];
          }
        }

        chol(nkk, B[k1][l1]);

        for (int i1 = 0; i1 < lendata; i1++) {
          lpdatagivenl[i1][l1] = lnormprob(k1, nkk, l1, mu, B, data[i1]);
        }
        l1++;

      } else {
        if (fmod(Lkk, 5) < 0.05) {
          printf("\n");
        }
        printf("%d(%d-n) ", Lkk, count);
        natann = 1;
        if (l1 < (Lkk - 1)) {
          for (int l2 = l1; l2 < (Lkk - 1); l2++) {
            lambda[k1][l2] = lambda[k1][l2 + 1];
            for (int j1 = 0; j1 < nkk; j1++) {
              mu[k1][l2][j1] = mu[k1][l2 + 1][j1];
              for (int j2 = 0; j2 <= j1; j2++) {
                BBT[l2][j1][j2] = BBT[l2 + 1][j1][j2];
                B[k1][l2][j1][j2] = B[k1][l2 + 1][j1][j2];
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
          sumlambda += lambda[k1][l2];
        }
        for (int l2 = 0; l2 < Lkk; l2++) {
          lambda[k1][l2] /= sumlambda;
        }
      }

      lpn = 0.0;
      for (int i1 = 0; i1 < lendata; i1++) {
        double sum = 0.0;
        for (int l2 = 0; l2 < Lkk; l2++) {
          logw[l2] = log(lambda[k1][l2]) + lpdatagivenl[i1][l2];
          w[i1][l2] = exp(logw[l2]);
          sum += w[i1][l2];
        }
        if (sum > 0) {
          for (int l2 = 0; l2 < Lkk; l2++) {
            w[i1][l2] /= sum;
          }
          lpn += log(sum);
        } else {
          /* if no component fits point well make equally likely */
          for (int l2 = 0; l2 < Lkk; l2++) {
            w[i1][l2] = 1.0 / Lkk;
          }
          lpn += (-500.0);
        }
      }
    }

    double sum = 0.0;
    for (int l1 = 0; l1 < Lkk; l1++) {
      sum += log(lendata * lambda[k1][l1] / 12.0);
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
        lambdamin[l1] = lambda[k1][l1];
        for (int j1 = 0; j1 < nkk; j1++) {
          mumin[l1][j1] = mu[k1][l1][j1];
          for (int j2 = 0; j2 <= j1; j2++) {
            Bmin[l1][j1][j2] = B[k1][l1][j1][j2];
          }
        }
      }
    }
    if ((fabs(costfn - costfnnew) < min(tol * fabs(costfn), 0.01)) &&
        (count > 1)) {
      if (Lkk == 1) {
        stop = 1;
      } else {
        if (fmod(Lkk, 5) < 0.05) {
          printf("\n");
        }
        printf("%d(%d-f) ", Lkk, count);
        forceann = 2;
        double minlambda = lambda[k1][0];
        int ldel = 0;
        for (int l1 = 1; l1 < Lkk; l1++) {
          if (minlambda > lambda[k1][l1]) {
            minlambda = lambda[k1][l1];
            ldel = l1;
          }
        }
        if (ldel < (Lkk - 1)) {
          for (int l1 = ldel; l1 < (Lkk - 1); l1++) {
            lambda[k1][l1] = lambda[k1][l1 + 1];
            for (int j1 = 0; j1 < nkk; j1++) {
              mu[k1][l1][j1] = mu[k1][l1 + 1][j1];
              for (int j2 = 0; j2 <= j1; j2++) {
                BBT[l1][j1][j2] = BBT[l1 + 1][j1][j2];
                B[k1][l1][j1][j2] = B[k1][l1 + 1][j1][j2];
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
          sumlambda += lambda[k1][l1];
        }
        for (int l1 = 0; l1 < Lkk; l1++) {
          lambda[k1][l1] /= sumlambda;
        }

        lpn = 0.0;
        for (int i1 = 0; i1 < lendata; i1++) {
          sum = 0.0;
          for (int l2 = 0; l2 < Lkk; l2++) {
            logw[l2] = log(lambda[k1][l2]) + lpdatagivenl[i1][l2];
            w[i1][l2] = exp(logw[l2]);
            sum += w[i1][l2];
          }
          if (sum > 0) {
            for (int l2 = 0; l2 < Lkk; l2++) {
              w[i1][l2] /= sum;
            }
            lpn += log(sum);
          } else {
            /* if no component fits point well make equally likely */
            for (int l2 = 0; l2 < Lkk; l2++) {
              w[i1][l2] = 1.0 / Lkk;
            }
            lpn += (-500.0);
          }
        }

        sum = 0.0;
        for (int l1 = 0; l1 < Lkk; l1++) {
          sum += log(lendata * lambda[k1][l1] / 12.0);
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

  for (int j1 = 0; j1 < nkk; j1++) {
    free(M1[j1]);
  }
  free(M1);
  free(datamean);
  for (int i1 = 0; i1 < lendata; i1++) {
    free(data[i1]);
    free(w[i1]);
    free(lpdatagivenl[i1]);
  }
  free(w);
  free(lpdatagivenl);
  free(data);
  free(logw);
  free(sumw);
  free(init);
  Lk[k1] = Lkkmin;
  for (int l1 = 0; l1 < Lkkmin; l1++) {
    lambda[k1][l1] = lambdamin[l1];
    for (int j1 = 0; j1 < nkk; j1++) {
      mu[k1][l1][j1] = mumin[l1][j1];
    }
    for (int j1 = 0; j1 < nkk; j1++) {
      for (int j2 = 0; j2 <= j1; j2++) {
        B[k1][l1][j1][j2] = Bmin[l1][j1][j2];
      }
    }
  }
}

void usage(char *invocation) {
  char *name = strrchr(invocation, '/');
  if (name == NULL) {
    name = invocation;
  } else {
    name += 1;
  }
  printf("Usage: %s [-m int] [-N int] [-n int] [-a bool] \
           [-p bool] [-s int] [-t int] [-f string] [-h, --help]\n",
         name);
  printf(
      "-m int: Specify mode. 0 if mixture fitting, 1 if user supplied mixture params, \
           2 if AutoRJ.\n");
  printf("-N int: Number of reversible jump sweeps in stage 3.\n");
  printf("-n int: max(n, 10000 * nk, 100000) sweeps in within-model RWM in "
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
