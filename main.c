// Driver main function for the AutoMix RJMCMC sampler
#define VERSION "1.3"

#include "automix.h"
#include "logwrite.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define max(A, B) ((A) > (B) ? (A) : (B))

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
  int nburn = max(10000, (int)(nsweep / 10));

  // Initialize the Proposal (jumping) Distribution
  proposalDist jd;
  initJD(&jd);
  // Struct to hold run statistic variables
  runStats st;
  initializeRunStats(&st, nsweep, nsweep2, nburn, jd);

  // vector of adapted RWM scale parameters for each model k
  double **sig = (double **)malloc(jd.nmodels * sizeof(double *));
  for (int k = 0; k < jd.nmodels; k++) {
    int mdim = jd.model_dims[k];
    sig[k] = (double *)malloc(mdim * sizeof(double));
  }

  // --- Section 5.1 - Read in mixture parameters if mode 1 (m=1) ---
  if (mode == 1) {
    // Read AutoMix parameters from file if mode = 1
    int ok = read_mixture_params(fname, jd, sig);
    if (ok == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
  } else {
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
      rwm_within_model(model_k, jd.model_dims, nsweep2, st, sig[model_k], dof,
                       samples);
      printf("\nMixture Fitting: Model %d", model_k + 1);
      if (mode == 0) {
        // Section 5.2.2 - Fit Mixture to within-model sample, (stage 2)
        // Fit a Normal mixture distribution to the conditional target
        // distributions pi(theta_k|k). See theis, p 144.
        fit_mixture_from_samples(model_k, jd, samples, nsamples, &st);
      }
      if (mode == 2) {
        //--- Section 5.2.3 - Fit AutoRJ single mu vector and B matrix --
        fit_autorj(model_k, jd, samples, nsamples);
      }
      free(samples[0]);
      free(samples);
    }
  }

  // Write adaptation statistics to file
  write_adapt_to_file(fname, mode, jd, st);
  // Write mixture parameters to file
  write_mix_to_file(fname, jd, sig);
  // Write cf statistics to file
  write_cf_to_file(fname, mode, jd, st);

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
  int xr_i = 0;
  printf("Start of main sample:");
  ch.isBurning = 0;
  for (int sweep = 1; sweep <= nsweep; sweep++) {
    // Every 10 sweeps to block RWM
    ch.doBlockRWM = (sweep + nburn % 10 == 0);
    ch.gamma_sweep = pow(1.0 / (sweep + nburn + 1), (2.0 / 3.0));

    reversible_jump_move(&ch, jd, dof, &st, sig);

    (st.ksummary[ch.current_model_k])++;
    st.k_which_summary[sweep - 1] = ch.current_model_k + 1;
    st.logp_summary[sweep - 1][0] = ch.log_posterior;
    st.logp_summary[sweep - 1][1] = ch.log_likelihood;
    for (int k1 = 0; k1 < jd.nmodels; k1++) {
      st.pk_summary[sweep - 1][k1] = ch.pk[k1];
    }
    write_theta_to_file(fname, ch.current_model_k, ch.mdim, ch.theta);

    if (sweep > st.keep && ((sweep - st.keep) % st.nsokal == 0)) {
      st.xr[xr_i++] = ch.current_model_k;
    }
    // Print what's below about 10 times
    if (sweep % (nsweep / 10) == 0) {
      printf("\nNo. of iterations remaining: %d", nsweep - sweep);
    }
    fflush(NULL);
  }
  printf("\n");
  freeChain(&ch);

  clock_t endtime = clock();
  double timesecs = (endtime - starttime) / ((double)CLOCKS_PER_SEC);

  // --- Section 10 - Write statistics to files ---------
  write_stats_to_file(fname, seed, mode, adapt, doperm, nsweep2, nsweep, jd,
                      sig, st, timesecs);
  freeRunStats(st, jd);
  freeJD(jd);

  return EXIT_SUCCESS;
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
