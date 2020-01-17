// Driver main function for the AutoMix RJMCMC sampler
#define VERSION "1.3"

#include "automix.h"
#include "logwrite.h"
#include "user.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define max(A, B) ((A) > (B) ? (A) : (B))

void usage(char *invocation);
int parse_cmdline_args(int argc, char *argv[], char **fname, int *nsweep,
                       int *nsweep2, unsigned long *seed, int *doperm,
                       int *adapt, int *mode, int *dof);
// wrapper function to fit with AutoMix
double logposterior(int model_k, int mdim, double *x) {
  double lpost = 0;
  double likelihood;
  logpost(model_k, mdim, x, &lpost, &likelihood);
  return lpost;
}

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
  int parsing = parse_cmdline_args(argc, argv, &fname, &nsweep, &nsweep2, &seed,
                                   &doperm, &adapt, &mode, &dof);
  if (parsing == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }
  sdrni(&seed);
  int nburn = max(10000, (int)(nsweep / 10));

  // Initialize the AutoMix sampler
  int nmodels = get_nmodels();
  int *model_dims = (int *)malloc(nmodels * sizeof(int));
  load_model_dims(nmodels, model_dims);
  amSampler am;
  initAMSampler(&am, nmodels, model_dims, logposterior, get_rwm_init);

  // --- Section 5.1 - Read in mixture parameters if mode 1 (m=1) ---
  if (mode == 1) {
    // Read AutoMix parameters from file if mode = 1
    int ok = read_mixture_params(fname, &am);
    if (ok == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
  } else {
    initCondProbStats(&(am.cpstats), am.jd, nsweep2);
    estimate_conditional_probs(&am, nsweep2);
    report_cond_prob_estimation(fname, am);
    freeCondProbStats(am.cpstats, am.jd);
  }

  // Initialization of the MC Markov Chain parameters
  initChain(&(am.ch), am.jd, adapt, logposterior, get_rwm_init);
  // Struct to hold run statistic variables
  runStats st;
  initRunStats(&st, nsweep, am.jd);

  // -----Start of main loop ----------------
  // Burn some samples first
  burn_samples(&am, nburn, &st);
  // Collect nsweep RJMCMC samples
  rjmcmc_samples(&am, nsweep, &st);
  // --- Section 10 - Write statistics to files ---------
  write_stats_to_file(fname, am.ch, seed, mode, nsweep2, nsweep, am.jd, st);

  freeChain(&(am.ch));
  freeRunStats(st, am.jd);

  clock_t endtime = clock();
  double timesecs = (endtime - starttime) / ((double)CLOCKS_PER_SEC);
  printf("Total time elapsed: %f sec.\n", timesecs);

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

int parse_cmdline_args(int argc, char *argv[], char **fname, int *nsweep,
                       int *nsweep2, unsigned long *seed, int *doperm,
                       int *adapt, int *mode, int *dof) {
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-f")) {
      *fname = argv[++i];
      continue;
    } else if (!strcmp(argv[i], "-N")) {
      *nsweep = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-n")) {
      *nsweep2 = atoi(argv[++i]);
      *nsweep2 = max(*nsweep2, 1E5);
      continue;
    } else if (!strcmp(argv[i], "-s")) {
      *seed = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-p")) {
      *doperm = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-m")) {
      *mode = atoi(argv[++i]);
      if (*mode > 2 || *mode < 0) {
        printf("Error: Invalid mode entered. Mode must be {0, 1, 2}.\n");
        usage(argv[0]);
        return EXIT_FAILURE;
      }
      continue;
    } else if (!strcmp(argv[i], "-a")) {
      *adapt = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "-t")) {
      *dof = atoi(argv[++i]);
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
  return EXIT_SUCCESS;
}
