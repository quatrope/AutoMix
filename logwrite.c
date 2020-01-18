#include "logwrite.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void write_ac_to_file(char *fname, int m, double *xr);
void write_log_to_file(char *fname, amSampler am, unsigned long seed, int mode,
                       chainState ch, int nsweep2, int nsweep, proposalDist jd,
                       runStats st);
void write_pk_to_file(char *fname, int nsweep, int nmodels,
                      double **pk_summary);
void write_k_to_file(char *fname, int nsweep, int *k_which_summary);
void write_lp_to_file(char *fname, int nsweep, double **logp_summary);
void write_cf_to_file(char *fname, amSampler am, condProbStats cpstats);
void write_adapt_to_file(char *fname, proposalDist jd, condProbStats cpstats);
void write_mix_to_file(char *fname, proposalDist jd);
void write_theta_to_file(char *fname, runStats st, proposalDist jd);

void report_cond_prob_estimation(char *fname, amSampler am) {
  // Write adaptation statistics to file
  write_adapt_to_file(fname, am.jd, am.cpstats);
  // Write mixture parameters to file
  write_mix_to_file(fname, am.jd);
  // Write cf statistics to file
  write_cf_to_file(fname, am, am.cpstats);
}

void write_cf_to_file(char *fname, amSampler am, condProbStats cpstats) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_cf.data", fname);
  FILE *fp_cf = fopen(datafname, "w");
  free(datafname);
  int nmodels = am.jd.nmodels;
  if (am.am_mixfit == FIGUEREIDO_MIX_FIT) {
    for (int model_k = 0; model_k < nmodels; model_k++) {
      fprintf(fp_cf, "RWM for Model %d\n", model_k + 1);
      for (int i = 0; i < cpstats.nfitmix[model_k]; i++) {
        fprintf(fp_cf, "%d %lf %lf %d\n", cpstats.fitmix_Lkk[model_k][i],
                cpstats.fitmix_lpn[model_k][i],
                cpstats.fitmix_costfnnew[model_k][i],
                cpstats.fitmix_annulations[model_k][i]);
      }
      fflush(NULL);
    }
  }
  fclose(fp_cf);
}

void write_adapt_to_file(char *fname, proposalDist jd, condProbStats cpstats) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_adapt.data", fname);
  FILE *fp_adapt = fopen(datafname, "w");
  free(datafname);
  for (int model_k = 0; model_k < jd.nmodels; model_k++) {
    int mdim = jd.model_dims[model_k];
    fprintf(fp_adapt, "RWM for Model %d\n", model_k + 1);
    for (int j = 0; j < cpstats.rwm_summary_len; j++) {
      for (int i = 0; i < mdim; i++) {
        fprintf(fp_adapt, "%lf %lf ", cpstats.sig_k_rwm_summary[model_k][j][i],
                cpstats.nacc_ntry_rwm[model_k][j][i]);
      }
      fprintf(fp_adapt, "\n");
    }
  }
  fclose(fp_adapt);
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

void write_theta_to_file(char *fname, runStats st, proposalDist jd) {
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  for (int model_k = 0; model_k < jd.nmodels; model_k++) {
    sprintf(datafname, "%s_theta%d.data", fname, model_k + 1);
    FILE *fp_theta = fopen(datafname, "a");
    int sample_len = st.theta_summary_len[model_k];
    int mdim = jd.model_dims[model_k];
    for (int sample_i = 0; sample_i < sample_len; sample_i++) {
      double *theta = st.theta_summary[model_k][sample_i];
      for (int j = 0; j < mdim; j++) {
        fprintf(fp_theta, "%lf ", theta[j]);
      }
      fprintf(fp_theta, "\n");
    }
    fclose(fp_theta);
  }
  free(datafname);
}

void write_stats_to_file(char *fname, amSampler am, int mode, int nsweep2,
                         int nsweep) {
  runStats st = am.st;
  write_pk_to_file(fname, nsweep, (am.jd).nmodels, st.pk_summary);
  write_k_to_file(fname, nsweep, st.k_which_summary);
  write_lp_to_file(fname, nsweep, st.logp_summary);
  sokal(st.nkeep, st.xr, &(st.var), &(st.tau), &(st.m));
  write_log_to_file(fname, am, am.seed, mode, am.ch, nsweep2, nsweep, am.jd,
                    st);
  write_ac_to_file(fname, st.m, st.xr);
  write_theta_to_file(fname, st, am.jd);
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

void write_mix_to_file(char *fname, proposalDist jd) {
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
      fprintf(fp_mix, "%lf\n", jd.sig[k1][j1]);
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

void write_log_to_file(char *fname, amSampler am, unsigned long seed, int mode,
                       chainState ch, int nsweep2, int nsweep, proposalDist jd,
                       runStats st) {

  // Print user options to log file
  unsigned long fname_len = strlen(fname);
  char *datafname = (char *)malloc((fname_len + 50) * sizeof(*datafname));
  sprintf(datafname, "%s_log.data", fname);
  FILE *fp_log = fopen(datafname, "w");
  free(datafname);
  fprintf(fp_log, "seed: %ld\n", seed);
  fprintf(fp_log, "m: %d\n", mode);
  fprintf(fp_log, "a: %d\n", am.doAdapt);
  fprintf(fp_log, "p: %d\n", am.doPerm);
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
      fprintf(fp_log, "%lf ", jd.sig[k1][j1]);
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
  fprintf(fp_log, "nkeep:%d, nsokal:%d, var:%lf, tau:%lf\n", st.nkeep,
          st.nsokal, st.var, st.tau);
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
  double timesecs = st.timesecs_burn + st.timesecs_rjmcmc;
  fprintf(fp_log, "Time: %lf\n", timesecs);
  fclose(fp_log);
}
