#include "automix.h"

void write_stats_to_file(char *fname, chainState ch, unsigned long seed, int mode, int nsweep2, int nsweep, proposalDist jd,
                         double **sig, runStats st, double timesecs);
void write_cf_to_file(char *fname, int mode, proposalDist jd, runStats st);
void write_adapt_to_file(char *fname, int mode, proposalDist jd, runStats st);
void write_mix_to_file(char *fname, proposalDist jd, double **sig);
void write_theta_to_file(char *fname, int current_model_k, int mdim,
                         double *theta);
