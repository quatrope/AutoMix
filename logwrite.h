#include "automix.h"

void write_stats_to_file(char *fname, chainState ch, unsigned long seed,
                         int mode, int nsweep2, int nsweep, proposalDist jd,
                         runStats st);
void report_cond_prob_estimation(char *fname, int mode, proposalDist jd,
                                 condProbStats cpstats);
