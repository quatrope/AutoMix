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

/* Routines calculating quantities related to the gamma function */
double rgamma(double s);
double loggamma(double s);

/* Functions to estimates integrated autocorrelation time using
 method of Sokal. Taken from PJG function sokal.f
 Note that the definition is the sum from minus infinity to
 infinity of the autocorrelation function, hence twice Sokal's
 definition. */
void sokal(int n, double *xreal, double *var, double *tau, int *m);
