/* --- User supplied functions ------------ */

/* Functions must be supplied in user***.c file (see e.g. usertoy1.c,
   usercpt.c etc bundled with this software).
 */

// Evaluates log posterior (up to an additive constant) at (k,theta).
// The function can also return the likelihood at this point in llh.
double lpost(int k, int nkk, double *theta, double *llh);

// Returns the number of models
void getkmax(int *kmax);

// Returns the dimensions nk for model k=1,...,kmax.
void getnk(int kmax, int nk[kmax]);

// Returns the possibly random starting point for the rwm in stage 1 of the
// AutoMix sampler
void getic(int k, int nkk, double *rwm);
