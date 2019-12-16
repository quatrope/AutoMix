/* --- User supplied functions ------------ */

/* Functions must be supplied in user***.c file (see e.g. usertoy1.c,
   usercpt.c etc bundled with this software).

   Descriptions:
   1. lpost(&k,theta,&llh)
   This should be a c function written by the user that evaluates
   log posterior (up to an additive constant) at (k,theta). The function
   can also return the likelihood at this point in llh.
   2. getkmax(&kmax)
   This should be a c function written by the user that returns the
   number of models kmax.
   3. getnk(kmax, nk)
   This should be a c function written by the user that returns the dimensions
   nk for model k=1,...,kmax.
   4. getic(k,nkk,rwm)
   This should be a c function written by the user that returns the
   possibly random starting point for the rwm in stage 1 of the AutoMix
   sampler */

double lpost(int k, int nkk, double *theta, double *llh);
void getkmax(int *kmax);
void getnk(int kmax, int nk[kmax]);
void getic(int k, int nkk, double *rwm);
