/* --- User supplied functions ------------ */

/* Functions must be supplied in user***.c file (see e.g. usertoy1.c,
   usercpt.c etc bundled with this software).
 */

// Evaluates log posterior (up to an additive constant) at (k,theta).
// The function can also return the likelihood at this point in llh.
double logpost(int model_k, int mdim, double *theta, double *llh);

// Returns the number of models
int get_nmodels(void);

// Loads the dimensions of each model in model_dims.
void load_model_dims(int nmodels, int *model_dims);

// Returns the possibly random starting point for the rwm in stage 1 of the
// AutoMix sampler
void get_rwm_init(int model_k, int mdim, double *rwm);
