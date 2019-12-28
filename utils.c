#include "utils.h"
#include <math.h>
#define DOUB
/* Combined congruential and Tauseworthe generators from SuperDuper
 * package. Should work on machines with unsigned long of at least 32
 * bits. JC and JT must be initialized to values with 0 < JC < 2^32 and
 * 0 < JT < 2^32. JC must be odd.
 * References: Marsaglia, Ananthanarayanan & Paul, 1973,
 * Learmonth & Lewis, 1973, and Dudewicz, 1976)
 */

/* Compilation flags:
   -DSUBR sdrand is used as a subroutine sdrand(u)
          instead of a function u = sdrand()
   -DSUNF generate Sun Fortran compatible entry names
          (with trailing _ symbol)
   -DDOUB sdrand or its argument are double precision
   -DBSD  solves problem with Sun Fortran, pre-Solaris versions
          (i.e. BSD), with single precision function option.
          Do not use this option for Solaris.
   -DRETS returns effective seed in argument to sdrni;
          if used, it is essential to use a variable, not
          a constant, in calling sdrni.
   -DLOG  prints seed value in file "rnilog".

   Examples:
        cc -c -o sd.o -DSUNF -DRETS sd.c
   (single precision Sun Fortran function, returning seed value from sdrni)

        cc -c -o sdc.o sd.c
   (single precision C function)
*/

#include <stdio.h>
#include <time.h>

#ifdef SUNF
#ifndef DOUB
#include <math.h>
#endif
#endif

static unsigned long JC, JT;
static double Norm = 4.656612873E-10;

#ifdef SUNF
#ifdef SUBR
void sdrand_(u)
#else
#ifdef DOUB
double sdrand_()
#else
#ifdef BSD
FLOATFUNCTIONTYPE sdrand_()
#else
float sdrand_()
#endif
#endif
#endif
#else
#ifdef SUBR
void sdrand(u)
#else
#ifdef DOUB
double sdrand()
#else
float sdrand()
#endif
#endif
#endif

#ifdef SUBR
#ifdef DOUB
    double *u;
#else
    float *u;
#endif
#endif

{
  JC = (JC * 69069) & 037777777777; /* congruential part */
  JT ^= JT >> 15;                   /* tausworthe part */
  JT ^= (JT << 17) & 037777777777;

#ifdef SUBR
#ifdef DOUB
  *u = ((JT ^ JC) >> 1) * Norm;
#else
  *u = (float)(((JT ^ JC) >> 1) * Norm);
#endif
#else
#ifdef DOUB
  return (((JT ^ JC) >> 1) * Norm);
#else
#ifdef SUNF
#ifdef BSD
  RETURNFLOAT((float)((JT ^ JC) >> 1) * Norm);
#else
  return ((float)((JT ^ JC) >> 1) * Norm);
#endif
#else
  return ((float)((JT ^ JC) >> 1) * Norm);
#endif
#endif
#endif
}

#ifdef SUNF
void sdrni_(unsigned long *i);
#else
void sdrni(unsigned long *i);
#endif

#ifdef SUNF
void sdrni_(i)
#else
void sdrni(i)
#endif

    unsigned long *i;
{
#ifdef LOG
  FILE *stream;
#endif
  unsigned long k = *i;
  if (k == 0)
    k = time(0);
  JT = k / 65536;
  JC = k - 65536 * JT;
  JT = 65536 * JT + 1;
  JC = 32768 * JC + 1;
#ifdef LOG
  stream = fopen("rnilog", "a+");
  fprintf(stream, "%12d\n", k);
  fclose(stream);
#endif
#ifdef RETS
  *i = k;
#endif
}

#ifdef SUNF
void sdpseed_()
#else
void sdpseed()
#endif
{
  printf("%lu %lu \n", JC, JT);
}

#ifdef SUNF
void sdset_(i, j)
#else
void sdset(i, j)
#endif
    unsigned long *i,
    *j;
{
  JC = *i;
  JT = *j;
}

#ifdef SUNF
void sdget_(i, j)
#else
void sdget(i, j)
#endif
    unsigned long *i,
    *j;
{
  *i = JC;
  *j = JT;
}
/* Functions to estimates integrated autocorrelation time using
   method of Sokal. Taken from PJG function sokal.f
   Note that the definition is the sum from minus infinity to
   infinity of the autocorrelation function, hence twice Sokal's
   definition. */

#include <math.h>
#include <stdio.h>

#define pi 3.141592653589793

void fastfr(int nin, double *xreal, double *ximag);

void sokal(int n, double *xreal, double *var, double *tau, int *m) {

  int i1;
  double ximag[n], c, sum;

  if (n > pow(2.0, 20.0)) {
    printf("\nAuto-correlation length exceeded");
    return;
  }

  for (i1 = 0; i1 < n; i1++) {
    ximag[i1] = 0.0;
  }

  /* Use FFT to compute autocorrelations (in xr)  and variance (in var). */

  fastfr(n, xreal, ximag);

  for (i1 = 0; i1 < n; i1++) {
    xreal[i1] = pow(xreal[i1], 2.0) + pow(ximag[i1], 2.0);
    ximag[i1] = 0.0;
  }

  xreal[0] = 0.0;

  fastfr(n, xreal, ximag);
  *var = xreal[0] / ((double)n * (n - 1));
  c = 1.0 / xreal[0];

  for (i1 = 0; i1 < n; i1++) {
    xreal[i1] = xreal[i1] * c;
  }

  /* Use Sokal's adaptive truncated periodogram method to estimate
     integrated autocorrelation time (in tau). */

  sum = -(1.0 / 3.0);
  for (i1 = 0; i1 < n; i1++) {
    sum += xreal[i1] - (1.0 / 6.0);
    if (sum < 0.0) {
      goto CALC_TAU;
    }
  }

CALC_TAU:
  *tau = 2.0 * (sum + i1 / 6.0);
  *m = i1 + 1;

  return;
}

void fastfr(int nin, double *xreal, double *ximag) {

  /* radix 4 complex discrete fast Fourier transform
     without usual normalisation
     Eric Renshaw -> PJG 20 March 1987 -> DIH, March 2004

     xreal = array which on input contains real part of data
     for transformation and on output gives real part of the
     result,type real, dimension iabs(isize) or greater.
     ximag = array which on input contains imaginary part
     for transformation and on output gives imaginary part of
     the result, type real, dimension iabs(isize) or greater.
     isize = integer variable or constant specifying length and
     type of transform. The length is iabs(isize) which must be
     a power of 2, minimum length 4, maximum 2**20. Illegal length
     leads to warning and return. If isize is positive the forward
     transform is calculated. If negative the inverse transform is found.

     The transform is defined by
     out(r) = sum(j = 0 to n-1) in(j)*exp(-2*pi*i*r*j/n)
                        if isize = n > 0,
     out(r) = sum(j = 0 to n-1) in(j)*exp(2*pi*i*r*j/n)
                        if isize = -n < 0,
     for r = 0,1,2,...,(n-1),
       where i = sqrt(-1), and both in(j) and out(j)
       are stored in xreal(j+1)+i*ximag(j+1)              */

  double z, bcos, bsin, temp, cw1, sw1;
  double cw2 = 0.0;
  double sw2 = 0.0;
  double cw3 = 0.0;
  double sw3 = 0.0;
  double xs0, xs1, xs2, xs3, ys0, ys1, ys2, ys3, x1, x2, x3, y1, y2, y3;

  int i0, i1, i2, i3, ul[20], n, counta, countb, time, test;
  int indic = 0;
  int j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16,
      j17, j18, j19;

  n = nin < 0 ? -nin : nin;

  if (n < 4) {
    printf("\nFast Fourier length too short");
    return;
  }

  /* if this is to be an inverse transform, conjugate the data */

  if (nin < 0) {
    for (i1 = 0; i1 < n; i1++) {
      ximag[i1] = -ximag[i1];
    }
  }

  test = n;
  while (test > 1) {
    if (test & 1) {
      printf("\nFast Fourier length must be power of 2");
      return;
    }
    test >>= 1;
  }

  counta = n / 4;
  time = 0;

  while (counta > 0) {
    countb = counta * 4;
    time += 2;

    /* do the transforms required by this stage */

    z = pi / countb;
    bcos = -2.0 * pow(sin(z), 2.0);
    bsin = sin(2.0 * z);

    cw1 = 1.0;
    sw1 = 0.0;

    /* this is the main calculation of radix 4 transforms */

    for (j1 = 0; j1 < counta; j1++) {
      for (j2 = j1; j2 < n; j2 += countb) {
        i0 = j2;
        i1 = i0 + counta;
        i2 = i1 + counta;
        i3 = i2 + counta;
        xs0 = xreal[i0] + xreal[i2];
        xs1 = xreal[i0] - xreal[i2];
        ys0 = ximag[i0] + ximag[i2];
        ys1 = ximag[i0] - ximag[i2];
        xs2 = xreal[i1] + xreal[i3];
        xs3 = xreal[i1] - xreal[i3];
        ys2 = ximag[i1] + ximag[i3];
        ys3 = ximag[i1] - ximag[i3];

        xreal[i0] = xs0 + xs2;
        ximag[i0] = ys0 + ys2;

        x1 = xs1 + ys3;
        y1 = ys1 - xs3;
        x2 = xs0 - xs2;
        y2 = ys0 - ys2;
        x3 = xs1 - ys3;
        y3 = ys1 + xs3;

        if (j1 == 0) {
          xreal[i2] = x1;
          ximag[i2] = y1;
          xreal[i1] = x2;
          ximag[i1] = y2;
          xreal[i3] = x3;
          ximag[i3] = y3;
        } else {

          /* multiply by twiddle factors if required */

          xreal[i2] = x1 * cw1 + y1 * sw1;
          ximag[i2] = y1 * cw1 - x1 * sw1;
          xreal[i1] = x2 * cw2 + y2 * sw2;
          ximag[i1] = y2 * cw2 - x2 * sw2;
          xreal[i3] = x3 * cw3 + y3 * sw3;
          ximag[i3] = y3 * cw3 - x3 * sw3;
        }
      }
      if (j1 < (counta - 1)) {

        /* calculate a new set of twiddle factors */

        z = cw1 * bcos - sw1 * bsin + cw1;
        sw1 = bcos * sw1 + bsin * cw1 + sw1;
        temp = 1.5 - 0.5 * (z * z + sw1 * sw1);
        cw1 = z * temp;
        sw1 = sw1 * temp;
        cw2 = cw1 * cw1 - sw1 * sw1;
        sw2 = 2.0 * cw1 * sw1;
        cw3 = cw1 * cw2 - sw1 * sw2;
        sw3 = cw1 * sw2 + cw2 * sw1;
      }
    }

    indic = 0;
    if (counta > 1) {
      /* set up the transform split for the next stage */
      counta /= 4;
      indic = 1;
    } else {
      counta = 0;
    }
  }
  if (indic) {

    /* this is the calculation of a radix two stage */
    for (j1 = 0; j1 < n; j1 += 2) {
      temp = xreal[j1] + xreal[j1 + 1];
      xreal[j1 + 1] = xreal[j1] - xreal[j1 + 1];
      xreal[j1] = temp;
      temp = ximag[j1] + ximag[j1 + 1];
      ximag[j1 + 1] = ximag[j1] - ximag[j1 + 1];
      ximag[j1] = temp;
    }
    time++;
  }

  /* if this was an inverse transform, conjugate the result */

  if (nin < 0) {
    for (j1 = 0; j1 < n; j1++) {
      ximag[j1] = -ximag[j1];
    }
  }

  /* unscramble the result */

  if (time > 20) {
    printf("\nFast Fourier length too long");
    return;
  }

  i1 = 20 - time;
  for (j1 = 0; j1 < i1; j1++) {
    ul[j1] = 1;
  }
  if (i1 == 0) {
    ul[0] = 2;
    i1++;
  }
  for (j1 = i1; j1 < 20; j1++) {
    ul[j1] = 2 * ul[j1 - 1];
  }

  i0 = 0;
  for (j0 = 0; j0 < ul[0]; j0++) {
    for (j1 = j0; j1 < ul[1]; j1 += ul[0]) {
      for (j2 = j1; j2 < ul[2]; j2 += ul[1]) {
        for (j3 = j2; j3 < ul[3]; j3 += ul[2]) {
          for (j4 = j3; j4 < ul[4]; j4 += ul[3]) {
            for (j5 = j4; j5 < ul[5]; j5 += ul[4]) {
              for (j6 = j5; j6 < ul[6]; j6 += ul[5]) {
                for (j7 = j6; j7 < ul[7]; j7 += ul[6]) {
                  for (j8 = j7; j8 < ul[8]; j8 += ul[7]) {
                    for (j9 = j8; j9 < ul[9]; j9 += ul[8]) {
                      for (j10 = j9; j10 < ul[10]; j10 += ul[9]) {
                        for (j11 = j10; j11 < ul[11]; j11 += ul[10]) {
                          for (j12 = j11; j12 < ul[12]; j12 += ul[11]) {
                            for (j13 = j12; j13 < ul[13]; j13 += ul[12]) {
                              for (j14 = j13; j14 < ul[14]; j14 += ul[13]) {
                                for (j15 = j14; j15 < ul[15]; j15 += ul[14]) {
                                  for (j16 = j15; j16 < ul[16]; j16 += ul[15]) {
                                    for (j17 = j16; j17 < ul[17];
                                         j17 += ul[16]) {
                                      for (j18 = j17; j18 < ul[18];
                                           j18 += ul[17]) {
                                        for (j19 = j18; j19 < ul[19];
                                             j19 += ul[18]) {
                                          if ((i0 - j19) < 0) {
                                            temp = xreal[i0];
                                            xreal[i0] = xreal[j19];
                                            xreal[j19] = temp;
                                            temp = ximag[i0];
                                            ximag[i0] = ximag[j19];
                                            ximag[j19] = temp;
                                          }
                                          i0++;
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return;
}
/* Routines calculating quantities related to the gamma function */

#include <math.h>

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

/* Taken from algama.f (PJG) - converted to C using appropriate machine
   constants by DIH 04/11/03 */

double loggamma(double X) {

  /* This routine calculates the LOG GAMMA function for a positive real
     argument X.  Computation is based on an algorithm outlined in
     references 1 and 2.  The program uses rational functions that
     theoretically approximate LOG GAMMA to at least 18 significant
     decimal digits.  The approximation for X > 12 is from reference
     3, while approximations for X < 12.0 are similar to those in
     reference 1, but are unpublished.  The accuracy achieved depends
     on the arithmetic system, the compiler, the intrinsic functions,
     and proper selection of the machine-dependent constants.

     Values taken from float.h and algama.f (for XBIG)

     ---------------------------------------------------------------

     Explanation of machine-dependent constants

  beta   - radix for the floating-point representation
  maxexp - the smallest positive power of beta that overflows
  XBIG   - largest argument for which LN(GAMMA(X)) is representable
           in the machine, i.e., the solution to the equation
                   LN(GAMMA(XBIG)) = beta**maxexp
  XINF   - largest machine representable floating-point number;
           approximately beta**maxexp.
  EPS    - The smallest positive floating-point number such that
           1.0+EPS > 1.0
  FRTBIG - Rough estimate of the fourth root of XBIG

     ---------------------------------------------------------------

     Error returns

     The program returns the value XINF for X <= 0.0 or when
     overflow would occur.  The computation is believed to
     be free of underflow and overflow.

     ---------------------------------------------------------------
     References:

     1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
     1967, pp. 198-203.

     2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
     1969.

     3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
     York, 1968.

  -----------------------------------------------------------------
  Start of code
  -----------------------------------------------------------------*/

  int I;
  double CORR, D1, D2, D4, EPS, FRTBIG, FOUR, HALF, ONE, PNT68;
  double RES, SQRTPI, THRHAL, TWELVE, TWO, XBIG, XDEN, XINF;
  double XM1, XM2, XM4, XNUM, Y, YSQ, ZERO;
  double C[7], P1[8], P2[8], P4[8], Q1[8], Q2[8], Q4[8];

  /*---------------------------------------------------------------------
  Mathematical constants
  ---------------------------------------------------------------------*/

  ONE = 1.0E0;
  HALF = 0.5E0;
  TWELVE = 12.0E0;
  ZERO = 0.0E0;
  FOUR = 4.0E0;
  THRHAL = 1.5E0;
  TWO = 2.0E0;
  PNT68 = 0.6796875E0;
  SQRTPI = 0.9189385332046727417803297E0; /* eh? */

  /*---------------------------------------------------------------------
    Machine dependent parameters
    -------------------------------------------------------------------*/

  XBIG = 2.55E305;
  XINF = 1.79E308;
  EPS = 2.22E-16;
  FRTBIG = 2.25E76;

  /*--------------------------------------------------------------------
    Numerator and denominator coefficients for rational minimax
    approximation over (EPS,1.5).
    -------------------------------------------------------------------*/
  D1 = -5.772156649015328605195174E-1;
  P1[0] = 4.945235359296727046734888E0;
  P1[1] = 2.018112620856775083915565E2;
  P1[2] = 2.290838373831346393026739E3;
  P1[3] = 1.131967205903380828685045E4;
  P1[4] = 2.855724635671635335736389E4;
  P1[5] = 3.848496228443793359990269E4;
  P1[6] = 2.637748787624195437963534E4;
  P1[7] = 7.225813979700288197698961E3;
  Q1[0] = 6.748212550303777196073036E1;
  Q1[1] = 1.113332393857199323513008E3;
  Q1[2] = 7.738757056935398733233834E3;
  Q1[3] = 2.763987074403340708898585E4;
  Q1[4] = 5.499310206226157329794414E4;
  Q1[5] = 6.161122180066002127833352E4;
  Q1[6] = 3.635127591501940507276287E4;
  Q1[7] = 8.785536302431013170870835E3;

  /*---------------------------------------------------------------------
     Numerator and denominator coefficients for rational minimax
     Approximation over (1.5,4.0).
     ------------------------------------------------------------------*/

  D2 = 4.227843350984671393993777E-1;
  P2[0] = 4.974607845568932035012064E0;
  P2[1] = 5.424138599891070494101986E2;
  P2[2] = 1.550693864978364947665077E4;
  P2[3] = 1.847932904445632425417223E5;
  P2[4] = 1.088204769468828767498470E6;
  P2[5] = 3.338152967987029735917223E6;
  P2[6] = 5.106661678927352456275255E6;
  P2[7] = 3.074109054850539556250927E6;
  Q2[0] = 1.830328399370592604055942E2;
  Q2[1] = 7.765049321445005871323047E3;
  Q2[2] = 1.331903827966074194402448E5;
  Q2[3] = 1.136705821321969608938755E6;
  Q2[4] = 5.267964117437946917577538E6;
  Q2[5] = 1.346701454311101692290052E7;
  Q2[6] = 1.782736530353274213975932E7;
  Q2[7] = 9.533095591844353613395747E6;

  /*--------------------------------------------------------------------
    Numerator and denominator coefficients for rational minimax
    Approximation over (4.0,12.0).
    -------------------------------------------------------------------*/

  D4 = 1.791759469228055000094023E0;
  P4[0] = 1.474502166059939948905062E4;
  P4[1] = 2.426813369486704502836312E6;
  P4[2] = 1.214755574045093227939592E8;
  P4[3] = 2.663432449630976949898078E9;
  P4[4] = 2.940378956634553899906876E10;
  P4[5] = 1.702665737765398868392998E11;
  P4[6] = 4.926125793377430887588120E11;
  P4[7] = 5.606251856223951465078242E11;
  Q4[0] = 2.690530175870899333379843E3;
  Q4[1] = 6.393885654300092398984238E5;
  Q4[2] = 4.135599930241388052042842E7;
  Q4[3] = 1.120872109616147941376570E9;
  Q4[4] = 1.488613728678813811542398E10;
  Q4[5] = 1.016803586272438228077304E11;
  Q4[6] = 3.417476345507377132798597E11;
  Q4[7] = 4.463158187419713286462081E11;

  /*---------------------------------------------------------------------
    Coefficients for minimax approximation over (12, INF).
    -------------------------------------------------------------------*/
  C[0] = -1.910444077728E-03;
  C[1] = 8.4171387781295E-04;
  C[2] = -5.952379913043012E-04;
  C[3] = 7.93650793500350248E-04;
  C[4] = -2.777777777777681622553E-03;
  C[5] = 8.333333333333333331554247E-02;
  C[6] = 5.7083835261E-03;

  /*----------------------------------------------------------------------
    0 < X <= EPS
    --------------------------------------------------------------------*/
  Y = X;
  if ((Y > 0) && (Y <= XBIG)) {
    if (Y <= EPS) {
      RES = -log(Y);
    } else if (Y <= THRHAL) {
      /*-----------------------------------------------------------------------
        EPS < X <= 1.5
        ---------------------------------------------------------------------*/
      if (Y < PNT68) {
        CORR = -log(Y);
        XM1 = Y;
      } else {
        CORR = ZERO;
        XM1 = (Y - HALF) - HALF;
      }

      if ((Y <= HALF) || (Y >= PNT68)) {
        XDEN = ONE;
        XNUM = ZERO;
        for (I = 0; I < 8; I++) {
          XNUM = XNUM * XM1 + P1[I];
          XDEN = XDEN * XM1 + Q1[I];
        }
        RES = CORR + (XM1 * (D1 + XM1 * (XNUM / XDEN)));
      } else {
        XM2 = (Y - HALF) -
              HALF; /*Is XM2 symbol used to agree with other 2 symbols*/
        XDEN = ONE;
        XNUM = ZERO;
        for (I = 0; I < 8; I++) {
          XNUM = XNUM * XM2 + P2[I];
          XDEN = XDEN * XM2 + Q2[I];
        }
        RES = CORR + XM2 * (D2 + XM2 * (XNUM / XDEN));
      }
    } else if (Y <= FOUR) {

      /*---------------------------------------------------------------------
        1.5 < X <= 4.0
        -------------------------------------------------------------------*/
      XM2 = Y - TWO;
      XDEN = ONE;
      XNUM = ZERO;
      for (I = 0; I < 8; I++) {
        XNUM = XNUM * XM2 + P2[I];
        XDEN = XDEN * XM2 + Q2[I];
      }
      RES = XM2 * (D2 + XM2 * (XNUM / XDEN));
    } else if (Y <= TWELVE) {

      /*----------------------------------------------------------------------
        4.0 < X <= 12.0
        --------------------------------------------------------------------*/
      XM4 = Y - FOUR;
      XDEN = -ONE;
      XNUM = ZERO;
      for (I = 0; I < 8; I++) {
        XNUM = XNUM * XM4 + P4[I];
        XDEN = XDEN * XM4 + Q4[I];
      }
      RES = D4 + XM4 * (XNUM / XDEN);
    } else {

      /*---------------------------------------------------------------------
        X > 12.0,
        -------------------------------------------------------------------*/
      RES = ZERO;
      if (Y <= FRTBIG) {
        RES = C[6];
        YSQ = Y * Y;
        for (I = 0; I < 6; I++) {
          RES = RES / YSQ + C[I];
        }
      }
      RES = RES / Y;
      CORR = log(Y);
      RES = RES + SQRTPI - HALF * CORR;
      RES = RES + Y * (CORR - ONE);
    }
  } else {

    /*----------------------------------------------------------------------
      Return for bad arguments
      --------------------------------------------------------------------*/
    RES = XINF;
  }

  /*----------------------------------------------------------------------
    Final return
    --------------------------------------------------------------------*/
  return (RES);
}

/* Function for generating a Gamma(s,1) random variable
   Converted to C from PJG rgamma FORTRAN function
   Note: (1/t)*GAMMA(s,1) is GAMMA(s,t) */

double rgamma(double s) {

  double exp1, b, c1, c2, c3, c4, c5, u1, u2, w, bu, out;

  exp1 = exp(1.0);

  if (s < 1) {
    b = (s + exp1) / exp1;
    c1 = 1.0 / s;
  LAB_1:
    bu = b * sdrand();
    if (bu <= 1.0) {
      out = exp(max(-30.0, c1 * log(bu)));
      if (sdrand() >= exp(-out)) {
        goto LAB_1;
      }
    } else {
      out = -log((b - bu) / s);
      if (sdrand() >= pow(out, (s - 1.0))) {
        goto LAB_1;
      }
    }
  } else if (s == 1.0) {
    out = -log(sdrand());
  } else {
    c1 = s - 1.0;
    c2 = (s - 1.0 / (6.0 * s)) / c1;
    c3 = 2.0 / c1;
    c4 = c3 + 2.0;
    c5 = 1.0 / sqrt(s);
  LAB_2:
    u1 = sdrand();
    u2 = sdrand();
    if (s > 2.5) {
      u1 = u2 + c5 * (1.0 - 1.86 * u1);
    }
    if (u1 <= 0.0 || u1 >= 1.0) {
      goto LAB_2;
    }
    w = c2 * u2 / u1;
    if ((c3 * u1 + w + 1.0 / w) <= c4) {
      goto LAB_3;
    }
    if ((c3 * log(u1) - log(w) + w) >= 1.0) {
      goto LAB_2;
    }
  LAB_3:
    out = c1 * w;
  }

  return out;
}

void gauss(double *z, int n) {
  /* Uses Box mueller method to simulate n N(0,1) variables and stores them
     in z */

  int n1;
  double u, v;

  n1 = n - 1;
  for (int i = 0; i < n1; i += 2) {
    u = sqrt(-2.0 * log(sdrand()));
    v = 2.0 * M_PI * sdrand();
    z[i] = u * sin(v);
    z[i + 1] = u * cos(v);
  }
  if (fmod(n, 2) < 0.5) {
    return;
  } else {
    u = sqrt(-2.0 * log(sdrand()));
    v = 2.0 * M_PI * sdrand();
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

void chol(int n, double **A) {
  /* Performs cholesky decompositon of A and returns result in the
     same matrix - adapted from PJG Fortran function*/

  for (int j1 = 0; j1 < n; j1++) {
    double sum = A[j1][j1];
    for (int j2 = 0; j2 < j1; j2++) {
      sum -= pow(A[j1][j2], 2);
    }
    A[j1][j1] = sqrt(sum);

    for (int j2 = j1 + 1; j2 < n; j2++) {
      sum = A[j2][j1];
      for (int j3 = 0; j3 < j1; j3++) {
        sum -= A[j2][j3] * A[j1][j3];
      }
      A[j2][j1] = sum / A[j1][j1];
    }
  }
}

void perm(double *work, int n) {
  /* Randomly permutes the n-dim work vector */

  for (int i = 0; i < (n - 1); i++) {
    int j = i + (int)((n - i) * sdrand());
    if (j != i) {
      double temp = work[j];
      work[j] = work[i];
      work[i] = temp;
    }
  }
  return;
}

double ltprob(int dof, double z) {
  /* Evaluates the log of p.d.f. of a t variable with dof degrees of freedom
     at point z */

  double constt =
      loggamma(0.5 * (dof + 1)) - loggamma(0.5 * dof) - 0.5 * log(dof * M_PI);
  double out = constt - 0.5 * (dof + 1) * log(1.0 + pow(z, 2.0) / dof);
  return out;
}

double lnormprob(int n, int l, double **mu_k, double ***B_k, double *datai) {
  /* Evaluates log of p.d.f. for a multivariate normal for model
     k, of dimension n, component l. The summary of means and
     sqrt of cov matrices (for all models and all component)
     are supplied in mu and B */

  double work[n];

  for (int i = 0; i < n; i++) {
    work[i] = datai[i] - mu_k[l][i];
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      (work[i]) -= B_k[l][i][j] * work[j];
    }
    (work[i]) /= B_k[l][i][i];
  }
  double out = 0.0;
  for (int i = 0; i < n; i++) {
    out += (work[i] * work[i]);
  }
  out = -0.5 * out - (n / 2.0) * log(2.0 * M_PI) - log(det(n, l, B_k));
  return out;
}

double det(int n, int l, double ***B_k) {

  /* Evaluates the determinant of a matrix in B corresponding to model k,
     component l. */
  double out = 1.0;
  for (int i = 0; i < n; i++) {
    out *= B_k[l][i][i];
  }
  return out;
}
