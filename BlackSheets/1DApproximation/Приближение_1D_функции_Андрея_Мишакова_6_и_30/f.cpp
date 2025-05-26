#include <math.h>
#include "f.h"

#define EPS (1e-15)

double f_0 (double x)
{
  (void) x;
  return 1;
}

double f_1 (double x)
{
  return x;
}

double f_2 (double x)
{
  return x * x;
}

double f_3 (double x)
{
  return x * x * x;
}

double f_4 (double x)
{
  return x * x * x * x;
}

double sin_my (double x0, double eps)
  {
    double res = 0, pr;
    int zn = 1, i;
    if (x0 < 0) // [0, inf)
      {
        zn *= -1;
        x0 = -x0;
      }
    x0 = fmod (x0, 2*M_PI); // [0, 2*pi)
    if (x0 > M_PI) // [0, pi)
      {
        zn *= -1;
        x0 -= M_PI;
      }
    if (x0 > M_PI/2) /* [0, pi/2) */
      x0 = M_PI - x0;
    if (x0 > M_PI/4) /* [pi/4, pi/2) */
      res = cos_my (M_PI/2-x0, eps);
    else /* [0, pi/4) */
      {
        eps = fmax(eps, EPS);
        pr = x0; /* x/1! */
        for (i = 1; fabs (pr) >= eps; i++)
          { 
            res += pr;
            pr *= - x0 / (2 * i) * x0 / (2 * i + 1);
            /*printf("%lf %lf\n", pr, res);*/
          }
      }
    /*printf("%f\n", x0);*/
    if (zn == 1)
      return res;
    else
      return -res;
  }
double cos_my (double x0, double eps)
  {
    double res = 0, pr;
    int zn = 1, i;
    if (x0 < 0) /* [0, inf) */
      x0 = -x0;
    x0 = fmod(x0, 2 * M_PI); /* [0, 2*pi) */
    if (x0 > M_PI) /* [0, pi) */
      {
        zn *= -1;
        x0 -= M_PI;
      }
    if (x0 > M_PI / 2) /* [0, pi/2) */
      {
        zn *= -1;
        x0 = M_PI - x0;
      }
    if (x0 > M_PI / 4) /* [pi/4, pi/2) */
      res = sin_my (M_PI / 2 - x0, eps);
    else /* [0, pi/4) */
      {
        eps = fmax(eps, EPS);
        pr = 1; /* x^0/0! */
        for (i = 1; fabs(pr) >= eps; i++)
          { 
            res += pr;
            pr *= - x0 /(2 * i) * x0 / (2 * i - 1);
            /*printf("%lf %lf\n", pr, res);*/
          }
      }
    /*printf("%f\n", x0);*/
    if (zn == 1)
      return res;
    else
      return -res;
  }

double exp_my (double x0, double eps)
  {
    double res_dr, res_z = 0, pr, x_z, x_dr;
    int zn = 1, i;
    if (x0 < 0) { /* [0, inf) */
      zn *= -1;
      x0 = -x0;
    }
    x_z = floor(x0); 
    x_dr =  x0 - x_z; /* [0, 1) */
    
    res_z = pow(M_E, x_z);
    
    res_dr = 0;
    eps = fmax(eps, EPS);
    pr = 1; /* x^1/0! */
    for (i = 1; fabs(pr) >= eps; i++)
      { 
        res_dr += pr;
        pr *= x_dr/i;
        /*printf("%lf %lf\n", pr, res);*/
      }
    /*printf("%f\n", x0);*/
    
    if (zn == 1)
      return res_dr*res_z;
    else
      return 1/res_dr/res_z;
  }

double f_5 (double x)
{
  return exp_my(x);
}

double f_6 (double x)
{
  return 1 / (x * x * 25 + 1);
}


double dd_0 (double x)
{
  (void) x;
  return 0;
}

double dd_1 (double x)
{
  (void) x;
  return 0;
}

double dd_2 (double x)
{
  (void) x;
  return 2;
}

double dd_3 (double x)
{
  return 6 * x;
}

double dd_4 (double x)
{
  return 12 * x * x;
}

double dd_5 (double x)
{
  return exp_my(x);
}

double dd_6 (double x)
{
  return (3750 * x * x - 50) / (25 * x * x + 1) / (25 * x * x + 1) / (25 * x * x + 1);
}


