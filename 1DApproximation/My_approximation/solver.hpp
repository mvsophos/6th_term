#include <cstdio>
#include <cmath>




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
  return exp(x);
}

double dd_6 (double x)
{
  return (3750 * x * x - 50) / (25 * x * x + 1) / (25 * x * x + 1) / (25 * x * x + 1);
}



void solve_1 (double x[], double f_x[], int n, double alpha[], double z[])
{
  //printf("SOLVE_01\n");
  (void) x;
  double g_i_2, g_i_1, g_i;
  int i, j;
  if (n == 1)
    {
      alpha[0] = f_x[0];
    }
  else
  if (n == 2)
    {
      alpha[0] = 0.5 * (f_x[0] + f_x[1]);
      alpha[1] = 0.5 * z[0] * f_x[0] + 0.5 * z[1] * f_x[1]; // 2/n * ...
    }
  else
    {
      for (i = 0; i < n; i++)
        alpha[i] = 0;
      
      for (j = 0; j < n; j++)
        {
          g_i_2 = f_x[j];
          alpha[0] += f_x[j];
          g_i_1 = 0.5 * z[j] * f_x[j];
          alpha[1] += 0.5 * z[j] * f_x[j];
          
          for (i = 2; i < n; i++)
            {
              g_i = z[j] * g_i_1 - g_i_2;
              alpha[i] += g_i;
              
              g_i_2 = g_i_1;
              g_i_1 = g_i;
            }
          
        }
        
      alpha[0] /= n;
      for (i = 1; i < n; i++)
        alpha[i] = alpha[i] * 2 / n;
      
    }
  /*if (n <= 50)
    {
      for (i = 0; i < n; i++)
          printf(" %10.3e", alpha[i]);
      printf("\n");
    }*/
}

double Pf_1 (double x, double a, double b, int n, double x_mas[], double alpha[])
{
  (void) x_mas;
  double sum, t_i_2, t_i_1, t_i;
  double z = 2 * (2 * x - (b + a)) / (b - a);
  int i;
  if (n == 1)
    {
      return alpha[0]; // ... * 1
      
    }
  else
  if (n == 2)
    {
      return alpha[0] + alpha[1] * z * 0.5;
      
    }
  else
    {
      sum = 0;
      t_i_2 = 1;
      sum += alpha[0] * t_i_2;
      t_i_1 = z * 0.5;
      sum += alpha[1] * t_i_1;
      
      for (i = 2; i < n; i++)
        {
          t_i = z * t_i_1 - t_i_2;
          sum += alpha[i] * t_i;
          t_i_2 = t_i_1;
          t_i_1 = t_i;
          
        }
          
      return sum;
    }
}

void solve_2 (double x[], double f_x[], double mas_4n[], int n, double dd0, double ddn)
{
  //printf("SOLVE_02\n");
  (void) x;
  double f_xk_xk1, f_xk1_xk2, f_xk1_xk, dk1, dk;
  int k;
  if (n == 2)
    {
      k = 0;
      f_xk_xk1 = (f_x[k+1] - f_x[k]) / (x[k+1] - x[k]);
      //printf("f_xk_xk1 %10.3e %10.3e\n", f_x[k+1], f_x[k]);
      //printf("f_xk_xk1 %10.3e %10.3e %10.3e\n", x[k+1], x[k], f_xk_xk1);
      double a = 3 * f_xk_xk1 - 0.5 * dd0 * (x[k+1] - x[k]);
      double b = 3 * f_xk_xk1 + 0.5 * ddn * (x[n-1] - x[n-2]); // ddn1
      dk = (2 * a - b) / 3;
      dk1 = (2 * b - a) / 3;
      
      mas_4n[0] = f_x[0];
      mas_4n[1] = dk;
      mas_4n[2] = (3 * f_xk_xk1 - 2 * dk - dk1) / (x[k+1] - x[k]);
      mas_4n[3] = (dk + dk1 - 2 * f_xk_xk1) / (x[k+1] - x[k]) / (x[k+1] - x[k]);
      
    }
  else
  if (n > 2)
    {
      k = 0;
      f_xk_xk1 = (f_x[k+1] - f_x[k]) / (x[k+1] - x[k]);
      f_xk1_xk2 = (f_x[k+2] - f_x[k+1]) / (x[k+2] - x[k+1]);
      //printf("f_xk_xk1 %10.3e %10.3e\n", f_x[k+1], f_x[k]);
      //printf("f_xk_xk1 %10.3e %10.3e %10.3e\n", x[k+1], x[k], f_xk_xk1);
      dk1 = (fabs(f_xk_xk1) < fabs(f_xk1_xk2) ? f_xk_xk1 : f_xk1_xk2);
      dk = 0.5 * (- dk1 + 3 * f_xk_xk1 - 0.5 * dd0 * (x[k+1] - x[k]));
      
      mas_4n[0] = f_x[0];
      mas_4n[1] = dk;
      mas_4n[2] = (3 * f_xk_xk1 - 2 * dk - dk1) / (x[k+1] - x[k]);
      mas_4n[3] = (dk + dk1 - 2 * f_xk_xk1) / (x[k+1] - x[k]) / (x[k+1] - x[k]);
      
      //printf(" %d %10.3e %10.3e %10.3e %10.3e\n", k, mas_4n[4*k], mas_4n[4*k+1], mas_4n[4*k+2], mas_4n[4*k+3]);
      f_xk_xk1 = f_xk1_xk2;
      f_xk1_xk = f_xk_xk1;
      dk = dk1;
      for (k = 1; k < n - 2; k++) // #### n-1
        {
          f_xk1_xk2 = (f_x[k+2] - f_x[k+1]) / (x[k+2] - x[k+1]);
          if (f_xk_xk1 * f_xk1_xk2 > 0)
            {
              dk1 = (fabs(f_xk_xk1) < fabs(f_xk1_xk2) ? f_xk_xk1 : f_xk1_xk2);
            }
          else
            {
              dk1 = 0;
            }
          
          mas_4n[4*k + 0] = f_x[k];
          mas_4n[4*k + 1] = dk;
          mas_4n[4*k + 2] = (3 * f_xk_xk1 - 2 * dk - dk1) / (x[k+1] - x[k]);
          mas_4n[4*k + 3] = (dk + dk1 - 2 * f_xk_xk1) / (x[k+1] - x[k]) / (x[k+1] - x[k]);
          //printf(" %d %10.3e %10.3e %10.3e %10.3e\n", k, mas_4n[4*k], mas_4n[4*k+1], mas_4n[4*k+2], mas_4n[4*k+3]);
          //printf("f_xk_xk1 %10.3e f_xk1_xk2 %10.3e\n", f_xk_xk1, f_xk1_xk2);
          f_xk_xk1 = f_xk1_xk2;
          f_xk1_xk = f_xk_xk1;
          dk = dk1;
        }
      k = n - 2;
      dk1 = 0.5 * (- dk + 3 * f_xk1_xk + 0.5 * ddn * (x[n-1] - x[n-2])); // ddn1
      //printf("f_xk1_xk %10.3e\n", f_xk1_xk);
      mas_4n[4*k + 0] = f_x[n-2];
      mas_4n[4*k + 1] = dk;
      mas_4n[4*k + 2] = (3 * f_xk1_xk - 2 * dk - dk1) / (x[k+1] - x[k]);
      mas_4n[4*k + 3] = (dk + dk1 - 2 * f_xk1_xk) / (x[k+1] - x[k]) / (x[k+1] - x[k]);
      //printf(" %d %10.3e %10.3e %10.3e %10.3e\n", k, mas_4n[4*k], mas_4n[4*k+1], mas_4n[4*k+2], mas_4n[4*k+3]);
      //printf("DD %10.3e %10.3e\n", dd0, ddn);
    }
  
      
  //for (i = 0; i < 4*n; i++)
  //  printf(" %10.3e", alpha[i]);
  //printf("\n");
}

int bin_search_to_add (double x, double * a, int n)
{
  int l = 0, r = n, s = (l + r) / 2;
  while (l != r)
    {
      if (a[s] < x)
        {
          l = s + 1;
        }
      else
        {
          r = s;
        }
      s = (l + r) / 2;
    }
  return l;
}

double Pf_2 (double x, double a, double b, int n, double x_mas[], double mas_4n[])
{
  //printf("PF_02\n");
  (void) x_mas;
  (void) a;
  (void) b;
  int i;
  if (n == 2)
    {
      //return 0;//alpha[0] + alpha[1] * z * 0.5;
      i = 0;
      return mas_4n[4*i] + mas_4n[4*i+1]*(x-x_mas[i]) 
                         + mas_4n[4*i+2]*(x-x_mas[i])*(x-x_mas[i])
                         + mas_4n[4*i+3]*(x-x_mas[i])*(x-x_mas[i])*(x-x_mas[i]);
    }
  else
  if (n > 2)
    {
      i = bin_search_to_add(x, x_mas, n) - 1;
      if (i < 0)
        i = 0;
      if (i == n - 1)
        i = n - 2;
      return mas_4n[4*i] + mas_4n[4*i+1]*(x-x_mas[i]) 
                         + mas_4n[4*i+2]*(x-x_mas[i])*(x-x_mas[i])
                         + mas_4n[4*i+3]*(x-x_mas[i])*(x-x_mas[i])*(x-x_mas[i]);
      for (i = 0; i < n-2; i++)
        if (x_mas[i] <= x && x <= x_mas[i+1])
          {
          //printf("PF_02 %10.3e %10.3e %10.3e\n", x_mas[i], x, x_mas[i+1]);
          return mas_4n[4*i] + mas_4n[4*i+1]*(x-x_mas[i]) 
                             + mas_4n[4*i+2]*(x-x_mas[i])*(x-x_mas[i])
                             + mas_4n[4*i+3]*(x-x_mas[i])*(x-x_mas[i])*(x-x_mas[i]);
          }
      for (i = 0; i < n-2; i++)
        if (x_mas[i] <= x && x <= x_mas[i+1])
          {
          //printf("PF_02 %10.3e %10.3e %10.3e\n", x_mas[i], x, x_mas[i+1]);
          return mas_4n[4*i] + mas_4n[4*i+1]*(x-x_mas[i]) 
                             + mas_4n[4*i+2]*(x-x_mas[i])*(x-x_mas[i])
                             + mas_4n[4*i+3]*(x-x_mas[i])*(x-x_mas[i])*(x-x_mas[i]);
          }
          
      return mas_4n[4*i] + mas_4n[4*i+1]*(x-x_mas[i]) + mas_4n[4*i+2]*(x-x_mas[i])*(x-x_mas[i]) + mas_4n[4*i+3]*(x-x_mas[i])*(x-x_mas[i])*(x-x_mas[i]);
    }
  return 0;
}

