#ifndef FUNC_HPP
#define FUNC_HPP

#include <pthread.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
#include <fenv.h>



double function_0(double, double);
double function_1(double x, double);
double function_2(double, double y);
double function_3(double x, double y);
double function_4(double x, double y);
double function_5(double x, double y);
double function_6(double x, double y);
double function_7(double x, double y);

void print_array(double *array, int n);
void reduce_sum(int p, double *a = nullptr, int n = 0);
void reduce_sum_int(int p, int *a = nullptr, int n = 0);
void reduce_max(int p, double *a = nullptr, int n = 0);
int init_reduce_sum(int p);
double reduce_sum_det(int p, int k, double s);
void one_deletion_of_results();
double get_full_time();
double get_cpu_time();
void thread_rows(int n, int p, int q, int &l1, int &l2);
void ij2l(int nx, int ny, int i, int j, int &l);
void l2ij(int nx, int ny, int &i, int &j, int l);
int get_len_msr(int nx, int ny);
int IA_ij(int nx, int ny, double hx, double hy, int i, int j, int is, int js, int s, int *I = nullptr, double *A = nullptr);
int get_off_diag(int nx, int ny, double hx, double hy, int i, int j, int *I, double *A);
double F_IJ(int nx, int ny, double hx, double hy, double a, double c, int i, int j, double (*f)(double, double), int parameter, double norma);
int get_len_msr_off_diag(int nx, int ny, double *A, int *I);
void fill_I(int nx, int ny, int *I);
int fill_IA(int nx, int ny, double hx, double hy, int *I, double *A, int p, int k);
void fill_B(int N, int nx, int ny, double hx, double hy, double *b, double a, double c, int p, int k, double (*f)(double, double), int parameter, double norma);
double sign(double x);
double Pf(double *res, double x, double y, double a, double c, double hx, double hy, int nx, int ny);
double r1_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k);
double r2_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k);
double r3_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k);
double r4_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k);
void all_residuals(double &r1, double &r2, double &r3, double &r4, double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double));
void fill_RD(int nx, int ny, double *A, int *I, double *R, double *D, int p, int k);
double scalar_product(int n, double *x, double *y, int p, int k);
void apply_precondition_msr_matrix(int n, double *A, int *I, double *R, double *D, double *vec, double *v, double *r, int p, int k);
void mult_sub_vector(int n, double *x, double *y, double t, int p, int k);
void matrix_mult_vector_msr(double *A, int *I, int n, double *x, double *y, int p, int k);
int minimal_error_msr_matrix(int N, double *A, int *I, double *B, double *R, double *D, double *x, double *vec, double *r, double *u, double *v, double eps, int maxit, int p, int k);
int minimal_error_msr_matrix_full(int n, double *A, int *I, double *B, double *R, double *D, double *x, double *vec, double *r, double *u, double *v, double eps, int maxit, int maxstep, int p, int k);

#endif // FUNC_HPP