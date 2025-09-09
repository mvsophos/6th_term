#ifndef thread_hpp
#define thread_hpp

#include <pthread.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sys/sysinfo.h>



#define DEF_EPS 1e-13
#define EPS 1e-12

#define abs_maxiters 200

void *thread_func(void *arg);

class Dataset {
public:
	double *A = nullptr;
	int    *I = nullptr;
	double *D = nullptr;
	double *R = nullptr;
	double *B = nullptr;
	double *x = nullptr;
	double *bufer = nullptr;
	double *u = nullptr;
	double *v = nullptr;
	double *r = nullptr;

	double a, b, c, d;
	double eps;
    int func_id;
	int nx, ny, mx, my;
	int m;
	int maxit;
	int p;
	double (*f)(double, double) = nullptr;
	
	// это для увеличения погрешности в центре
	int parameter;
	double norma;

	void repeat_new_data();
	double find_min_max(double &abs_min, double &abs_max);
	double pf(double, double);
	void pfind_min_max(double &, double &);
	void residual_min_max(double &, double &);

	~Dataset() {
		if (A) delete [] A;
		if (I) delete [] I;
		if (R) delete [] R;
		if (D) delete [] D;
		if (B) delete [] B;
		if (x) delete [] x;
		if (bufer) delete [] bufer;
		if (u) delete [] u;
		if (v) delete [] v;
		if (r) delete [] r;
	}
	//~Dataset();
};

class Args {
public:
	void copy_data(const Dataset &data);
	
	int nx = 0, ny = 0, N = 0, len_msr = 0, func_id = 0;
	double a = 0, b = 0, c = 0, d = 0, hx = 0, hy = 0;
	int p = 0, k = 0;

	int parameter = 0;
	double norma = 0;

	double  *A = nullptr;
	int     *I = nullptr;
	double  *D = nullptr;
	double  *R = nullptr;
	double  *B = nullptr;
	double  *x = nullptr;
	double  *bufer = nullptr;
	double  *r = nullptr;
	double  *u = nullptr;
	double  *v = nullptr;

	double (*f)(double, double) = nullptr;

	int its = 0, maxit = 0;
	double eps = 0;
	double t1 = 0, t2 = 0, r1 = -1, r2 = -1, r3 = -1, r4 = -1;

	bool close_window = false;
	bool ready = true;

	pthread_t tid = 0;
	pthread_cond_t *cond = nullptr;
	pthread_mutex_t *mutex = nullptr;
	double cpu_time = 0;
	double cpu_time_of_all_threads = 0;
	double astr_time = 0.0;
	double res = 0;

	double (*set_function(int ooo))(double, double);
};

#endif