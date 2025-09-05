#include <pthread.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
#include <fenv.h>				// для обработки исключений с помощью feenableexcept
#include "func.hpp"

double function_0(double, double) 		{ return 1.0; }
double function_1(double x, double) 	{ return x; }
double function_2(double, double y) 	{ return y; }
double function_3(double x, double y) 	{ return x + y; }
double function_4(double x, double y) 	{ return sqrt(x*x + y*y); }
double function_5(double x, double y) 	{ return x*x + y*y; }
double function_6(double x, double y) 	{ return exp(x*x - y*y); }
double function_7(double x, double y) 	{ return 1.0 / (1 + 25*(x*x + y*y)); }





// pf это как функция у меня Pf
/* 
class Gui_data{
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
	int nx, ny, mx, my;
	int m;
	int maxit;
	int p;
	double (*f)(double, double) = nullptr;

	void realloc_data();
	double find_min_max(double &abs_min, double &abs_max);
	double pf(double, double);
	void pfind_min_max(double &, double &);
	void residual_min_max(double &, double &);

	~Gui_data() {
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
};



class Args {
public:
	void copy_data(const Gui_data &data){
		A = data.A;
		I = data.I;
		R = data.R;
		D = data.D;
		B = data.B;
		x = data.x;
		bufer = data.bufer;
		u = data.u;
		v = data.v;
		r = data.r;
		a = data.a;
		b = data.b;
		c = data.c;
		d = data.d;
		eps = data.eps;
		nx = data.nx;
		ny = data.ny;
		f = data.f;
		p = data.p;
		maxit = data.maxit;
	}

	int nx = 0, ny = 0, N = 0, len_msr = 0, func_id = 0;
	double a = 0, b = 0, c = 0, d = 0, hx = 0, hy = 0;
	int p = 0, k = 0;

	double 	*A = nullptr;
	int 	*I = nullptr;
	double	*D = nullptr;		// это для Холецкого
	double	*R = nullptr;		// это для Холецкого
	double 	*B = nullptr;
	double 	*x = nullptr;
	double	*bufer = nullptr;
	double 	*r = nullptr;
	double 	*u = nullptr;
	double 	*v = nullptr;

	double (*f)(double, double) = nullptr;

	int its = 0, maxit = 0;
	double eps = 0;
	double t1 = 0, t2 = 0,r1 = -1, r2 = -1, r3 = -1, r4 = -1;

	bool close_window = false;
	bool ready = true;

	pthread_t tid = 0;
	pthread_cond_t *cond = nullptr;
	pthread_mutex_t *mutex = nullptr;
	double cpu_time = 0;
	double cpu_time_of_all_threads = 0;
	double astr_time = 0.0;					//////////////////////
	double res = 0;

	double (*set_function(int ooo))(double, double) {
		switch (ooo) {
		case 0: return function_0;
		case 1: return function_1;
		case 2: return function_2;
		case 3: return function_3;
		case 4: return function_4;
		case 5: return function_5;
		case 6: return function_6;
		case 7: return function_7;
		}
		printf("Ошибка в set_function\n");
		return 0;
	}
};
*/



void print_array(double *array, int n) {
	for (int i = 0; i < n; i++) printf(" %lf ", array[i]);
	printf("\n");
}




//////////////// important_func.hpp ////////

// МЕТОДЫ В ЭТОМ ФАЙЛЕ БЫЛИ НАПИСАНЫ НА ЛЕКЦИИ, ВСЕ КРОМЕ ПРЕДОБУСЛАВЛИВАТЕЛЯ, ХОЛЕЦКОГО И Pf (ошибок нет)

void reduce_sum(int p, double *a, int n) {
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;
	static int t_out = 0;
	static double *r = nullptr;
	int i;
	if (p <= 1) return;
	pthread_mutex_lock(&mutex);
	if (r == nullptr) r = a;
	else for (i = 0; i < n; i++) r[i] += a[i];

	t_in++;

	if (t_in >= p) {
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else while (t_in < p) pthread_cond_wait(&c_in, &mutex);

	if(r != a) for (i = 0; i < n; i++) a[i] = r[i];

	t_out++;
	
	if (t_out >= p) {
		t_in = 0;
		r = nullptr;
		pthread_cond_broadcast(&c_out);
	}
	else while (t_out < p) pthread_cond_wait(&c_out, &mutex);
	pthread_mutex_unlock(&mutex);
}

void reduce_sum_int(int p, int *a, int n) {
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;
	static int t_out = 0;
	static int *r = nullptr;
	int i;
	if (p <= 1) return;
	pthread_mutex_lock(&mutex);
	if (r == nullptr) r = a;
	else for (i = 0; i < n; i++) r[i] += a[i];

	t_in++;

	if (t_in >= p) {
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else while (t_in < p) pthread_cond_wait(&c_in, &mutex);

	if(r != a) for (i = 0; i < n; i++) a[i] = r[i];

	t_out++;
	
	if (t_out >= p) {
		t_in = 0;
		r = nullptr;
		pthread_cond_broadcast(&c_out);
	}
	else while (t_out < p) pthread_cond_wait(&c_out, &mutex);
	pthread_mutex_unlock(&mutex);
}

void reduce_max(int p, double *a, int n) {
	static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;
	static int t_out = 0;
	static double* r = nullptr;
	int i;

	if (p <= 1)	return;
	pthread_mutex_lock(&my_mutex);
	if (r == nullptr) r = a;
	else for (i = 0; i < n; i++) r[i] = (r[i] > a[i] ? r[i] : a[i]);

	t_in++;

	if (t_in >= p) {
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else while (t_in < p) pthread_cond_wait(&c_in, &my_mutex);

	if (r != a) for (i = 0; i < n; i++) a[i] = r[i];

	t_out++;
	if (t_out >= p) {
		t_in = 0;
		r = nullptr;
		pthread_cond_broadcast(&c_out);
	}
	else while (t_out < p) pthread_cond_wait(&c_out, &my_mutex);
	pthread_mutex_unlock(&my_mutex);
}





// эти штуки написаны верно, они нужны

static double *results = nullptr;
int init_reduce_sum(int p) {
	results = new double[p];
	if (results == nullptr) return -1;
	return 0;
}
double reduce_sum_det(int p, int k, double s) {
	double sum = 0.;
	results[k] = s;
	reduce_sum(p);
	for (int l = 0; l < p; l++) sum += results[l];
	reduce_sum(p);
	return sum;
}
void one_deletion_of_results() {
	delete [] results;
}

double get_full_time() {
	struct timeval buf;
	gettimeofday(&buf, 0);
	return buf.tv_sec + buf.tv_usec * 1e-6;
}

double get_cpu_time() {
	struct rusage buf;
	getrusage(RUSAGE_THREAD, &buf);
	return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec * 1e-6;
}


// функции для работы функций для работы в thread_func

void thread_rows(int n, int p, int q, int &l1, int &l2) {
	l1 = n * q; l1 = l1 / p;
	l2 = n * (q + 1); l2 = l2 / p;
}

void ij2l(int nx, int /* ny */, int i, int j, int &l) {
	nx += 1;
	l = i + j * nx;
}

void l2ij(int nx, int /* ny */, int &i, int &j, int l) {
	nx += 1;
	j = l / nx;
	i = l - j * nx;
}

// общая длина (количество диагональных, внедиагональных и плюс один дополнительный)
int get_len_msr(int nx, int ny) {
	nx -= 1; ny -= 1;
	return 6 * nx * ny + 4 * 2 * (nx + ny) + 3 * 2 + 2 * 2;
}



int IA_ij(int nx, int ny, double hx, double hy, int i, int j, int is, int js, int s, int *I, double *A) {
	int l, ls;
	ij2l(nx, ny, i,  j,  l);
	ij2l(nx, ny, is, js, ls);

	if (I != nullptr) I[s] = ls;

	if (A != nullptr) {
		// то есть один и тот же элемент
		if (l == ls) {
			A[s] = 0; 
			if (i < nx && j > 0)	A[s] += hx * hy / 12;
			if (i > 0 && j > 0)		A[s] += 2 * hx * hy / 12;
			if (i > 0 && j < ny)	A[s] += hx * hy / 12;
			if (i < nx && j < ny)	A[s] += 2 * hx * hy / 12;
		}
		else {
			A[s] = 0;
			if (is == i + 1 && js == j) {
				if (j < ny) A[s] += hx * hy / 24;
				if (j > 0)  A[s] += hx * hy / 24;
			}
			else if (is == i && js == j - 1) {
				if (i < nx) A[s] += hx * hy / 24;
				if (i > 0)  A[s] += hx * hy / 24;
			}
			else if (is == i - 1 && js == j - 1) A[s] = 2 * hx * hy / 24;

			else if (is == i - 1 && js == j) {
				if (j < ny) A[s] += hx * hy / 24;
				if (j > 0)  A[s] += hx * hy / 24;
			}
			else if (is == i && js == j + 1) {
				if (i > 0)  A[s] += hx * hy / 24;
				if (i < nx) A[s] += hx * hy / 24;
			}
			else if (is == i + 1 && js == j + 1) A[s] = 2 * hx * hy / 24;
			
			else return -1;	// чистая формальность, чтобы компилятор не истерил
		}
	}
	return 0;
}

#define F(IS, JS, S) IA_ij(nx, ny, hx, hy, i, j, (IS), (JS), (S), I, A)
int get_off_diag(int nx, int ny, double hx, double hy, int i, int j, int *I, double *A) {
	int s = 0;
	if (i < nx) 			{ F(i + 1, j, s); s++; }
	if(j > 0) 				{ F(i, j-1,s); s++; }
	if(i > 0 && j > 0)		{ F(i-1, j-1, s); s++; }
	if(i > 0)				{ F(i-1, j, s); s++; }
	if(j < ny)				{ F(i, j+1, s); s++; }
	if(i < nx && j < ny)	{ F(i+1, j+1, s); s++; }
	return s;
}

#define Fb(I, J) (f(a + (I) * hx, c + (J) * hy))
double F_IJ(int nx, int ny, double hx, double hy, double a, double c, int i, int j, double (*f)(double, double)) {
	double w = hx * hy / 192;
	if (i > 0 && i < nx && j > 0 && j < ny) {
		return w * (
			36 * Fb(i, j)
			+ 20 * ( Fb(i + 0.5, j) + Fb(i, j - 0.5) + Fb(i - 0.5, j - 0.5) + Fb(i - 0.5, j) + Fb(i, j + 0.5) + Fb(i + 0.5, j + 0.5) )
			+ 4 *  ( Fb(i + 0.5, j - 0.5) + Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5) + Fb(i - 0.5, j + 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5) )
			+ 2 *  ( Fb(i + 1, j) + Fb(i, j - 1) + Fb(i - 1, j - 1) + Fb(i - 1, j) + Fb(i, j + 1) + Fb(i + 1, j + 1) )  );
	}
	if (i > 0 && i < nx && j == 0) {
		return w * (
			18 * Fb(i, j)
			+ 10 * ( Fb(i + 0.5, j) + Fb(i - 0.5, j) )
			+ 20 * ( Fb(i, j + 0.5) + Fb(i + 0.5, j + 0.5) )
			+ 4 *  ( Fb(i - 0.5, j + 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5) )
			+ 1 *  ( Fb(i - 1, j) + Fb(i + 1, j) ) + 2 * ( Fb(i, j + 1) + Fb(i + 1, j + 1) )  );
	}
	if (i > 0 && i < nx && j == ny) {
		return w * (
			18 * Fb(i, j)
			+ 10 * ( Fb(i - 0.5, j) + Fb(i + 0.5, j) )
			+ 20 * ( Fb(i, j - 0.5) + Fb(i - 0.5, j - 0.5) )
			+ 4 *  ( Fb(i + 0.5, j - 0.5) + Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5) )
			+ 1 *  ( Fb(i - 1, j) + Fb(i + 1, j) )
			+ 2 *  ( Fb(i, j - 1) + Fb(i - 1, j - 1) )  );
	}
	if (i == 0 && j > 0 && j < ny) 
	{
		return w * (
			18 * Fb(i, j)
			+ 10 * ( Fb(i, j - 0.5) + Fb(i, j + 0.5) )
			+ 20 * ( Fb(i + 0.5, j) + Fb(i + 0.5, j + 0.5) )
			+ 4 *  ( Fb(i + 0.5, j - 0.5) + Fb(i + 0.5, j + 1) + Fb(i + 1, j + 0.5) )
			+ 1 *  ( Fb(i, j - 1) + Fb(i, j + 1) )
			+ 2 *  ( Fb(i + 1, j) + Fb(i + 1, j + 1) )  );
	}
	if (i == nx && j > 0 && j < ny) 
	{
		return w * (
			18 * Fb(i, j)
			+ 10 * ( Fb(i, j - 0.5) + Fb(i, j + 0.5) )
			+ 20 * ( Fb(i - 0.5, j) + Fb(i - 0.5, j - 0.5) )
			+ 4 *  ( Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5) + Fb(i - 0.5, j + 0.5) )
			+ 1 *  ( Fb(i, j - 1) + Fb(i, j + 1) )
			+ 2 *  ( Fb(i - 1, j) + Fb(i - 1, j - 1) )  );
	}
	if (i == 0 && j == 0) {
		return w * (
			12 * Fb(i, j)
			+ 10 * ( Fb(i + 0.5, j) + Fb(i, j + 0.5) )
			+ 20 * ( Fb(i + 0.5, j + 0.5) )
			+ 4 *  ( Fb(i + 1, j + 0.5) + Fb(i + 0.5, j + 1) )
			+ 1 *  ( Fb(i + 1, j) + Fb(i, j + 1) )
			+ 2 *  ( Fb(i + 1, j + 1) )  );
	}
	if (i == nx && j == ny) {
		return w * (
			12 * Fb(i, j)
			+ 10 * ( Fb(i - 0.5, j) + Fb(i, j - 0.5) )
			+ 20 * ( Fb(i - 0.5, j - 0.5) )
			+ 4 *  ( Fb(i - 0.5, j - 1) + Fb(i - 1, j - 0.5) )
			+ 1 *  ( Fb(i, j - 1) + Fb(i - 1, j) )
			+ 2 *  ( Fb(i - 1, j - 1) )  );
	}
	if (i == 0 && j == ny) {
		return w * (
			6 * Fb(i, j)
			+ 10 * ( Fb(i + 0.5, j) + Fb(i, j - 0.5) )
			+ 4 *  ( Fb(i + 0.5, j - 0.5) )
			+ 1 *  ( Fb(i + 1, j) + Fb(i, j - 1) )  );
	}
	if (i == nx && j == 0) {
		return w * (
			6 * Fb(i, j)
			+ 10 * ( Fb(i - 0.5, j) + Fb(i, j + 0.5) )
			+ 4 *  ( Fb(i - 0.5, j + 0.5) )
			+ 1 *  ( Fb(i - 1, j) + Fb(i, j + 1) )  );
	}
	return 1e308;		// никогда не выполнится, просто чтобы компилятор не истерил
}



// вроде эту функцию вообще не требуется использовать, хотя мы ее написали
int get_len_msr_off_diag(int nx, int ny, double *A, int *I) {
	double hx = 0, hy = 0;
	int i, j, res = 0;
	for(i = 0; i<=nx; i++) {
		for (j = 0; j < ny; j++) res += get_off_diag(nx, ny, hx, hy, i, j, I, A);
	}
	return res;
}

// функции для работы в thread_func
void fill_I(int nx, int ny, int *I) {
	int i, j, N = (nx + 1) * (ny + 1);
	int r = N + 1;
	double hx = 0, hy = 0;

	for (int l = 0; l < N; l++) {
		l2ij(nx, ny, i, j, l);
		int s = get_off_diag(nx, ny, hx, hy, i, j, 0, 0);
		I[l] = r;
		r += s;
	}
	I[N] = r;
}

int fill_IA(int nx, int ny, double hx, double hy, int *I, double *A, int p, int k) {
	int i, j, l, l1, l2, N = (nx + 1) * (ny + 1), r, s, t;
	double error = 0;
	int len = 0;
	thread_rows(N, p, k, l1, l2);

	for (l = l1; l < l2; l++) {
		r = I[l];
		s = I[l+1] - I[l];		// количество не равных 0 элементов в строке номер l (в нумерации от 0)
		l2ij(nx, ny, i, j, l);

		if (IA_ij(nx, ny, hx, hy, i, j, i, j, 0, nullptr, A + l) != 0) {
			error = 1;
			break;
		}
		t = get_off_diag(nx, ny, hx, hy, i, j, I + r, A + r);
		if (t != s) {
			error = -1;
			break;
		}
		len += s;
	}
	reduce_sum(p, &error, 1);
	if(error < 0) return -1;
	reduce_sum_int(p, &len, 1);
	if (I[N] != N + 1 + len) return -2;
	return 0;
}

void fill_B(int N, int nx, int ny, double hx, double hy, double *b, double a, double c, int p, int k, double (*f)(double, double)) {
	int i1 = 0, i2 = 0, l, i = 0, j = 0;
	thread_rows(N, p, k, i1, i2);
	for(l = i1; l < i2; l++) {
		l2ij(nx, ny, i, j, l);
		b[l] = F_IJ(nx, ny, hx, hy, a, c, i, j, f); 
	}
	reduce_sum(p);
}



//////////////// for_precond.hpp ////////

double sign(double x) {
	if (x > 0) return  1;
	if (x < 0) return -1;
	else return 0;
}





double Pf(double *res, double x, double y, double a, double c, double hx, double hy, int nx, int ny) {
	int i, j, l0, l1, l2;
	double x0, y0, z0;
	double x1, y1, z1;
	double x2, y2, z2;
	double dx0, dy0, dx1, dy1, dz1, dx2, dy2, dz2;

	double det_1, det_2, condition;

	i = (x - a) / hx;	j = (y - c) / hy;

	ij2l(nx, ny, i,     j,     l0);
	ij2l(nx, ny, i + 1, j + 1, l2);

	x0 = a + i * hx;	y0 = c + j * hy;	z0 = res[l0];
	condition = hy * (x - x0) / hx;		condition += y0 - y;

	if (condition >= 0) {
		ij2l(nx, ny, i + 1, j, l1);
		x1 = a + (i + 1) * hx;
		y1 = c + j * hy;
	}
	else {
		ij2l(nx, ny, i, j + 1, l1);
		x1 = a + i * hx;
		y1 = c + (j + 1) * hy;
	}
	z1 = res[l1];
	ij2l(nx, ny, i + 1, j + 1, l2);

	x2 = a + (i + 1) * hx;		y2 = c + (j + 1) * hy;		z2 = res[l2];

	dx0 = x  - x0;	dy0 = y  - y0;
	dx1 = x1 - x0;	dy1 = y1 - y0;	dz1 = z1 - z0;
	dx2 = x2 - x0;	dy2 = y2 - y0;	dz2 = z2 - z0;

	det_1 = dy1*dx2*z0 + dy0*dz1*dx2 + dx0*dy1*dz2 - dx1*dy2*z0 - dx0*dz1*dy2 - dy0*dx1*dz2;
	det_2 = dy1*dx2 - dx1*dy2;

	return det_1 / det_2;
}

// i на самом деле это ряд (а не столб), j наоборот является итератором по столбам

double r1_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k) {
	double result = 0, buf = 0;
	int i, j, begin, end;
	thread_rows(ny /* - 1 */, p, k, begin, end);

	for (i = begin; i < end; ++i) {			// по рядам
		for (j = 0; j < nx; j++) {			// по столбам
			buf = fabs( f(a + (j + 2.0 / 3) * hx, c + (i + 1.0 / 3) * hy) - Pf(x, a + (j + 2.0 / 3) * hx, c + (i + 1.0 / 3) * hy, a, c, hx, hy, nx, ny) );
			result = fmax(buf, result);
			buf = fabs( f(a + (j + 1.0 / 3) * hx, c + (i + 2.0 / 3) * hy) - Pf(x, a + (j + 1.0 / 3) * hx, c + (i + 2.0 / 3) * hy, a, c, hx, hy, nx, ny) );
			result = fmax(buf, result);
		}
	}
	reduce_max(p, &result, 1);
	return result;
}

double r2_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k) {
	double result = 0;
	int i, j, begin, end;
	double Sq = 0.5 * hx * hy;
	thread_rows(ny /* - 1 */, p, k, begin, end);

	for (i = begin; i < end; i++) {			// по рядам
		for (j = 0; j < nx; j++) {			// по столбам
			result += fabs( f(a + (j + 2.0 / 3) * hx, c + (i + 1.0 / 3) * hy) - Pf(x, a + (j + 2.0 / 3) * hx, c + (i + 1.0 / 3) * hy, a, c, hx, hy, nx, ny) );
			result += fabs( f(a + (j + 1.0 / 3) * hx, c + (i + 2.0 / 3) * hy) - Pf(x, a + (j + 1.0 / 3) * hx, c + (i + 2.0 / 3) * hy, a, c, hx, hy, nx, ny) );
		}
	}
	result *= Sq;
	reduce_sum(p, &result, 1);
	return result;
}

double r3_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k) {
	double result = 0;
	int i, j, begin, end, l;
	thread_rows(ny + 1, p, k, begin, end);

	for (i = begin; i < end; i++) {			// по рядам
		for (j = 0; j <= nx; j++) {			// по столбам
			ij2l(nx, ny, j, i, l);
			result = fmax( result, fabs( f(a + j * hx, c + i * hy) - x[l] ) );
		}
	}
	reduce_max(p, &result, 1);
	return result;
}

double r4_nev(double (*f)(double, double), double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k) {
	double result = 0;
	int i, j, begin, end, l;
	double Sq = hx * hy;
	thread_rows(ny+1, p, k, begin, end);  // Исправлено: ny+1 вместо ny

	for (i = begin; i < end; i++) {
		for (j = 0; j <= nx; j++) {
			ij2l(nx, ny, j, i, l);
			result += fabs(f(a + j * hx, c + i * hy) - x[l]);  // Суммирование разностей
		}
	}
	result *= Sq;
	reduce_sum(p, &result, 1);  // Правильная редукция для суммы
	return result;
}

void all_residuals(double &r1, double &r2, double &r3, double &r4, double *x, double a, double c, double hx, double hy, int nx, int ny, int p, int k, double (*f)(double, double)) {
	r1 = r1_nev(f, x, a, c, hx, hy, nx, ny, p, k);
	r2 = r2_nev(f, x, a, c, hx, hy, nx, ny, p, k);
	r3 = r3_nev(f, x, a, c, hx, hy, nx, ny, p, k);
	r4 = r4_nev(f, x, a, c, hx, hy, nx, ny, p, k);
}





// нахождение матрицы R и D (разложение Холецкого)
void fill_RD (int nx, int ny, double *A, int *I, double *R, double *D, int p, int k) {
	int N = (int)(nx + 1) * (ny + 1), l1, l2, l, len_s, len_j, len_i;
	int i, j, s, q, m, iter_i, iter_j;			// итераторы
	int count = 0;
	double temp, sum;

	thread_rows(N, p, k, l1, l2);

	for (i = l1; i < l2; i++) {
		sum = 0;
		for (s = l1; s < i; s++) {
			len_s = I[s + 1] - I[s];
			for (l = 0; l < len_s; l++) if (I[I[s] + l] == i) break;

			if (l < len_s) sum += D[s] * R[I[s] + l] * R[I[s] + l];
			else count += 1;
		}

		temp = A[i] - sum;
		R[i] = sqrt(fabs(temp));		// диагональные элементы
		D[i] = sign(temp);

		len_i = I[i + 1] - I[i];
		for (m = 0; m < len_i; m++) {
			j = I[I[i] + m];

			if (j > i && l1 <= j && j < l2) {
				sum = 0;
				for (s = l1; s < i; s++) {
					len_s = I[s + 1] - I[s];
					for (iter_i = 0; iter_i < len_s; iter_i++) if (I[I[s] + iter_i] == i) break;
					for (iter_j = 0; iter_j < len_s; iter_j++) if (I[I[s] + iter_j] == j) break;

					if (iter_i < len_s && iter_j < len_s) sum += R[I[s] + iter_i] * D[s] * R[I[s] + iter_j];
				}
				R[I[i] + m] = (A[I[i] + m] - sum) / (R[i] * D[i]);

				len_j = I[j + 1] - I[j];
				for (q = 0; q < len_j; q++) if (I[I[j] + q] == i) break;

				if (q < len_j) R[I[j] + q] = R[I[i] + m];
			}
		}
	}
	reduce_sum(p);
}



/////////////// method.hpp //////////////

// ВСЕ КРОМЕ ПРЕДПОСЛЕДНЕГО МЕТОДА И ПРЕДОБУСЛАВЛИВАТЕЛЯ БЫЛО НАПИСАНО НА ЛЕКЦИИ

double scalar_product(int n, double *x, double *y, int p, int k) {
	int i, i1, i2;  double s = 0;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; i++) s += x[i] * y[i];
	// было написано reduce_sum(p, &s, 1);
	return reduce_sum_det(p, k, s);
}

// применить предобуславливатель Холецкого (сам)
void apply_precondition_msr_matrix(int n, double * /* A */, int *I, double *R, double *D, double *vec, double *v, double *r, int p, int k) {
	int i1, i2, len_i, len;
	double sum;
	thread_rows (n, p, k, i1, i2);
	// для предобуславливателя Якоби : for (i = i1; i < i2; i++) v[i] = r[i] / A[i];
	// для Якоби все работает и так

	for (int i = i1; i < i2; i++) {
		sum = 0;
		len_i = I[i + 1] - I[i];

		for (int s = 0; s < len_i; s++)	{
			int j = I[I[i] + s];
			if ((j < i) && (i1 <= j && j < i2)) sum += R[I[i] + s] * vec[j];
		}
		vec[i] = (r[i] - sum) / R[i];
	}
	for (int i = i1; i < i2; i++) vec[i] *= D[i];

	reduce_sum(p);

	len = i2 - i1;
	for (int j = 0, i = (i2 - 1); j < len; i--, j++) {
		sum = 0;
		len_i = I[i + 1] - I[i];

		for (int s = 0; s < len_i; s++) {
			int j = I[I[i] + s];
			if ((j > i) && (i1 <= j && j < i2)) sum += R[I[i] + s] * v[j];
		}
		v[i] = (vec[i] - sum) / R[i];
	}

	reduce_sum(p);
}



void mult_sub_vector(int n, double *x, double *y, double t, int p, int k) {
	int l1, l2;
	thread_rows(n, p, k, l1, l2);
	for (int l = l1; l < l2; l++) x[l] -= t * y[l];
	reduce_sum(p);
}

void matrix_mult_vector_msr(double *A, int *I, int n, double *x, double *y, int p, int k) {
	int i, j, l, J, i1, i2;
	double s;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; i++) {
		s = A[i] * x[i];		// A_ii * x_i = диагональный элемент
		l = I[i+1] - I[i];		// количество ненулевых внедиагональных элементов
		J = I[i];				// начало строки i
		for (j = 0; j < l; j++)	s += A[J + j] * x[I[J + j]]; // I[J + j] - номер столбца для A[j + j]
		y[i] = s;
	}
	reduce_sum(p);
}

// было в начале лекции №6, но переделано из минимальных невязок в минимальные ошибки
// отличие от метода минимальных невязок состоит только в том какие вектора там скалярно умножать, методы похожи
int minimal_error_msr_matrix(int N, double *A, int *I, double *B, double *R, double *D, double *x /* ответ в вектор x */, double *vec /* дополнительный вектор */, double *r, double *u, double *v, double eps, int maxit, int p, int k) {
	double prec = 0, b_norm2 = 0, tau = 0;
	int it = 0;
	double c1 = 0, c2 = 0;
	b_norm2 = scalar_product(N, B, B, p, k);
	prec = b_norm2 * eps * eps;
	
	matrix_mult_vector_msr(A, I, N, x, r, p, k);		// r = Ax
	mult_sub_vector(N, r, B, 1., p, k);					// r -= 1 * b

	// до этого момента все сходится

	for (it = 0; it < maxit; it++) {
		apply_precondition_msr_matrix(N, A, I, R, D, vec, v, r, p, k);	// Mv = r
		//print_array(v, N);
		matrix_mult_vector_msr(A, I, N, v, u, p, k);			// u = Av, u = AM^(-1) r
		//print_array(u, N);
		c1 = scalar_product(N, v, r, p, k);			// c1 = (v, r)
		c2 = scalar_product(N, u, v, p, k);			// c2 = (u, v)
		if (c1 < prec || c2 < prec) break;				//// если критерий выполнен ////
		tau = c1 / c2;
		mult_sub_vector(N, x, v, tau, p, k);		// x -= \tau * v;
		mult_sub_vector(N, r, u, tau, p, k);		// r -= \tau * u
	}
	if (it >= maxit) return -1;
	return it;
}


//int minimal_error_msr_matrix_full(int n, double *A, int *I, double *R, double *D, double *B, double *x /* ответ в вектор x */,	double *vec, double *r,	double *u, double *v, double eps, int maxit, int maxstep, int p, int k) {
int minimal_error_msr_matrix_full(int n, double *A, int *I, double *B, double *R, double *D, double *x /* ответ в вектор x */,	double *vec, double *r,	double *u, double *v, double eps, int maxit, int maxstep, int p, int k) {
	int step, ret, its = 0;
	for (step = 0; step < maxstep; step++) {
		ret = minimal_error_msr_matrix(n, A, I, B, R, D, x, vec, r, u, v, eps, maxit, p, k); 
		if (ret >= 0) {
			its += ret;
			break;
		}
		its += maxit;
	}
	if (step >= maxstep) return -1;
	return its;
}

