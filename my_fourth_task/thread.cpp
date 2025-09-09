#include "thread.hpp"
#include "func.hpp"


// написал эти функции тут, чтобы определить их, чтобы можно было модульно подключать все .hpp файлы

void Dataset::repeat_new_data() {
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

	I = nullptr;
	A = nullptr;
	
	int N = (nx + 1) * (ny + 1);
	int len = N + get_len_msr(nx, ny) + 1;
	A = new double[len];
	I = new int[len];
	D = new double[len];
	R = new double[len];		for (int i = 0; i < len; i++) R[i] = 0;

	init_reduce_sum(p);
	B = new double[N];
	x = new double[N];
	bufer = new double[N];
	r = new double[N];
	u = new double[N];
	v = new double[N];

	fill_I(nx, ny, I);
	memset(x, 0, N * sizeof(double));
}

double Dataset::find_min_max(double &abs_min, double &abs_max) {
	double hx = (b - a) / mx, hy = (d - c) / my;
	double value1, value2;
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			value1 = f(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);
			value2 = f(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);
			abs_min = fmin(abs_min, value1);
			abs_max = fmax(abs_max, value1);
			abs_min = fmin(abs_min, value2);
			abs_max = fmax(abs_max, value2);
		}
	}
	return fmax(abs_min, abs_max);
}

double Dataset::pf(/* double *res,  */double x_coord, double y_coord /* , double a, double c, double hx, double hy, int nx, int ny */) {
	int i, j, l0, l1, l2;
	double hx = (b - a) / nx, hy = (d - c) / ny;
	double x0, y0, z0;
	double x1, y1, z1;
	double x2, y2, z2;
	double dx0, dy0, dx1, dy1, dz1, dx2, dy2, dz2;

	//printf("ДАННЫЕ: \n%d %d\n\n", nx, ny);

	double det_1, det_2, condition;

	i = (x_coord - a) / hx;	j = (y_coord - c) / hy;

	ij2l(nx, ny, i,     j,     l0);
	ij2l(nx, ny, i + 1, j + 1, l2);

	x0 = a + i * hx;	y0 = c + j * hy;	z0 = x[l0];
	condition = hy * (x_coord - x0) / hx;	condition += y0 - y_coord;

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
	z1 = x[l1];
	ij2l(nx, ny, i + 1, j + 1, l2);

	x2 = a + (i + 1) * hx;		y2 = c + (j + 1) * hy;		z2 = x[l2];

	dx0 = x_coord  - x0;	dy0 = y_coord  - y0;
	dx1 = x1 - x0;			dy1 = y1 - y0;			dz1 = z1 - z0;
	dx2 = x2 - x0;			dy2 = y2 - y0;			dz2 = z2 - z0;

	det_1 = dy1*dx2*z0 + dy0*dz1*dx2 + dx0*dy1*dz2 - dx1*dy2*z0 - dx0*dz1*dy2 - dy0*dx1*dz2;
	det_2 = dy1*dx2 - dx1*dy2;

	/* printf("ЗначениR = %le %le %le\n", z0, z1, z2);
	printf("\n%d %d %d\n\n", l0, l1, l2); */

	return det_1 / det_2;
}

void Dataset::pfind_min_max(double &abs_min, double &abs_max) {
	double hx = (b - a) / mx, hy = (d - c) / my;
	double value1, value2;
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			value1 = pf(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);
			value2 = pf(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);
			abs_min = fmin(abs_min, value1);
			abs_max = fmax(abs_max, value1);
			abs_min = fmin(abs_min, value2);
			abs_max = fmax(abs_max, value2);
		}
	}
}

void Dataset::residual_min_max(double &abs_min, double &abs_max) {
	double hx = (b - a) / mx, hy = (d - c) / my;
	double f_val, pf_val, value1, value2;
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			f_val = f(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);
			pf_val = pf(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);
			value1 = fabs(f_val - pf_val);

			f_val = f(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);
			pf_val = pf(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);
			value2 = fabs(f_val - pf_val);

			abs_min = fmin(abs_min, value1);
			abs_max = fmax(abs_max, value1);
			abs_min = fmin(abs_min, value2);
			abs_max = fmax(abs_max, value2);
		}
	}
}





// Реализация Args::copy_data
void Args::copy_data(const Dataset &data) {
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
	func_id = data.func_id;
	p = data.p;
	maxit = data.maxit;
	parameter = data.parameter;
	norma = data.norma;
}

double (*Args::set_function(int ooo))(double, double) {
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
	return nullptr;  // или return 0;
}

////////////////// НАДО ПЕРЕПИСАТЬ ЭТУ ФУНКЦИЮ ПРАВИЛЬНО, ПО СУТИ ЭТО РЕШЕНИЕ

// основная функция
void *thread_func(void *arg) {
	Args *ptr = (Args *) arg;
	pthread_cond_t *cond = ptr->cond;
	pthread_mutex_t *mutex = ptr->mutex;

	int q = ptr->k, p = ptr->p;



	cpu_set_t cpu;
	CPU_ZERO(&cpu);
	int n_cpus = get_nprocs();
	int cpu_id = n_cpus - 1 - (q % n_cpus);
	CPU_SET(cpu_id, &cpu);
	pthread_t tid = pthread_self();
	pthread_setaffinity_np(tid, sizeof(cpu), &cpu);



	while (!(ptr->close_window)) {
		if (!ptr->ready) {
			int 	N = ptr->N;		// бессмысленно --- это число тут нулевое
			int 	nx = ptr->nx;
			int 	ny = ptr->ny;
			double 	a = ptr->a;
			double 	b = ptr->b;
			double 	c = ptr->c;
			double 	d = ptr->d;
			int 	k = ptr->k;
			double 	eps = ptr->eps;
			int		func_id = ptr->func_id;
			int 	maxit = ptr->maxit;

			//printf("this %d\n", func_id);
			if (k == 0) printf("\n");

			int it;
			
			double 	*A = ptr->A;
			double 	*B = ptr->B;
			double 	*R = ptr->R;
			double 	*D = ptr->D;
			double 	*x = ptr->x;
			double 	*myl = ptr->bufer;
			double 	*r = ptr->r;
			double 	*u = ptr->u;
			double 	*v = ptr->v;
			int 	*I = ptr->I;

			int parameter = ptr->parameter;	//if (k == 0) printf("param  %d\n", parameter);
			double norma = ptr->norma;		//if (k == 0) printf("normal %lf\n", norma);

			ptr->f = ptr->set_function(func_id);

			

			double (*f)(double, double) = ptr->f;

			//printf("k = %d\n", k);
			N = (nx + 1) * (ny + 1);

			if (k == 0) init_reduce_sum(p);
			reduce_sum(p);

			


			
			//ptr->set_function(func_id);
			//printf("function = %lf\n", (*f)(1,1));

			double hx = (b - a) / nx, hy = (d - c) / ny;
			ptr->hx = hx, ptr->hy = hy;

			// это нужно, чтобы найти норму функции (а она нужна для (я думаю можны бы сделать))
			
			double t1 = -1, t2 = -1, r1 = -1, r2 = -1, r3 = -1, r4 = -1, err = 0;

			memset(x, 0, N * sizeof(double));



			t1 = get_full_time();

			if (k == 0) fill_I(nx, ny, I);
			/* reduce_sum(p); */

			err = fill_IA(nx, ny, hx, hy, I, A, p, k);
			
			// /* вывод просто так */ printf("\n МАТРИЦА А:   "); for (int kol = 0; kol < N; kol++) printf(" [ %lf %d ] ", A[kol], I[kol]); printf("\n\n");
			
			reduce_sum(p, &err, 1);
			if (fabs(err) > 0) return nullptr;
			fill_B(N, nx, ny, hx, hy, B, a, c, p, k, f, parameter, norma);		// f используется здесь для задания B, поэтому надо чтобы она была нормальной
			fill_RD(nx, ny, A, I, R, D, p, k);

			//printf("Values: %lf %lf %lf\n", B[N - 2], B[4], B[1]);
			// /* вывод просто так */ printf("\n"); for (int kol = 0; kol < N; kol++) printf(" [ %lf %lf ] ", B[kol], x[kol]); printf("\n\n");

			/* reduce_sum(p); */

			//printf("it = %d\n", maxit);

			it = minimal_error_msr_matrix_full(N, A, I, B, R, D, x, myl, r, u, v, eps, maxit, abs_maxiters, p, k);
			//printf("it = %d\n", it);
			
			t1 = get_full_time() - t1;

			//printf("Values: %lf %lf %lf\n", x[0], myl[4], x[11]);
			// проблема с массивом x[], он полностью нулевой, а не должен быть (там -1, 0 и -0.5)



			/* reduce_sum(p); */
			t2 = get_full_time();
			all_residuals(r1, r2, r3, r4, x, a, c, hx, hy, nx, ny, p, k, f);
			t2 = get_full_time() - t2;
			/* reduce_sum(p); */

			ptr->t1 = t1;
			ptr->t2 = t2;
			ptr->r1 = r1;
			ptr->r2 = r2;
			ptr->r3 = r3;
			ptr->r4 = r4;
			ptr->its = it;

			ptr->ready = true;
			reduce_sum(p);

			pthread_mutex_lock(mutex);
			while (!(ptr->ready)) {
				pthread_cond_wait(cond, mutex);
			}
			pthread_mutex_unlock(mutex);
		}
	}

	ptr->ready = true;
	return nullptr;
}


