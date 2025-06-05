#include "func.hpp"

#define DEF_EPS 1e-16

#define abs_maxiters 100



// основная функция
void *thread_func(void *arg) {
	Args *ptr = (Args *) arg;

	int 	N = ptr->N;
	int 	nx = ptr->nx;
	int 	ny = ptr->ny;
	double 	a = ptr->a;
	double 	b = ptr->b;
	double 	c = ptr->c;
	double 	d = ptr->d;
	int 	p = ptr->p;
	int 	k = ptr->k;
	double 	eps = ptr->eps;
	int		func_id = ptr->func_id;
	int 	maxit = ptr->maxit;

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

	ptr->f = ptr->set_function(func_id);

	double (*f)(double, double) = ptr->f;

	cpu_set_t cpu;		// четыре строки ставят каждому потоку свое ядро процесса
	CPU_ZERO(&cpu);
	pthread_t tid = ptr->tid;
	CPU_SET(p - 1 - (k % p), &cpu);

	pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
	ptr->set_function(func_id);

	double hx = (b - a) / nx, hy = (d - c) / ny;
	ptr->hx = hx, ptr->hy = hy;
	
	double t1 = -1, t2 = -1, r1 = -1, r2 = -1, r3 = -1, r4 = -1, err = 0;

	memset(x, 0, N * sizeof(double));



	t1 = get_full_time();

	if (k == 0) fill_I(nx, ny, I);
	reduce_sum(p);

	err = fill_IA(nx, ny, hx, hy, I, A, p, k);
	reduce_sum(p, &err, 1);
	if(fabs(err) > 0) return nullptr;
	fill_B(N, nx, ny, hx, hy, B, a, c, p, k, f);
	fill_RD(nx, ny, A, I, R, D, p, k);

	reduce_sum(p);

	it = minimal_error_msr_matrix_full(N, A, I, B, R, D, x, myl, r, u, v, eps, maxit, abs_maxiters, p, k); 
	
	t1 = get_full_time() - t1;



	reduce_sum(p);
	t2 = get_full_time();
	all_residuals(r1, r2, r3, r4, x, a, c, hx, hy, nx, ny, p, k, f);
	t2 = get_full_time() - t2;
	reduce_sum(p);

	ptr->t1 = t1;
	ptr->t2 = t2;
	ptr->r1 = r1;
	ptr->r2 = r2;
	ptr->r3 = r3;
	ptr->r4 = r4;
	ptr->all_iters = it;

	return nullptr;
}



int main(int argc, char *argv[]) {
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

	double a = 0, b = 0, c = 0, d = 0, eps = 0;				// параметры
	int nx = 0, ny = 0, func_id = 0, maxit = 0, p = 0;		// параметры

	int i, k;		// итераторы

	if ( !(argc == 11 && 
		sscanf(argv[1], "%lf", &a) == 1 && 
		sscanf(argv[2], "%lf", &b) == 1 && 
		sscanf(argv[3], "%lf", &c) == 1 && 
		sscanf(argv[4], "%lf", &d) == 1 && 
		sscanf(argv[5], "%d", &nx) == 1 && nx > 0 &&
		sscanf(argv[6], "%d", &ny) == 1 && ny > 0 &&
		sscanf(argv[7], "%d", &func_id) == 1 && func_id >= 0 && func_id <= 7 &&
		sscanf(argv[8], "%lf", &eps) == 1 && eps > 0 && 
		sscanf(argv[9], "%d", &maxit) == 1 && maxit > 0 && 
		sscanf(argv[10], "%d", &p) == 1 && p > 0) ) {
		printf("Usage : %s  a  b  c  d  nx  ny  f_id  eps  maxit  p \n", argv[0]);
		return -1;
	}

	Args *duck = nullptr;
	duck = new Args[p];

	double 	*A = nullptr;
	int 	*I = nullptr;
	double 	*R = nullptr;
	double 	*D = nullptr;
	double 	*B = nullptr;
	double 	*x = nullptr;
	double 	*bufer = nullptr;
	double 	*r = nullptr;
	double 	*u = nullptr;
	double 	*v = nullptr;

	int N = (nx + 1) * (ny + 1);
	int len_msr = N + 1 + get_len_msr(nx, ny);

	double t1 = 0, t2 = 0, r1 = -1, r2 = -1, r3 = -1, r4 = -1;

	A = new double[len_msr];
	I = new int[len_msr];
	R = new double[len_msr];		for (i = 0; i < len_msr; i++) R[i] = 0;
	D = new double[len_msr];
	B = new double[N];
	x = new double[N];
	bufer = new double[N];
	r = new double[N];
	u = new double[N];
	v = new double[N];

	init_reduce_sum(p);

	for (k = 0; k < p; k++) {
		duck[k].p = p;
		duck[k].k = k;
		duck[k].A = A;
		duck[k].I = I;
		duck[k].R = R;
		duck[k].D = D;
		duck[k].B = B;
		duck[k].x = x;
		duck[k].bufer = bufer;
		duck[k].r = r;
		duck[k].u = u;
		duck[k].v = v;
		duck[k].N = N;
		duck[k].len_msr = len_msr;
		duck[k].nx = nx;
		duck[k].ny = ny;
		duck[k].a = a;
		duck[k].b = b;
		duck[k].c = c;
		duck[k].d = d;
		duck[k].hx = fabs(b - a) / nx;
		duck[k].hy = fabs(d - c) / ny;
		duck[k].func_id = func_id;
		duck[k].eps = eps;
		duck[k].maxit = maxit;
	}


	for (k = 1; k < p; k++) {
		if (pthread_create(&duck[k].tid, nullptr, thread_func, duck + k)) {
			printf("Ошибочка при создании потока %d\n", k);
			// вообще-то надо чистить память
			return -2;
		}
	}
	duck[0].tid = pthread_self();

	thread_func(duck);

	for (k = 1; k < p; k++) pthread_join(duck[k].tid, nullptr);

	int it = duck->all_iters;
	t1 = duck->t1;
	t2 = duck->t2;
	r1 = duck->r1;
	r2 = duck->r2;
	r3 = duck->r3;
	r4 = duck->r4;

	printf("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n", argv[0], 8, r1, r2, r3, r4, t1, t2, it, eps, func_id, nx, ny, p);

    one_deletion_of_results();
    
	delete [] A;
	delete [] I;
	delete [] R;
	delete [] D;
	delete [] B;
	delete [] x;
    delete [] bufer;
	delete [] r;
	delete [] u;
	delete [] v;
	delete [] duck;
	return 0;
}
