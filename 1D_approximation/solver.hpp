#include <cstdio>
#include <cmath>


double d_0(double /* x */) {
    return 0;
}
double d_1(double /* x */) {
    return 1.;
}
double d_2(double x) {
    return 2 * x;
}
double d_3(double x) {
    return 3 * x * x;
}
double d_4(double x) {
    return 4 * x * x * x;
}
double d_5(double x) {
    return exp(x);
}
double d_6(double x) {
    return -50 * x / ((25 * x * x + 1) * (25 * x * x + 1));
}





// Pf_1 вызывается отдельно на основе массива альфа, в самом начале вызывается калькулятор, а он вызывает аппроксимацию Чебышева
// калькуляторе мы считаем точки в которых считаем массивы, в этой проге init_memory находит точки в которых надо аппроксимировать
void solve_1(double *f_x, int n, double *alpha, double *z) {
    double g_i_2, g_i_1, g_i;
    int i, j;
    if (n == 1) {
        alpha[0] = f_x[0];
    }
    else if (n == 2) {
        alpha[0] = 0.5 * (f_x[0] + f_x[1]);
        alpha[1] = 0.5 * z[0] * f_x[0] + 0.5 * z[1] * f_x[1]; // 2/n * ...
    }
    else {
        for (i = 0; i < n; i++) alpha[i] = 0;
        
        for (j = 0; j < n; j++) {
            g_i_2 = f_x[j];
            alpha[0] += f_x[j];
            g_i_1 = 0.5 * z[j] * f_x[j];
            alpha[1] += 0.5 * z[j] * f_x[j];
            
            for (i = 2; i < n; i++) {
                g_i = z[j] * g_i_1 - g_i_2;
                alpha[i] += g_i;
                
                g_i_2 = g_i_1;
                g_i_1 = g_i;
            }
        }
        alpha[0] /= n;
        for (i = 1; i < n; i++) alpha[i] = alpha[i] * 2 / n;
    }
}

// аппроксимация Чебышева
double Pf_1(double x, double a, double b, int n, double *alpha) {
    double sum, T_i_2, T_i_1, T_i;
    double z = 2 * (2 * x - (b + a)) / (b - a);
    int i;
    if (n == 1)      return alpha[0];
    else if (n == 2) return alpha[0] + alpha[1] * z * 0.5;
    else {
        sum = 0;
        T_i_2 = 1;
        sum += alpha[0] * T_i_2;
        T_i_1 = z * 0.5;
        sum += alpha[1] * T_i_1;
        
        for (i = 2; i < n; i++) {
            T_i = z * T_i_1 - T_i_2;
            sum += alpha[i] * T_i;
            T_i_2 = T_i_1;
            T_i_1 = T_i;    
        }
            
        return sum;
    }
}







// чисто для второго приближения
void solver_system(int n, double *a, double *c, double *d, double *b, double *x) {
	int i;

	c[0] /= a[0];
	for (i = 1; i < n - 1; i++) {
		a[i] -= d[i - 1] * c[i - 1];
		c[i] /= a[i];
	}
	a[n - 1] -= d[n - 2] * c[n - 2];

	x[0] = b[0] / a[0];
	for (i = 1; i < n; i++) x[i] = (b[i] - d[i - 1] * x[i - 1]) / a[i];

	for (i = n - 2; i >= 0; i--) x[i] -= c[i] * x[i + 1];
}

int bin_search(double x, double *a, int n) {
    int result = 0;
    int r = n;
    int s = (result + r) / 2;
    while (result != r) {
        if (a[s] < x) result = s + 1;
        else r = s;
        s = (result + r) / 2;
    }
    return result;
}



void solve_2(int n, double a, double b, double *x, double *f_x, 
            double *c, double *v, double *ksi, double *a1, double *c1, 
            double *d1, double d0, double dn_1) {
	int i, j = 0;
	double tmp1;

	for (i = 1; i < n; i++)
		ksi[i] = 0.5 * (x[i - 1] + x[i]);

	ksi[0] = a - (ksi[2] - ksi[1]);
	ksi[n] = b + (ksi[n - 1] - ksi[n - 2]);

	for (i = 1; i < n; i++) {
		a1[i] = 1.0/(ksi[i] - x[i - 1]) + 1.0/(ksi[i] - ksi[i - 1]) + 1.0/(x[i] - ksi[i]) + 1.0/(ksi[i + 1] - ksi[i]);
		c1[i] = 1.0/(ksi[i + 1] - x[i]) - 1.0/(ksi[i + 1] - ksi[i]);
		d1[i - 1] = 1.0/(x[i - 1] - ksi[i - 1]) - 1.0/(ksi[i] - ksi[i - 1]);
		c[i] = f_x[i - 1] * (1.0/(x[i - 1] - ksi[i - 1]) + 1.0/(ksi[i] - x[i - 1])) + f_x[i] * (1.0/(x[i] - ksi[i]) + 1.0/(ksi[i + 1] - x[i]));
	}

	a1[0]     = 1. / (ksi[1] -     ksi[0])   - 1. / (x[0] -   ksi[0]);
	c1[0]     = 1. / (ksi[1] -       x[0])   - 1. / (ksi[1] - ksi[0]);
	a1[n]     = 1. / (ksi[n] -   x[n - 1])   - 1. / (ksi[n] - ksi[n - 1]);
	d1[n - 1] = 1. / (ksi[n] - ksi[n - 1])   - 1. / (x[n - 1] - ksi[n - 1]);

    c[0] = d0   - f_x[0]     * (1. / (x[0]     - ksi[0])     - 1. / (ksi[1] - x[0]));
	c[n] = dn_1 - f_x[n - 1] * (1. / (x[n - 1] - ksi[n - 1]) - 1. / (ksi[n] - x[n - 1]));
    // уравнения замыкающие систему, первая производная в граничных условиях

	solver_system(n + 1, a1, c1, d1, c, v);         // j = 0;

	for (i = 0; i < n; i ++) {
		c[j + 0] = v[i];
		tmp1 = ((v[i + 1] - f_x[i]) / (ksi[i + 1] - x[i]) - (f_x[i] - v[i]) / (x[i] - ksi[i])) / (ksi[i + 1] - ksi[i]);
		c[j + 1] = (f_x[i] - v[i]) / (x[i] - ksi[i]) - (x[i] - ksi[i]) * tmp1;
		c[j + 2] = tmp1;
		j += 3;
	}
}

// x = массив чисел x_2
double Pf_2(double t, double *c, double *ksi, int n) {
	int i;
	//for (i = 0; i < n - 1; i++) if (t <= ksi[i + 1]) break;
    i = bin_search(t, ksi, n);
    i -= 1;
	return  c[3 * i] 
            + c[3 * i + 1] * (t - ksi[i]) 
            + c[3 * i + 2] * (t - ksi[i]) * (t - ksi[i]);
}

