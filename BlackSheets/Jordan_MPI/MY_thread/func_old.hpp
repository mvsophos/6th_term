#include "mpi.h"
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

enum various {
	SUCCESS = 0, ERROR_READ = -1, ERROR_OPEN = -2, SING_MATRIX = -3, ERROR_MEMORY = -4, ERROR = -5
};

// -1 = ошибка чтения фалйа
// -2 = ошибка открытия файла
// -3 = аглоритм неприменим
// -5 = неизвестная ошибка


double f(int i, int j, int n, int s) {
	if (s == 1) return n - (i > j ? i : j);
	if (s == 2) return (i > j ? i : j) + 1;
	if (s == 3) return fabs(i - j);
	if (s == 4) return 1. / (i + j + 1);
	return -1;	// это не выполнится никогда
}

int l2g(int m, int i, int p, int q) { // тут вместо q было k
	return (i / m) * m * p + q * m + (i % m);
}

void set_block(double *A, int i, int j, int n, int m, double *buf, int height, int width) {
	for (int kh = 0; kh < height; kh++) {
		for (int kw = 0; kw < width; kw++) {
			A[(i * m + kh) * n + j * m + kw] = buf[kh * width + kw];
		}
	}
}

void get_block(double *A, int i, int j, int n, int m, double *buf, int height, int width) {
	for (int kh = 0; kh < height; kh++) {
		for (int kw = 0; kw < width; kw++) {
			buf[kh * width + kw] = A[(i * m + kh) * n + j * m + kw];
		}
	}
}

void set_vec(double *B, double *buf, int index, int m, int size) {
	for (int i = 0; i < size; i++) B[index * m + i] = buf[i];
}

void get_vec(double *B, double *buf, int index, int m, int size) {
	for (int i = 0; i < size; i++) buf[i] = B[index * m + i];
}

void change_block_lines(double *A, double *B, int begin, int q, int max, int m, int lines, double *V2, double *V3, double *vec1, double *vec2) {
	for (int i = begin; i < lines; i++) {
		get_block(A,   q, i, lines * m, m, V2, m, m);
		get_block(A, max, i, lines * m, m, V3, m, m);
		set_block(A,   q, i, lines * m, m, V3, m, m);
		set_block(A, max, i, lines * m, m, V2, m, m);
	}

	get_vec(B, vec1,   q, m, m);
	get_vec(B, vec2, max, m, m);
	set_vec(B, vec2,   q, m, m);
	set_vec(B, vec1, max, m, m);
}



// вычитаем блок размером h*w так: A -= C
void extract(double *A, double *C, int height, int width) {
	for (int i = 0; i < height * width; i++) A[i] -= C[i];
}

void set_vector(double *V1, double *V2, int j, int m, int h) {
	for (int i = 0; i < h; i++) V1[i] = V2[j * m + i];
}

void get_vector(double *V1, double *V2, int j, int m, int h) {
	for (int i = 0; i < h; i++) V2[j * m + i] = V1[i];
}

int read_array(FILE *file, double *A, int len) {
	for (int i = 0; i < len; i++) if (fscanf(file, "%lf", A + i) != 1) return -1;
	return 0;
}

int init_matrix(double *A, int n, int m, int k, int lines, int p, int q, int s) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < lines * m; j++) {
			int glob_j = (j / m) * m * p + q * m + (j % m);
			A[i * lines * m + j] = f(i, glob_j, n, s);
		}
	}
	return k + lines;		// k возвращаем просто, чтобы не было предупреждения
}

int init_vector(double *B, double *A, int n, int m, int k, int lines, int p, int q, MPI_Comm comm) {
	double some, s = 0;
	int i, j, j_glob;
	for (i = 0; i < n; i++) {
		s = 0;
		for (j = 0; j < lines * m; j++) { // этот j в локальной нумерации
			j_glob = l2g(m, j, p, q);
			if (j_glob % 2 == 0 && j_glob < n) s += A[i * lines * m + j];
		}
		B[i] = s;
		some = s;
		for (int process = 0; process < p; process++) {
			if (process != q) MPI_Send(&some, 1, MPI_DOUBLE, process, 0, comm);
		}
		MPI_Status status;
		for (int process = 0; process < p; process++) {
			if (process != q) {
				MPI_Recv(&some, 1, MPI_DOUBLE, process, 0, comm, &status);
				//s += some;
				B[i] += some;
			}
		}
		//B[i] = s;
	}
	return k;		// потом исправить, чтобы без предупреждений
}

int print_vector(double *B, int r) {
	for (int i = 0; i < r; i++) printf(" %10.3e", B[i]);
	printf("\n");
	return 0;
}

int read_matrix(double *A, int n, int m, int lines, int p, int q, double *buf, char *name, MPI_Comm comm) {
	int main_q = 0;
	int err = SUCCESS;
	FILE *fp = nullptr;
	
	if (q == main_q) {
		fp = fopen(name, "r");
		if (fp == nullptr) err = ERROR_OPEN;
	}
		
	MPI_Bcast(&err, 1, MPI_INT, main_q, comm);
	
	if (err != 0) return err;
		
	double *block_buf = nullptr;
	block_buf = new double[m * m];
	if (block_buf == nullptr) return ERROR_MEMORY;
		
	int b, max_b = (n + m - 1) / m;
	
	int row_size = (m <= n - (max_b - 1) * m ? m : n - (max_b - 1) * m);		// либо m, а если неровно поделилось, то l = n - k * m
	memset(buf, 0, n * row_size * sizeof(double));
	
	//int blocked_columns_cnt = (n / m) + (n % m != 0);
	//int col_cnt = blocked_columns_cnt / p + (q < blocked_columns_cnt % p);
	
	for (b = 0; b < max_b; b++) {
		int rows = (m <= n - b * m ? m : n - b * m);
		int height = m;
		if (n % m != 0 && b == max_b - 1) height = n % m;
		if (q == main_q) {
			err += read_array(fp, buf, n * rows);
			for (int j = 0; j < max_b; j++) {
				int owner = j % p;
				int width = m;
				int j_loc = j / p;
				if (n % m != 0 && j == max_b - 1) width = n % m;
				get_block(buf, 0, j, n, m, block_buf, height, width);
				
				if (owner == q) {
					set_block(A, b, j_loc, lines * m, m, block_buf, height, width);
				}
				else {
					MPI_Send(block_buf, height * width, MPI_DOUBLE, owner, 0, comm);
				}
			}
		}
		else {
			for (int j = 0; j < max_b; j++) {
				int owner = j % p;
				int width = m;
				int j_loc = j / p;
				if (n % m != 0 && j == max_b - 1)
					width = n % m;
				if (owner == q) {
					MPI_Status st;
					MPI_Recv (block_buf, height * width, MPI_DOUBLE, main_q, 0, comm, &st);
					set_block (A, b, j_loc, lines * m, m, block_buf, height, width);
				}
			}
		}
	}
	if (q == main_q) {
		fclose (fp);
		fp = nullptr;
	}
	MPI_Bcast (&err, 1, MPI_INT, main_q, comm);
	delete [] block_buf;
	if (err != SUCCESS) return err;
	return SUCCESS;
}

// вывод матрицы (main это нулевой поток, ибо он всегда существует)
int print_matrix(double *A, int n, int m, int r, int k, int lines, int p, int q, double *buf, MPI_Comm comm) {
	int main = 0, owner;
	int i, j, j_loc, iter;
	int l;
	for (i = 0; i < r; i++) {
		for (j = 0; j < k; j++) {
			j_loc = j / p;
			owner = j % p;
			l = ((n % m != 0 && j == k - 1) ? n % m : m);
			if (q == main) {
				if (owner == main) {
					for (iter = 0; iter < l && j * m + iter < r; iter++) printf(" %10.3e", A[i * m * lines + j_loc * m + iter]);
				}
				else {
					MPI_Status status;
					MPI_Recv(buf, l, MPI_DOUBLE, owner, 0, comm, &status);
					for (iter = 0; iter < l && j * m + iter < r; iter++) printf(" %10.3e", buf[iter]);
				}
			}
			else {
				if (owner == q) {
					MPI_Send(A + i * m * lines + j_loc * m, l, MPI_DOUBLE, main, 0, comm);
				}
			}
		}
		if (q == main) printf("\n");
	}
	printf("\n");
	return 0;
}



double norm_of_matrix(double *A, int n, int m, int lines, int p, int q, MPI_Comm comm) {
	double norma = 0, all, dod;
	int i, j;
	for (i = 0; i < m; i++) {
		all = 0;
		for (j = 0; j < lines * m; j++) {
			if (l2g(m, j, p, q) < n) all += fabs(A[i * lines * m + j]);
			else break; // удалить 
		}
		for (j = 0; j < p; j++) {		// отправка процессами
			if (j != q) MPI_Send(&all, 1, MPI_DOUBLE, j, 0, comm);
		}
		MPI_Status status;
		dod = 0;
		for (j = 0; j < p; j++) {
			if (j != q) {
				MPI_Recv(&dod, 1, MPI_DOUBLE, j, 0, comm, &status);
				all += dod;
			}
		}
		norma = fmax(norma, all);
	}
	return norma;
}

// ДОБАВИТЬ ФУНКЦИИ НЕВЯЗКИ

// C = результат умножения A * B
int mult_mat_vec(double *A, double *B, double *C, int n, int m, int lines, int p, int q, MPI_Comm comm) {
	double some, s;
	int i, j_loc, j_glob, process;
	for (i = 0; i < n; i++) {
		s = 0;
		for (j_loc = 0; j_loc < lines * m; j_loc++) {
			j_glob = l2g(m, j_loc, p, q);
			if (j_glob < n) s += A[i * m * lines + j_loc] * B[j_glob];
		}
		C[i] = s;
		some = s;
		for (process = 0; process < p; process++) if (process != q) MPI_Send(&some, 1, MPI_DOUBLE, process, 0, comm);
		MPI_Status status;
		for (process = 0; process < p; process++) {
			if (process != q) {
				MPI_Recv(&some, 1, MPI_DOUBLE, process, 0, comm, &status);
				C[i] += some;
			}
		}
	}
	return 0;
}

// B = найденное решение, X = правая часть, AccSol = точное решение, r1 и r2 = указатели на невязки
void residuals(double *r1, double *r2, double *A, double *B, double *X, double *AccSol, int n, int m, int lines, int p, int q, MPI_Comm comm) {
	int i;
	*r1 = 0; *r2 = 0;
	double *C = new double[n];
	double x_norm = 0, rpart_norm = 0;
	mult_mat_vec(A, B, C, n, m, lines, p, q, comm);
	for (i = 0; i < n; i++) {
		//printf(" %lf ", C[i]);
		rpart_norm += fabs(X[i]);
		x_norm += fabs(AccSol[i]);
	}
	//printf("\n");
	for (i = 0; i < n; i++) {
		*r1 += fabs(C[i] - X[i]);
		*r2 += fabs(B[i] - AccSol[i]);
	}

	//printf("%le, %le\n", *r1, *r2);

	*r1 = *r1 / rpart_norm;
	*r2 = *r2 / x_norm;

	delete [] C;
}



// дополнительные функции
int printf_matrix_of_process(double *A, int need, int q, int n1, int n2) {
	if (need == q) {
		for (int i = n1; i < n2; i++) printf("  %lf", A[i]);
		printf("\n");
	}
	return 0;
}




// ФУНКЦИИ НЕОБХОДИМЫЕ ДЛЯ РЕШЕНИЯ СИСТЕМЫ
void swap(double *lhs, double *rhs) {
	double buf = *lhs;
	*lhs = *rhs;
	*rhs = buf;
}

double norma_of_matrix(double *A, int h, int w) {		// n = ширина матрицы, m = высота матрицы
	double norma = 0;
	double s;
	for (int i = 0; i < h; i++) {
		s = 0;
		for (int j = 0; j < w; j++) {
			s += fabs(A[i * w + j]);
		}
		if (s > norma) norma = s;
	}
	return norma;
}

int jord_inverse(double *A, double *E, int n, double ERROR) {
	int i, j; // итераторы
    double max = 0.; // максимальный элемент по столбцу
    int m = 0; // строчка и столбец с максимальным элементом
    double c = 0.; // коэффициент
    
    for (j = 0; j < n; j++) { // идем по столбам
        max = fabs(A[j * n + j]);
        m = j;

        for (i = j; i < n; i++) {
            if (fabs(A[i * n + j]) > max) {
                m = i;
                max   = fabs(A[i * n + j]);
            }
        } // нашли максимальный элемент и строку в которой он находится

        if (max <= ERROR) { // ошибка нужна чтобы нормально найти максимальный элемент
            return -1;
        }

        if (m != j) {
            for (i = 0; i < j; i++) { // меняем строки местами
                swap(&(E[j * n + i]), &(E[m * n + i])); // элемент на дигонали мы не меняем так как там все равно будет 1 в (j,j), а все остальные 0
            }
            for (i = j; i < n; i++) { // менять их в матрице А не надо так как там нули
                swap(&(A[j * n + i]), &(A[m * n + i]));
                swap(&(E[j * n + i]), &(E[m * n + i]));
            }
        }
        //c = 1. / A(j, j);
        c = A[j * n + j]; // этот элемент на диагонали не равен 0
        //printf("%lf    ", c);
        //if (fabs(c - 1) > ERROR) {
            for (m = 0; m < j + 1; m++) E[j * n + m] /= c; // m больше не нужен так как строки уже переставили
            A[j * n + j] = 1.;
            for (m = j + 1; m < n; m++) {
                A[j * n + m] /= c;
                E[j * n + m] /= c;
            }
        //}

        // цикл ниже вычитает обнуляет почти весь столб
        for (i = 0; i < n; i++) { // все элементы в столбе j кроме одного равны 0
            if (i != j){
                c = A[i * n + j];
                //if (fabs(c) > ERROR) {
                    for (m = 0; m < j + 1; m++) {
                        E[i * n + m] -= c * E[j * n + m];
                    }
                    // A(i, j) = 0.;
                    for (m = j + 1; m < n; m++) {
                        A[i * n + m] -= c * A[j * n + m];
                        E[i * n + m] -= c * E[j * n + m];
                    }
                //}
            }
        }
    } // теперь в Rev находится обратная матрица решенная Жорданом
	return 0;
}





// умножение матрицы на вектор: A * vec = result_vec. (h = высота матрицы (и результата), w = ширина матрицы)
int multiply_mat_vec(double *A, double *vec, double *result_vec, int h, int w) {
	double s;
	for (int i = 0; i < h; i++) {
		s = 0;
		for (int j = 0; j < w; j++) {
			s += A[i * w + j] * vec[j];
		}
		result_vec[i] = s;
	}
	return 0;
}

// c = a * b
int multiply_blocks(double *a, int av, int ah, double *b, int bv, int bh, double *c) {
    int r = 0, t = 0, q, temp_r, temp_t;
    // для вышеописанной логики
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    if (ah != bv)
        return -1;
    // Зануляем c
    for (r = 0; r < av; r++)
        for (t = 0; t < bh; t++)
            c[r * bh + t] = 0.;

    for (r = 0; r + 2 < av; r += 3)
        for (t = 0; t + 2 < bh; t += 3) {
            c00 = 0., c01 = 0., c02 = 0.;
            c10 = 0., c11 = 0., c12 = 0.;
            c20 = 0., c21 = 0., c22 = 0.;
            for (q = 0; q < ah; q++) {
                c00 += a[(r + 0) * ah + q] * b[q * bh + (t + 0)];
                c01 += a[(r + 0) * ah + q] * b[q * bh + (t + 1)];
                c02 += a[(r + 0) * ah + q] * b[q * bh + (t + 2)];
                c10 += a[(r + 1) * ah + q] * b[q * bh + (t + 0)];
                c11 += a[(r + 1) * ah + q] * b[q * bh + (t + 1)];
                c12 += a[(r + 1) * ah + q] * b[q * bh + (t + 2)];
                c20 += a[(r + 2) * ah + q] * b[q * bh + (t + 0)];
                c21 += a[(r + 2) * ah + q] * b[q * bh + (t + 1)];
                c22 += a[(r + 2) * ah + q] * b[q * bh + (t + 2)];
            }
            c[(r + 0) * bh + (t + 0)] += c00;
            c[(r + 0) * bh + (t + 1)] += c01;
            c[(r + 0) * bh + (t + 2)] += c02;
            c[(r + 1) * bh + (t + 0)] += c10;
            c[(r + 1) * bh + (t + 1)] += c11;
            c[(r + 1) * bh + (t + 2)] += c12;
            c[(r + 2) * bh + (t + 0)] += c20;
            c[(r + 2) * bh + (t + 1)] += c21;
            c[(r + 2) * bh + (t + 2)] += c22;
        }
    temp_t = t;
    temp_r = r;
    // если не получилось разделить четно, то
    // повторяем процесс для последнего столбца и строчки
    // так, как делали раньше
    if ((av - temp_r) == 2) {
        for (t = 0; t + 1 < bh; t += 2) {
            c00 = 0., c01 = 0.;
            c10 = 0., c11 = 0.;
            for (q = 0; q < ah; q++) {
                c00 += a[(temp_r + 0) * ah + q] * b[q * bh + (t + 0)];
                c01 += a[(temp_r + 0) * ah + q] * b[q * bh + (t + 1)];
                c10 += a[(temp_r + 1) * ah + q] * b[q * bh + (t + 0)];
                c11 += a[(temp_r + 1) * ah + q] * b[q * bh + (t + 1)];
            }
            c[(temp_r + 0) * bh + (t + 0)] = c00;
            c[(temp_r + 0) * bh + (t + 1)] = c01;
            c[(temp_r + 1) * bh + (t + 0)] = c10;
            c[(temp_r + 1) * bh + (t + 1)] = c11;
        }
        if (t < bh) {
            for (r = temp_r; r < av; r++)
                for (t = bh & ~1; t < bh; t++) {
                    c00 = 0.;
                    for (q = 0; q < ah; q++) {
                        c00 += a[(r + 0) * ah + q] * b[q * bh + (t + 0)];
                    }
                    c[(r + 0) * bh + (t + 0)] = c00;
                }
        }
    }

    if ((bh - temp_t) == 2) {
        for (r = 0; r + 1 < av; r += 2) {
            c00 = 0., c01 = 0.;
            c10 = 0., c11 = 0.;
            for (q = 0; q < ah; q++) {
                c00 += a[(r + 0) * ah + q] * b[q * bh + (temp_t + 0)];
                c01 += a[(r + 0) * ah + q] * b[q * bh + (temp_t + 1)];
                c10 += a[(r + 1) * ah + q] * b[q * bh + (temp_t + 0)];
                c11 += a[(r + 1) * ah + q] * b[q * bh + (temp_t + 1)];
            }
            c[(r + 0) * bh + (temp_t + 0)] = c00;
            c[(r + 0) * bh + (temp_t + 1)] = c01;
            c[(r + 1) * bh + (temp_t + 0)] = c10;
            c[(r + 1) * bh + (temp_t + 1)] = c11;
        }
        for (r = av & ~1; r < av; r++) {
            for (t = temp_t; t < bh; t++) {
                c00 = 0.;
                for (q = 0; q < ah; q++) {
                    c00 += a[(r + 0) * ah + q] * b[q * bh + (t + 0)];
                }
                c[(r + 0) * bh + (t + 0)] = c00;
            }
        }
    }

    if (av - temp_r == 1) {
        for (t = 0; t < bh; t++) {
            c00 = 0.;
            for (q = 0; q < ah; q++)
                c00 += a[temp_r * ah + q] * b[q * bh + t];
            c[temp_r * bh + t] = c00;
        }
    }
    if (bh - temp_t == 1) {
        for (r = 0; r < av; r++) {
            c00 = 0.;
            for (q = 0; q < ah; q++)
                c00 += a[r * ah + q] * b[q * bh + temp_t];
			c[r * bh + temp_t] = c00;
        }
    }
    return 0;
}



int get_inversed_block(double *buf, int n, int m, int size, int i, int j, double *V1, double *V2, double norma) {
	get_block(buf, i, j, n, m, V2, m, m);
	for (int b = 0; b < size * size; b++) {
		V1[b] = ((b / size) == (b % size)  ?  1 : 0);
	}
	return jord_inverse(V2, V1, m, norma);
}

int get_main_block(std::pair <double, int> &glavny, /* double *glavny, */ double *buf, /* int n, */ int m, int size, /* int lines, */ int p, int q, int step, double *V1, double *V2, double norma) {
	/* glavny[0] = 1e304;
	glavny[1] = -1; */

	glavny.first = 1e304;
	glavny.second = -1;


	double minimum;
	int result;

	for (int i = step + q; i < size; i += p) {
		result = get_inversed_block(buf, m, m, m, i, 0, V1, V2, norma);

		/* if (result == -1) continue;
		minimum = norma_of_matrix(V1, m, m);
		if (minimum < glavny[0]) {
			glavny[0] = minimum;
			glavny[1] = i;
		} */

		if (result == 0) {
			minimum = norma_of_matrix(V1, m, m);
			if (minimum < glavny.first) {
				glavny.first = minimum;
				glavny.second = i;
			}

			/* if (minimum < glavny[0]) {
				glavny[0] = minimum;
				glavny[1] = i;
			} */
		}
	}
	return 0;
}

void subtract_with_multiply (double *a, int from_ind, int ind, int start_stl, double *mult, int col_cnt, int m, double *buf, double *result) {
	for (int w = start_stl; w < col_cnt; w++) {
		get_block(a, ind, w, col_cnt * m, m, buf, m, m);
		multiply_blocks(mult, m, m, buf, m, m, result);
		//multiply_blocks_m (mult, buf, res, m);
		get_block(a, from_ind, w, col_cnt * m, m, buf, m, m);
		for (int i = 0; i < m * m; i++) buf[i] -= result[i];
		//subtract_blocks (buf, res, m, m);
		set_block(a, from_ind, w, col_cnt * m, m, buf, m, m);
	}
}

void subtract_with_multiply_vector (double *b, int from_ind, int ind, double *mult, int m, double *buf, double *result) {
	get_vector(b, buf, ind, m, m);
	multiply_blocks(mult, m, m, buf, m, m, result);
	//multiply_matrix_and_vector (mult, buf, res, m, m);
	get_vector(b, buf, from_ind, m, m);
	for (int i = 0; i < m * 1; i++) buf[i] -= result[i];
	//subtract_blocks (buf, res, m, 1);
	set_vector(b, buf, from_ind, m, m); 
}



int solver(double *A, double *B, double *buf, double norma, int n, int m, int k, int lines, int p, int q, MPI_Comm comm) {
	int size = n / m;		// s = целые блочные строки (k когда нацело, и k + 1 когда не поделилось)
	int l = n % m;			// l = остаток (если разделилось нацело, то его нет и он не используется)
	int owner;
	int result;

	printf("---- 0\n");

	double *V1 = new double[m * m];
	double *V2 = new double[m * m];
	double *V3 = new double[m * m];
	double *V1_l = new double[l * l];
	double *V2_l = new double[l * l];

	double *buf_vec_1 = new double[m];
	double *buf_vec_2 = new double[m];



	for (int i = 0; i < k; i++) {
		owner = i % p;

		// отправляем всем столб, в котором надо искать главный элемент
		if (owner == q) {				// копируем столбец, чтобы отправить его остальным
			for (int d = 0; d < k * m; d++) {
				for (int e = 0; e < m; e++) buf[d * m + e] = A[d * lines * m + (i / p) * m + e];
			}
		}
		MPI_Bcast(buf, m * m * k, MPI_DOUBLE, owner, comm);

		// ищем главный элемент в отправленном буфере
		if (l != 0 && i == size && owner == q) {
			get_block(A, size, i / p, lines * m, m, V1, l, l);
			if (norma_of_matrix(V1, l, l) < norma) result = SING_MATRIX;
			else result = size;
		}
		else if (i < size || l == 0) {
			//double *sendbuf = new double[2], *recvbuf = new double[2];
			std::pair <double, int> send, get;
			get_main_block(send, buf, /* n, */ m, size, /* lines, */ p, q, i, V1, V2, norma);
			printf("---- 1\n");
			MPI_Allreduce(&send, &get, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
			printf("---- 2\n");
			result = get.second;
			//result = (int)recvbuf[1];
		}
		MPI_Bcast(&result, 1, MPI_INT, owner, comm);

		// если алгоритм неприменим, то чистим память
		if (result != 0) {
			// надо очистить память
			delete [] V1;
			delete [] V2;
			delete [] V3;
			delete [] V1_l;
			delete [] V2_l;
			delete [] buf_vec_1;
			delete [] buf_vec_2;
			return result;
		}


		if (owner == q) {
			if (l != 0 && i == size) {
				get_inversed_block(A, lines * m, m, l, result, size / p, V1_l, V2_l, norma);
				memset(V1, 0, m * m * sizeof(double));
				set_block(V1, 0, 0, m, m, V1_l, l, l);
			}
			else get_inversed_block(A, lines * m, m, m, result, i / p, V1, V2, norma);
		}

		// отправим 
		MPI_Bcast(V1, m * m, MPI_DOUBLE, owner, comm);

		// умножаем строку и правую часть на обратный к главному элементу
		// умножаем строку
		for (int j = (i / p) + (owner == q); j < lines; j++) {
			get_block(A, result, j, lines * m, m, V2, m, m);
			multiply_blocks(V1, m, m, V2, m, m, V3);
			set_block(A, result, j, lines * m, m, V3, m, m);
		}
		// умножаем правую часть
		for (int y = 0; y < m; y++) buf_vec_1[y] = *(B + result * m + y);
		for (int y = 0; y < m; y++) {
			double sum = 0;
			for (int u = 0; u < m; u++) {
				sum += V1[y * m + u] * buf_vec_1[u];
			}
			B[result * m + y] = sum;
		}



		// далее вычитаем строки
		for (int j = 0; j < k; j++) {
			if (j == result) continue;
			get_block(buf, j, 0, m, m, V1, m, m);
			//if (norma_of_matrix(V1, m, m) < norma) continue;
			subtract_with_multiply(A, j, result, i / p + (owner == q), V1, lines, m, V2, V3);
			subtract_with_multiply_vector(B, j, result, V1, m, buf_vec_1, buf_vec_2);
		}
		if (result != i) change_block_lines(A, B, i / p + (owner == q), i, result, m, lines, V2, V3, buf_vec_1, buf_vec_2);
	}

	//MPI_Barrier(comm);
	//printf("%lf %lf %lf, %le, %d %d %d %d\n", A[0], B[0], buf[0], norma, n, m, p, q);
	delete [] V1; delete [] V2; delete [] V3; delete [] V1_l; delete [] V2_l; delete [] buf_vec_1; delete [] buf_vec_2;
	return 0;
	return result;
}


