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







