#ifndef HEADER2  
#define HEADER2
//#include "algorithm.hpp"

class Args
{
public:
    int p = 0;
    int k = 0;
    MPI_Comm comm = MPI_COMM_WORLD;
    double cpu_time = 0;
    double cpu_time_of_all_threads = 0;
    double astr_time = 0.0;
    double res = 0;

    double* A = nullptr;
    double* B = nullptr; 
    double* buf = nullptr;
    bool* ZerosMatrix = nullptr;
    double* U = nullptr;
    double* ProductResult = nullptr;
    double* ZeroMatrix = nullptr;
    double* results = nullptr;

    double norm = 0;
    bool last_line_isnt_fool = false;

    std::string name = "";

    int n = 0; //Matrix dim
    int m = 0; //Block dim
    int M = 0; //Amount of block strings
    int r = 0; //Amount of printing elements
    int s = 0;
    int l = 0;
    int K = 0;
    int rows = 0;
    int max_rows = 0;
    int active_procceses = 0; 

    int cur_str = 0;
    int cur_global_str = 0;
    int nomer_v_okne = 0;
    int shag = 0;
    int send_size = 0;

/*     void PrintInfo() 
    {
        printf("Number of proces: %d \n", k);
        printf("SEND SIZE: %d \n", send_size);            
        printf("SHAG: %d \n", shag);            

        printf("n: %d \n", n);            
        printf("p: %d \n", p);
        printf("m: %d \n", m);
        //printf("r: %d \n", r);
        //printf("s: %d \n", s);
        printf("l: %d \n", l);
        //printf("K: %d \n", K);
        //printf("Norm: %lf \n", norm);
        printf("nomer_v_okne: %d \n", nomer_v_okne);
        //printf("cur_str: %d \n", cur_str);
        //printf("cpu time: %lf \n", cpu_time);
        //printf("cpu time of all_threads: %lf \n", cpu_time_of_all_threads);
        // printf("astr time: %lf \n",astr_time);

        if (!last_line_isnt_fool)
        {
            printf("Last row is fine\n");
        }
        else
        {
            printf("Last row is bad\n");
        }
    } */
};

double get_cpu_time();
double get_full_time();

int get_max_rows(int n, int m, int p);
int get_rows(int n, int m, int p, int k);
void FormulaMatrixInitialization(double* A, int n, int m, int p, int K, int s);
void PrintLocalMatrix(double* A,  int n, int m, int p, int K, int r);
void PrintAllData(double* A, int n, int m, int p, int K);
int ReadMatrixFromFile(double* A, int n, int m, int p, int K, char* file_name, double* buf, MPI_Comm comm);
void BuildE(double* B, int n, int m, int p, int K);
void PrintMatrix(double* matrix, int n, int m, int p, int K, int r, double* buf, MPI_Comm comm);

#endif
