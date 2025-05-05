#include "mpi.h"

#include <cstdio>
#include <cstring>
#define len 30

int main(int argc, char *argv[]) {
    int p = argc, k;
    printf("%d", p);
    MPI_Comm com = MPI_COMM_WORLD;
    int tag = 0;
    char buf[len];
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(com, &p); // количество потоков
    MPI_Comm_rank(com, &k); // номер потока, который запустил этот процесс
    snprintf(buf, len, "Hello from process %d!", k); // в буфере образуется строка, каждый поток формирует строку
    if (k != 0) { // не главный процесс отправляет строку главному процессу
        MPI_Send(buf, strlen(buf) + 1, MPI_CHAR, 0, tag, com);
    }
    else { // главный поток делает Hello от себя
        printf("%s\n", buf);
        for (int i = 1; i < p; i++) {
            MPI_Recv(buf, len, MPI_CHAR, i /* тот поток, который отправил свою строку */, tag, com, &status);
            printf("%s\n", buf);
        }
    }
    MPI_Finalize();
    return 0;
}