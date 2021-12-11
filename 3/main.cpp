#include <chrono>
#include <stdio.h>
#include <vector>
#include <cmath>

#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"

double f(double x) {
    return x;
}

int main(int argc, char **argv) {
    int current_rank = 0;
    int p = 0;
    const int n1 = 8;
    const int n2 = 4;
    int r2 = 2;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
    int r1 = n1 / p;
    MPI_Datatype block_type;
    MPI_Type_vector(r2, r1, n2, MPI_INT, &block_type);
    MPI_Type_commit(&block_type);
    MPI_Datatype block_line_type;
    MPI_Type_contiguous(r1, MPI_INT, &block_line_type);
    MPI_Type_commit(&block_line_type);
    if (current_rank == 0) {
        int a[n1][n2] = {{5, 7, 10, 3,}};
        for (int block_num = 0; block_num < n2 / r2; block_num++) {
            for (int i = 1; i < r1; i++) {
                for (int j = 0; j < r2; j++) {
                    a[i][block_num * r2 + j] = a[i - 1][block_num * r2 +  j] + 1;
                }
            }
            MPI_Send(&a[r1 - 1][block_num * r2], 1, block_line_type, 1, 1, MPI_COMM_WORLD);
        }
        for (int dest_rank = 1; dest_rank < p; ++dest_rank) {
            for (int block_num = 0; block_num < n2 / r2; block_num++) {
                MPI_Recv(&a[dest_rank * r1][block_num * r2], 1, block_type, dest_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        for (auto & row : a) {
            for (int el : row){
                std::cout << el << "\t";
            }
            std::cout << "\n";
        }
    } else {
        int a[r1 + 1][r2];
        for (int block_num = 0; block_num < n2 / r2; block_num++) {
            MPI_Recv(a[0], 1, block_line_type, current_rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 1; i < r1 + 1; i++) {
                for (int j = 0; j< r2; j++) {
                    a[i][j] = a[i - 1][j] + 1;
                }
            }
            if (current_rank != p - 1) {
                MPI_Send(a[r1], 1, block_line_type, current_rank + 1, 1, MPI_COMM_WORLD);
            }
            MPI_Send(&a[1][0], n2 / r2, block_line_type, 0, 1, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
}