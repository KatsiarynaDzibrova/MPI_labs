#include <chrono>
#include <stdio.h>
#include <vector>
#include <cmath>

#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"

double f(double x) {
    return x;
}

int main(int argc, char **argv) {
    int ping_tag = 55;
    int send_data_tag = 60;
    int end_data_tag = 61;
    int get_res_tag = 65;
    int current_rank = 0;
    int p = 0;
    int n = 10000;
    int r = 5;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
    if (current_rank != 0) {
        int dest = 0;
        while (true) {
            MPI_Send(&current_rank, 1, MPI_INT, dest, ping_tag, MPI_COMM_WORLD);
            double process_data[3];
            MPI_Recv(&process_data, 3, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            double begin = process_data[0];
            double end = process_data[1];
            double step = process_data[2];
            if (std::fabs(end - begin) < 0.0001) {
                break;
            }
            double res = 0;
            int iterations = (end - begin) / step;
            for (int i = 0; i < iterations; ++i) {
                double x1 = begin + step * i;
                double x2 = begin + step * (i + 1);
                res += (f(x1) + f(x2)) * (x2 - x1) / 2;
            }
            MPI_Send(&res, 1, MPI_DOUBLE, dest, get_res_tag, MPI_COMM_WORLD);
        }
    } else {
        double a = 0;
        double b = 2000;
        int blocks_number = n / r;
        double step = (b - a) / n;
        for (int i = 0; i < blocks_number; ++i) {
            int dest_rank;
            MPI_Recv(&dest_rank, 1, MPI_INT, MPI_ANY_SOURCE, ping_tag, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            double begin = a + step * r * i;
            double end = std::min(a + step * r * (i + 1), b);
            double process_data[3] = {begin, end, step};
            MPI_Send(&process_data, 3, MPI_DOUBLE, dest_rank, send_data_tag, MPI_COMM_WORLD);

        }
        for (int dest_rank = 1; dest_rank < p; ++dest_rank) {
            double end_msg[3] = {0, 0, 0};
            MPI_Send(&end_msg, 3, MPI_DOUBLE, dest_rank, end_data_tag, MPI_COMM_WORLD);
        }
        double res = 0;
        for (int i = 0; i < blocks_number; ++i) {
            double proc_res;
            MPI_Recv(&proc_res, 1, MPI_DOUBLE, MPI_ANY_SOURCE, get_res_tag, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            res += proc_res;
        }
        std::cout << "final result: " << res << std::endl;
    }
    MPI_Finalize();
}