#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"

int main(int argc, char **argv) {
    int tag = 55;
    int myrank = 0;
    int p = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank != 0) {
        int dest = 0;
        const char* msg = ("message " + std::to_string(myrank)).c_str();
        std::cout << "Hello from process " << myrank << std::endl;
        MPI_Send(msg, strlen(msg) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    } else {
        char* message[p];
        for (int i = 1; i < p; i++) {
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            int size;
            MPI_Get_count(&status, MPI_CHAR, &size);
            message[status.MPI_SOURCE - 1] = new char[status.MPI_SOURCE - 1];
            MPI_Recv(message[status.MPI_SOURCE - 1], size, MPI_CHAR, status.MPI_SOURCE, tag, MPI_COMM_WORLD,
                     &status);
        }
        for (int i = 0; i + 1 < p; ++i) {
            std::cout << message[i] << std::endl;
        }
    }
    MPI_Finalize();
}