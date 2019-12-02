#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#define SIZE 1024

int iterations = 50;
int tolerance = 0;
float* coeff;
float* terms;

void generate_diagonal_dominant_matrix(const ulong size, float min, float max) {
    int i, j;
    for (i = 0; i < size; ++i) {
      for (j = 0; j < size; ++j) {
        coeff[i * size + j] = min + rand() % abs(max - min);
      }
    }
    for (i = 0; i < size; ++i) {
        float sum = 0;
        for (j = 0; j < size; ++j) 
          sum += abs(coeff[i * size + j]);
        sum -= abs(coeff[i * size + i]);
        coeff[i * size + i] = abs(coeff[i* size + i]) + sum;
    }
}
                            
int main(int argc, char **argv) {
    int err, i, j;
    int sum = 0;
    int id, nrproc, n, nrecv;
    int* sendcounts = NULL;
    int* displs = NULL;
    MPI_Init(&argc, &argv);
    err = MPI_Comm_size(MPI_COMM_WORLD, &nrproc);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (id == 0) {
        // Scatter data
        // Each process takes an equal number of columns
        n = SIZE;
        coeff = (float*) malloc(n * n * sizeof(float));
        terms = (float*) malloc(n * sizeof(float));
        for (j = 0; j < n; ++j) {
          terms[j] = -10000 + rand() % 20000;
        }
        generate_diagonal_dominant_matrix(n, -10000, 10000);
        sendcounts = (int *) malloc(nrproc * sizeof(int));
        displs = (int *) malloc(nrproc * sizeof(int));
        for (i = 0; i < nrproc; ++i) {
            sendcounts[i] = (n / nrproc) * n ;
            displs[i] = sum;
            sum += sendcounts[i];
        }
        sendcounts[nrproc - 1] += (n % nrproc) * n;
    }
    // send specific data to processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tolerance, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(terms, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float* recv = (float *) malloc(n * sizeof(float));
    nrecv = (n / nrproc);
    if (id == nrproc - 1)
      nrecv += (n % nrproc);

    MPI_Scatterv(coeff, sendcounts, displs, MPI_FLOAT, 
                recv, nrecv * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float* solutions[2];
    solutions[0] = (float*) calloc(n, sizeof(float));
    solutions[1] = (float*) calloc(nrecv, sizeof(float));

    float zero = tolerance - tolerance;
    float error;

    if (id == 0) {
        for (i = 0; i < nrproc; ++i) {
            sendcounts[i] = (n / nrproc);
            displs[i] = sum;
            sum += sendcounts[i];
        }
        sendcounts[nrproc - 1] += (n % nrproc);
    }

    //Starting iterations
    int iteration;
    for (iteration = 0; iteration < iterations; ++iteration) {
        error = zero;
        int sol_it = iteration % 2;
        int osol_it = (iteration + 1) % 2;
        //calculate solutions
        int start, stop;
        start = (n / nrproc) * id;
        stop = (n / nrproc) * (id + 1);
        if (id == nrproc)
            stop = n;
        for (i = start; i < stop; ++i) {
            for (j = 0; j < n; ++j) {
                terms[i] -= (solutions[osol_it][j] * recv[i * n + j]);
            }
            solutions[sol_it][i] =  (terms[i] + (solutions[osol_it][i] * recv[i * n + i])) / recv[i * n + j];
        }

        MPI_Gatherv(solutions[sol_it] + start, (stop - start + 1), MPI_FLOAT,
                  solutions[sol_it], sendcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(solutions[sol_it], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //Gather solutions

    MPI_Finalize();
}