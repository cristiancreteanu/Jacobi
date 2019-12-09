#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
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
    double stime, etime;
    int err, i, j;
    int sum = 0;
    int id, nrproc, n, nrecv;
    int* sendcounts = NULL;
    int* displs = NULL;
    MPI_Init(&argc, &argv);
    err = MPI_Comm_size(MPI_COMM_WORLD, &nrproc);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    n = SIZE;
    sendcounts = (int *) malloc(nrproc * sizeof(int));
    displs = (int *) malloc(nrproc * sizeof(int));
    terms = (float*) malloc(n * sizeof(float));
    
    // Scatter data
    // Each process takes an equal number of columns
    coeff = (float*) malloc(n * n * sizeof(float));
    for (j = 0; j < n; ++j) {
      terms[j] = -10000 + rand() % 20000;
    }
    generate_diagonal_dominant_matrix(n, -10000, 10000);
    for (i = 0; i < nrproc; ++i) {
        sendcounts[i] = (n / nrproc) * n ;
        displs[i] = sum;
        sum += sendcounts[i];
    }
    sendcounts[nrproc - 1] += (n % nrproc) * n;
    
    printf("Init done\n");
    // send specific data to processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tolerance, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(terms, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    printf("Bcast done\n");
    float* recv = (float *) calloc(sendcounts[id], sizeof(float));
    nrecv = (n / nrproc);
    if (id == nrproc - 1)
      nrecv += (n % nrproc);

    printf("Scatterv\n");
    MPI_Scatterv(coeff, sendcounts, displs, MPI_FLOAT, 
                recv, nrecv * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    printf("Scatterv done\n");
    float* solutions[2];
    solutions[0] = (float*) calloc(nrecv, sizeof(float));
    solutions[1] = (float*) calloc(nrecv, sizeof(float));

    float zero = tolerance - tolerance;
    float error;

    printf("Start algo\n");
    stime = MPI_Wtime();
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
            float term = terms[i];
            int idx = i - start;
            for (j = 0; j < n; ++j) {
                term -= (solutions[osol_it][j] * recv[idx * n + j]);
            }
            solutions[sol_it][i] =  (term + (solutions[osol_it][i] * recv[idx * n + i])) / recv[idx * n + j];
        }
        printf("Gather\n");
        float *p = solutions[sol_it] + start;
        
        MPI_Gatherv(p, (stop - start + 1), MPI_FLOAT, solutions[sol_it], 
                sendcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
        printf("Broadcast\n");
        MPI_Bcast(solutions[sol_it], n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        printf("Barrier");
        MPI_Barrier(MPI_COMM_WORLD);
    }
    etime = MPI_Wtime();
    //Gather solutions
    printf("Total time is: %lf", etime - stime);
    MPI_Finalize();
    return 0;
}