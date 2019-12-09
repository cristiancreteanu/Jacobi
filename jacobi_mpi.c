#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#define SIZE 4096

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
    if (id == 0) {
        // Scatter data
        // Each process takes an equal number of columns
        coeff = (float*) malloc(n * n * sizeof(float));
        for (j = 0; j < n; ++j) {
          terms[j] = -10000 + rand() % 20000;
        }
        generate_diagonal_dominant_matrix(n, -10000, 10000);
    }

    for (i = 0; i < nrproc; ++i) {
        sendcounts[i] = (n / nrproc) * n ;
        displs[i] = sum;
        sum += sendcounts[i];
    }
    sendcounts[nrproc - 1] += (n % nrproc) * n;
    // send specific data to processes
    stime = MPI_Wtime();
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tolerance, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(terms, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float* recv = (float *) calloc(sendcounts[id], sizeof(float));
    nrecv = (n / nrproc);
    if (id == nrproc - 1)
      nrecv += (n % nrproc);

    MPI_Scatterv(coeff, sendcounts, displs, MPI_FLOAT, 
                recv, nrecv * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    float* solutions;
    float* old_solutions;
    solutions = (float*) calloc(n, sizeof(float));
    old_solutions = (float*) calloc(n, sizeof(float));

    float zero = tolerance - tolerance;
    float error;
    sum = 0;
    for (i = 0; i < nrproc; ++i) {
      sendcounts[i] = (n / nrproc);
      displs[i] = sum;
      sum += sendcounts[i];
    }
    sendcounts[nrproc - 1] += (n % nrproc);

    
    //Starting iterations
    int iteration;
    for (iteration = 0; iteration < iterations; ++iteration) {
        error = zero;
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
                term -= (old_solutions[j] * recv[idx * n + j]);
            }
            solutions[i] =  (term + (old_solutions[i] * recv[idx * n + i])) / recv[idx * n + i];
        }

        float *p = solutions + start;
        if (id != 0) {
            MPI_Send(p, stop - start, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        } else {
            for (i = 1; i < nrproc; ++i) {
                start = (n / nrproc) * i;
                MPI_Recv(old_solutions + start, n / nrproc, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("%d\n", id);
            }
        }
        printf("%d\n", id);
        
        MPI_Bcast(old_solutions, n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    etime = MPI_Wtime();

    printf("Total time is: %lf\n", etime - stime);
    MPI_Finalize();
    return 0;
}