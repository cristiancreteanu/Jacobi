

#ifndef JACOBI_JACOBI_H
#define JACOBI_JACOBI_H

#include <vector>
#include "utils.hpp"
#include <iostream>
#include <pthread.h>

using namespace std;

/**
 * Serial implementation of the Jacobi method simply iterates until reach the convergence or reach the max number of
 * iterations
 * @tparam T template, type of the values
 * @param coefficients coefficients matrix of the linear system
 * @param terms right vector of the system
 * @param iterations max number of iteration
 * @param tolerance error tolerated
 * @return solution vector
 */
template<typename T>
std::vector<T> serial_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                             const ulong iterations, const T tolerance) {

    start_time = Time::now();
    //allocate solution vectors
    std::vector<T> old_solutions __attribute__((aligned(64)));
    std::vector<T> solutions __attribute__((aligned(64)));

    T zero = tolerance - tolerance;
    T error;

    //initialize solution vectors
    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(zero);
        solutions.emplace_back(zero);
    }
    //Starting iterations
    init_time = Time::now();
    ulong iteration;
    for (iteration = 0; iteration < iterations; ++iteration) {
        error = zero;
        //calculate solutions
        for (ulong i = 0; i < solutions.size(); ++i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);

        }
        //compute the error
        for (ulong i = 0; i < solutions.size(); ++i)
            error += std::abs(solutions[i] - old_solutions[i]);

        // check the error
        error /= solutions.size();
        if (error <= tolerance) break;
        swap(solutions, old_solutions);
    }
    total_time = Time::now();
    print_metrics(iteration, error);

    return solutions;
}

template<typename T>
std::vector<T> openmp_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                             const ulong iterations, const T tolerance) {

    start_time = Time::now();
    //allocate solution vectors
    std::vector<T> solutions[2] __attribute__((aligned(64)));

    T zero = tolerance - tolerance;
    T error;

    //initialize solution vectors
    #pragma omp parallel for
    for (int i = 0; i < coefficients.size(); ++i) {
        solutions[0].emplace_back(zero);
        solutions[1].emplace_back(zero);
    }


    //Starting iterations
    init_time = Time::now();
    ulong iteration;
    for (iteration = 0; iteration < iterations; ++iteration) {
        error = zero;
        //calculate solutions
        #pragma omp parallel for
        for (ulong i = 0; i < solutions[0].size(); ++i) {
            solutions[(iteration + 1)%2][i] = solution_find(coefficients[i], solutions[iteration%2], terms[i], i);
        }

        //compute the error
        #pragma omp parallel for
        for (ulong i = 0; i < solutions[0].size(); ++i)
            error += std::abs(solutions[(iteration + 1)%2][i] - solutions[iteration%2][i]);

        // check the error
        error /= solutions[0].size();
        if (error <= tolerance) break;
    }


    total_time = Time::now();
    print_metrics(iteration, error);

    return solutions[iteration%2];
}
/*
template<typename T>
struct jacobi_data {
    std::vector<std::vector<T>> &coeff;
    std::vector<T> &terms;
    std::vector<T>&solutions;
    std::vector<T>&old_solutions;
    ulong iterations;
    ulong tolerance;
    int id;
    jacobi_data(std::vector<std::vector<T>> &coeff, std::vector<T> &terms, ulong iterations
                    ulong tolerance, int id, std::vector<T>&solutions, std::vector<T>&old_solutions) {
        this->coeff = coeff;
        this->terms = terms;
        this->iterations = iterations;
        this->tolerance = tolerance;
        this->id = id;
        this->solutions = solutions;
        this->old_solutions = old_solutions;
    }
};*/

/*
template<typename T> thread_calc(void* data) {
    jacobi_data* x = (jacobi_data *) data;
    ulong iteration;
    for (iteration = 0; iteration < x->iterations; ++iteration) {
        error = zero;
        //calculate solutions
        int start = TO DO 
        int stop = TO DO
        for (ulong i = start; i < stop; ++i) {
            x->solutions[i] = solution_find(x->coeff[i], x->old_solutions, x->terms[i], i);
        }
        //compute the error
        for (ulong i = start; i < stop; ++i)
            error += std::abs(x->solutions[i] - x->old_solutions[i]);

        // check the error
        error /= (stop - start + 1);
        if (error <= tolerance) break;

        for (ulong i = start; i < stop; ++i)
            swap(solutions[i], old_solutions[i]);
        //put barrier here
    }

}
*/

template<typename T>
std::vector<T> pthreads_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                             const ulong iterations, const T tolerance, int numthread) {
/*
    start_time = Time::now();
    //allocate solution vectors
    std::vector<T> old_solutions __attribute__((aligned(64)));
    std::vector<T> solutions __attribute__((aligned(64)));

    T zero = tolerance - tolerance;
    T error;

    //initialize solution vectors
    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(zero);
        solutions.emplace_back(zero);
    }
    init_time = Time::now();
    for (int i = 0; i < numthread; ++i) {
        struct jacobi_data data(coefficients, terms, iterations, tolerance, i, solutions, old_solutions);
        pthread_create(&threads[i], NULL, ,(void *) data)
    }

    total_time = Time::now();

    return solutions;
    */
    std::vector<T>sol;
    return sol;
}



/*
Need to port this function to C
std::vector<float> mpi_jacobi(const std::vector<float> coefficients, const std::vector<float> terms,
                             const ulong iterations, const T tolerance, std::ofstream &out) {
    int err;
    int id, nrproc, n, nrecv;
    int* sendcounts = NULL;
    int* displs = NULL;
    MPI_Init();
    err = MPI_Comm_size(MPI_COMM_WORLD, &nrproc);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (id == 0) {
        // Scatter data
        // Each process takes an equal number of columns
         n = terms.size();
        sendcounts = new int[nrproc];
        displs = new int[nrproc];
        int sum = 0;
        for (int i = 0; i < nrproc; ++i) {
            sendcounts[i] = (n / nrproc) * n ;
            displs[i] = sum;
            sum += sendcounts[i];
        }
        sendcounts[nrproc - 1] = (n % nrproc) * n;
    
    }
    // send specific data to processes
    MPI_Bcast(&n, 1, MPI_INT, , MPI_COMM_WORLD);
    MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tolerance, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(terms.data(), terms.size(), MPI_FLOAT, 0, MPI_COMM_WORLD);

    vector<float>recv(n, 0);
    nrecv = (n / nrproc);
    if (id == nrproc - 1)
        nrecv += (n % nrproc);
    MPI_Scatterv(coefficients.data(), sendcounts, displs, MPI_FLOAT, 
                recv.data(), nrecv * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

     Algorithm 

    std::vector<float> old_solutions(nrecv, 0) __attribute__((aligned(64)));
    std::vector<float> solutions(nrecv, 0) __attribute__((aligned(64)));

    float zero = tolerance - tolerance;
    float error;

    //Starting iterations
    ulong iteration;
    for (iteration = 0; iteration < iterations; ++iteration) {
        error = zero;
        //calculate solutions
        for (ulong i = 0; i < solutions.size(); ++i) {
            for (int j = 0; j < n; ++j) {
                terms[i] -= (old_solutions[j] * recv[i * n + j]);
            }
            solutions[i] =  (terms[i] + (old_solutions[i] * recv[i * n + i])) / recv[i * n + j];
        }

        //compute the error
        for (ulong i = 0; i < solutions.size(); ++i)
            error += std::abs(solutions[i] - old_solutions[i]);

        // check the error
        error /= solutions.size();
        if (error <= tolerance) break;
        swap(solutions, old_solutions);
    }

    //Gather solutions

    MPI_Finalize();
}*/


#endif //JACOBI_JACOBI_H
