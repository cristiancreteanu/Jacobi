#include <iostream>
#include <cstring>
#include "jacobi.hpp"
#include <unistd.h>

using namespace std;


vector<vector<float>> matrix __attribute__((aligned(64)));
vector<vector<float>> array __attribute__((aligned(64)));
vector<float> terms __attribute__((aligned(64)));


void parse_args(int argc, char *const argv[]);

ulong size = 1024;
ulong workers = 8;
ulong iterations = 50;
float tolerance = 0;

char *filename;

ofstream outfile;
auto to_csv = false;
auto debug = false;
auto flag = false;
unsigned int seed = 42;


int main(const int argc, char *const argv[]) {

    parse_args(argc, argv);

    srand((int) seed);

    vector<float> solution;
    generate_diagonal_dominant_matrix(size, matrix, -range, range);
    generate_vector(size, terms, -range, range);
    cout << matrix_size << ' ' << size << endl;
    cout << algorithm << ' ';

    switch (method) {
        case SEQUENTIAL:
            cout << sequential_string << endl;
            solution = serial_jacobi(matrix, terms, iterations, tolerance);
            break;
        case OPENMP:
            cout << openmp_string << endl;
            solution = openmp_jacobi(matrix, terms, iterations, tolerance);
        case PTHREADS:
            cout << pthreads_string << endl;
            solution = pthread_jacobi(matrix, terms, iterations, tolerance, workers);
    }

    if (debug) {
        cout << "solution: ";
        print_solution(solution);
        float error = check_error(matrix, terms, solution);
        cout << "error: " << error << endl;
    }

    return 0;
}


void print_helper() {
    cout << "Usage: " << "main " << "<algorithm> " << endl << endl;
    cout << "The required arguments are:" << endl;
    cout << '\t' << "algorithm " << "indicates the algorithm executed taken from the following list" << endl;
    cout << "\t\t" << sequential_string << ": sequential jacobi algorithm" << endl;
    cout << "\t\t" << openmp_string << ": openmp jacobi algorithm" << endl;
    if (flag) exit(0);
    exit(EINVAL);
}

void parse_args(const int argc, char *const argv[]) {
    if (argc < 2) {
        print_helper();
    }

    string arg = std::string(argv[1]);

    if (arg == sequential_string) method = SEQUENTIAL;
        else if (arg == openmp_string) method == OPENMP;
            else print_helper();
    
    errno = 0;
    if (errno) {
        print_helper();
    }
}

