build:
	g++ --std=c++11 main.cpp -fopenmp -lpthread -o jacobi

clean:
	rm -rf jacobi
