CC = g++
OPTC = -O3 -fopenmp -ffast-math -mtune=native
DEBUG = -O0
LIBTSG = -L ~/TasmanianSparseGrids/ -ltasmaniansparsegrid
ITSG = -I ~/TasmanianSparseGrids/ -I ~/eigen/
OMPENABLE = -fopenmp

all: main debug profile

main: main.cpp linalg.hpp
	$(CC) $(OPTC) $(ITSG) $(OMPENABLE) main.cpp -o main $(LIBTSG)

debug: main.cpp linalg.hpp
	$(CC) $(ITSG) $(OMPENABLE) main.cpp -o debug $(LIBTSG)

profile: main.cpp linalg.hpp
	$(CC) $(ITSG) $(OMPENABLE) -pg main.cpp -o profile $(LIBTSG)

test: test.cpp linalg.hpp
	$(CC) $(OPTC) $(ITSG) $(OMPENABLE) test.cpp -o test $(LIBTSG)

clean:
	rm main
	rm debug
	rm profile
	rm test
