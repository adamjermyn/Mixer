CC = g++
OPTC = -O3 -fopenmp -ffast-math -mtune=native -fpermissive
PROF = $(OPTC) -pg
DEBUG = -O0
LIBTSG = -L ~/TasmanianSparseGrids/ -ltasmaniansparsegrid
ITSG = -I ~/TasmanianSparseGrids/ -I ~/eigen/
LIBCUB = ~/cubature-1.0.2/hcubature.o
ICUB = -I ~/cubature-1.0.2/
OMPENABLE = -fopenmp

all: main debug profile

main: main.cpp linalg.hpp
	$(CC) $(OPTC)  $(ITSG) $(ICUB) $(OMPENABLE) main.cpp -o main    $(LIBTSG) $(LIBCUB) 

debug: main.cpp linalg.hpp
	$(CC) $(DEBUG) $(ITSG) $(ICUB) $(OMPENABLE) main.cpp -o debug   $(LIBTSG) $(LIBCUB)

test: test.cpp linalg.hpp
	$(CC) $(OPTC)  $(ITSG) $(ICUB) $(OMPENABLE) test.cpp -o test    $(LIBTSG) $(LIBCUB)

profile: main.cpp linalg.hpp
	$(CC) $(PROF)  $(ITSG) $(ICUB) $(OMPENABLE) main.cpp -o profile $(LIBTSG) $(LIBCUB)


clean:
	rm main
	rm debug
	rm profile
	rm test
#  