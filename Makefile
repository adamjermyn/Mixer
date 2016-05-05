CC = g++
OPTC = -O3 -fopenmp -ffast-math -mtune=native -fpermissive
PROF = $(OPTC) -pg
DEBUG = -O0
IEGN = -I ~/eigen/
LIBCUB = ~/cubature-1.0.2/hcubature.o
ICUB = -I ~/cubature-1.0.2/
OMPENABLE = -fopenmp

all: main debug profile

main: main.cpp linalg physics
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) linalg.o physics.o main.cpp -o main    $(LIBCUB) 

debug: main.cpp linalg physics
	$(CC) $(DEBUG) $(IEGN) $(ICUB) $(OMPENABLE) linalg.o physics.o main.cpp -o debug   $(LIBCUB)

test: test.cpp linalg physics
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) linalg.o physics.o test.cpp -o test    $(LIBCUB)

profile: main.cpp linalg physics
	$(CC) $(PROF)  $(IEGN) $(ICUB) $(OMPENABLE) linalg.o physics.o main.cpp -o profile $(LIBCUB)


physics: physics.cpp physics.hpp
	$(CC) $(OPTC) $(IEGN) $(OMPENABLE) -c physics.cpp -o physics.o	

linalg: linalg.cpp linalg.hpp
	$(CC) $(OPTC) $(IEGN) $(OMPENABLE) -c linalg.cpp -o linalg.o

clean:
	rm main
	rm debug
	rm profile
	rm test
	rm linalg.o
	rm physics.o