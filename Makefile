CC = g++
OPTC = -O3 -fopenmp -ffast-math -mtune=native -fpermissive -fPIC -shared
PROF = $(OPTC) -pg -g
DEBUG = -O0
IEGN = -I ~/eigen/
LIBCUB = ~/cubature-1.0.2/hcubature.so ~/cubature-1.0.2/pcubature.so
ICUB = -I ~/cubature-1.0.2/
OMPENABLE = -fopenmp
TOTALINC = integrator.o linalg.o physics.o

all: shared test profile

profile: main.cpp linalg.o linalg.hpp physics.o physics.hpp integrator.o integrator.hpp
	$(CC) $(PROF)  $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINC) main.cpp -o profile $(LIBCUB)

test: test.cpp catch.hpp linalg.o linalg.hpp physics.o physics.hpp integrator.o integrator.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINC) test.cpp -o test $(LIBCUB)

integrator.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) -c integrator.cpp -o integrator.o

physics.o: physics.cpp physics.hpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) -c physics.cpp -o physics.o	

linalg.o: linalg.cpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) -c linalg.cpp -o linalg.o

clean:
	rm main
	rm debug
	rm profile
	rm test
	rm linalg.o
	rm integrator.o
	rm physics.o
