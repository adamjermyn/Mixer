CC = g++
OPTC = -O3 -fopenmp -ffast-math -mtune=native -fpermissive
SHARED = -fPIC -shared
PROF = $(OPTC) -pg -g
DEBUG = -O0
IEGN = -I ~/eigen/
LIBCUB = ~/cubature-1.0.2/hcubature.so ~/cubature-1.0.2/pcubature.so
ICUB = -I ~/cubature-1.0.2/
OMPENABLE = -fopenmp
TOTALINC = ../Build/integrator.o ../Build/linalg.o ../Build/physics.o
TOTALINCS = ../Build/integrator.so ../Build/linalg.so ../Build/physics.so

all: shared test

shared: python.cpp linalg.so linalg.hpp physics.so physics.hpp integrator.so integrator.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINCS) python.cpp -o ../Build/turb.so $(LIBCUB)

profile: python.cpp linalg.so linalg.hpp physics.so physics.hpp integrator.so integrator.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINCS) $(PROF) python.cpp -o ../Build/turbProf.so $(LIBCUB)

test: test.cpp catch.hpp linalg.o linalg.hpp physics.o physics.hpp integrator.o integrator.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINC) test.cpp -o ../Build/test $(LIBCUB)

integrator.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) -c integrator.cpp -o ../Build/integrator.o

physics.o: physics.cpp physics.hpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) -c physics.cpp -o ../Build/physics.o	

linalg.o: linalg.cpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) -c linalg.cpp -o ../Build/linalg.o

integrator.so: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(OMPENABLE) -c integrator.cpp -o ../Build/integrator.so

physics.so: physics.cpp physics.hpp linalg.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(OMPENABLE) -c physics.cpp -o ../Build/physics.so	

linalg.so: linalg.cpp linalg.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(OMPENABLE) -c linalg.cpp -o ../Build/linalg.so

clean:
	rm ../Build/turb.so
	rm ../Build/test
	rm ../Build/linalg.o
	rm ../Build/integrator.o
	rm ../Build/physics.o
	rm ../Build/linalg.so
	rm ../Build/integrator.so
	rm ../Build/physics.so
