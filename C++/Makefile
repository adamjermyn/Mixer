CC = g++
OPTC = -O3 -fopenmp -ffast-math -fpermissive -std=gnu++11 -pipe
SHARED = -fPIC -shared
PROF = $(OPTC) -pg -g
DEBUG = -O0
IEGN = -I ~/eigen/
LIBCUB = ~/cubature-1.0.2/hcubature.so ~/cubature-1.0.2/pcubature.so
ICUB = -I ~/cubature-1.0.2/
OMPENABLE = -fopenmp
TOTALINC = ../Build/integrator.o ../Build/linalg.o ../Build/physics.o
TOTALINCS = ../Build/integrator.so ../Build/linalg.so ../Build/physics.so
TOTALINCP = ../Build/integratorP.o ../Build/linalgP.o ../Build/physicsP.o

all: shared profile test

shared: python.cpp ../Build/linalg.so linalg.hpp ../Build/physics.so physics.hpp ../Build/integrator.so integrator.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINCS) python.cpp -o ../Build/turb.so $(LIBCUB)

profile: python.cpp ../Build/linalgP.o linalg.hpp ../Build/physicsP.o physics.hpp ../Build/integratorP.o integrator.hpp
	$(CC) $(OPTC) $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINCP) $(PROF) profile.cpp -o ../Build/prof.e $(LIBCUB)

test: test.cpp catch.hpp ../Build/linalg.o linalg.hpp ../Build/physics.o physics.hpp ../Build/integrator.o integrator.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) $(TOTALINC) test.cpp -o ../Build/test.e $(LIBCUB)

../Build/integrator.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) -c integrator.cpp -o ../Build/integrator.o

../Build/physics.o: physics.cpp physics.hpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) -c physics.cpp -o ../Build/physics.o	

../Build/linalg.o: linalg.cpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) -c linalg.cpp -o ../Build/linalg.o

../Build/integratorP.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(OMPENABLE) $(PROF) -c integrator.cpp -o ../Build/integratorP.o

../Build/physicsP.o: physics.cpp physics.hpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) $(PROF) -c physics.cpp -o ../Build/physicsP.o	

../Build/linalgP.o: linalg.cpp linalg.hpp
	$(CC) $(OPTC)  $(IEGN) $(OMPENABLE) $(PROF) -c linalg.cpp -o ../Build/linalgP.o

../Build/integrator.so: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(OMPENABLE) -c integrator.cpp -o ../Build/integrator.so

../Build/physics.so: physics.cpp physics.hpp linalg.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(OMPENABLE) -c physics.cpp -o ../Build/physics.so	

../Build/linalg.so: linalg.cpp linalg.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(OMPENABLE) -c linalg.cpp -o ../Build/linalg.so

clean:
	rm ../Build/*.so
	rm ../Build/*.o
	rm ../Build/*.e

