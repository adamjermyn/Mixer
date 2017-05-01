CC = g++
OPTC = -O3 -ffast-math -fpermissive -std=gnu++11 -pipe
SHARED = -fPIC -shared
PROF = $(OPTC) -pg -g
DEBUG = -O0
IEGN = -I ../../eigen/
LIBCUB = ../../cubature-1.0.2/hcubature.so ../../cubature-1.0.2/pcubature.so
ICUB = -I ../../cubature-1.0.2/
TOTALINC = ../Build/integrator.o ../Build/linalg.o ../Build/physics.o ../Builds.basis.o
TOTALINCS = ../Build/integrator.so ../Build/linalg.so ../Build/physics.so ../Builds/basis.so
TOTALINCP = ../Build/integratorP.o ../Build/linalgP.o ../Build/physicsP.o ../Build/basisP.o

all: shared profile test

shared: python.cpp ../Build/linalg.so linalg.hpp ../Build/physics.so physics.hpp ../Build/integrator.so integrator.hpp ../Build/basis.so basis.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(TOTALINCS) python.cpp -o ../Build/turb.so $(LIBCUB)

profile: python.cpp ../Build/linalgP.o linalg.hpp ../Build/physicsP.o physics.hpp ../Build/integratorP.o integrator.hpp ../Build/basisP.o basis.hpp
	$(CC) $(OPTC) $(IEGN) $(ICUB) $(TOTALINCP) $(PROF) profile.cpp -o ../Build/prof.e $(LIBCUB)

test: test.cpp catch.hpp ../Build/linalg.o linalg.hpp ../Build/physics.o physics.hpp ../Build/integrator.o integrator.hpp ../Build/basis.o basis.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(TOTALINC) test.cpp -o ../Build/test.e $(LIBCUB)

../Build/basis.o: basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) -c basis.cpp -o ../Build/basis.o

../Build/integrator.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) -c integrator.cpp -o ../Build/integrator.o

../Build/physics.o: physics.cpp physics.hpp linalg.hpp basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) -c physics.cpp -o ../Build/physics.o	

../Build/linalg.o: linalg.cpp linalg.hpp basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) -c linalg.cpp -o ../Build/linalg.o
 
../Build/basisP.o: basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(PROF) -c basis.cpp -o ../Build/basisP.o

../Build/integratorP.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(PROF) -c integrator.cpp -o ../Build/integratorP.o

../Build/physicsP.o: physics.cpp physics.hpp linalg.hpp basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) $(PROF) -c physics.cpp -o ../Build/physicsP.o	

../Build/linalgP.o: linalg.cpp linalg.hpp basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) $(PROF) -c linalg.cpp -o ../Build/linalgP.o

../Build/basis.so: basis.cpp basis.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) -c basis.cpp -o ../Build/basisP.so

../Build/integrator.so: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) -c integrator.cpp -o ../Build/integrator.so

../Build/physics.so: physics.cpp physics.hpp linalg.hpp basis.cpp basis.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) -c physics.cpp -o ../Build/physics.so	

../Build/linalg.so: linalg.cpp linalg.hpp basis.cpp basis.hpp
	$(CC) $(OPTC) $(SHARED) $(IEGN) -c linalg.cpp -o ../Build/linalg.so

clean:
	rm ../Build/*.so
	rm ../Build/*.o
	rm ../Build/*.e

