MAKE_DIR := $(abspath $(lastword $(MAKEFILE_LIST)))
BUILD_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../Build)
EIGEN_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../eigen)
CUBE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../cubature-1.0.2)

CC = g++
OPTC = -O3 -ffast-math -march=native -fpermissive -std=gnu++11 -pipe
SHARED = -fPIC -shared
PROF = $(OPTC) -pg -g
DEBUG = -O0
IEGN = -I $(EIGEN_DIR)
LIBCUB = $(CUBE_DIR)/hcubature.so $(CUBE_DIR)/pcubature.so
ICUB = -I $(CUBE_DIR)/
TOTALINC = $(BUILD_DIR)/integrator.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/physics.o $(BUILD_DIR)/basis.o
TOTALINCS = $(BUILD_DIR)/integrator.so $(BUILD_DIR)/linalg.so $(BUILD_DIR)/physics.so $(BUILD_DIR)/basis.so
TOTALINCP = $(BUILD_DIR)/integratorP.o $(BUILD_DIR)/linalgP.o $(BUILD_DIR)/physicsP.o $(BUILD_DIR)/basisP.o

all: shared profile test

shared: $(TOTALINCS)
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(TOTALINCS) python.cpp -o $(BUILD_DIR)/turb.so $(LIBCUB)

profile: $(TOTALINCP)
	$(CC) $(OPTC) $(IEGN) $(ICUB) $(TOTALINCP) $(PROF) profile.cpp -o $(BUILD_DIR)/prof.e $(LIBCUB)

test: $(TOTALINC)
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(TOTALINC) test.cpp -o $(BUILD_DIR)/test.e $(LIBCUB)

../Build/basis.o: basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) -c basis.cpp -o $(BUILD_DIR)/basis.o

../Build/integrator.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) -c integrator.cpp -o $(BUILD_DIR)/integrator.o

../Build/physics.o:
	$(CC) $(OPTC)  $(IEGN) -c physics.cpp -o $(BUILD_DIR)/physics.o	

../Build/linalg.o:
	$(CC) $(OPTC)  $(IEGN) -c linalg.cpp -o $(BUILD_DIR)/linalg.o
 
../Build/basisP.o:
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(PROF) -c basis.cpp -o $(BUILD_DIR)/basisP.o

../Build/integratorP.o:
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(PROF) -c integrator.cpp -o $(BUILD_DIR)/integratorP.o

../Build/physicsP.o:
	$(CC) $(OPTC)  $(IEGN) $(PROF) -c physics.cpp -o $(BUILD_DIR)/physicsP.o	

../Build/linalgP.o:
	$(CC) $(OPTC)  $(IEGN) $(PROF) -c linalg.cpp -o $(BUILD_DIR)/linalgP.o

../Build/basis.so:
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) -c basis.cpp -o $(BUILD_DIR)/basis.so

../Build/integrator.so:
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) -c integrator.cpp -o $(BUILD_DIR)/integrator.so

../Build/physics.so:
	$(CC) $(OPTC) $(SHARED) $(IEGN) -c physics.cpp -o $(BUILD_DIR)/physics.so	

../Build/linalg.so:
	$(CC) $(OPTC) $(SHARED) $(IEGN) -c linalg.cpp -o $(BUILD_DIR)/linalg.so

clean:
	rm ../Build/*.so
	rm ../Build/*.o
	rm ../Build/*.e

