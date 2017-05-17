MAKE_DIR := $(abspath $(lastword $(MAKEFILE_LIST)))
BUILD_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../Build)
EIGEN_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../eigen)
CUBE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../../cubature-1.0.2)

OS := $(shell uname)
ifeq ($(OS),Darwin)
	CC = /usr/local/bin/g++-7
	LIBCUB = $(CUBE_DIR)/hcubature_mac.so
	OPTC = -O3 -ffast-math -fpermissive -std=gnu++11 -pipe
else	
	CC = g++
	LIBCUB = $(CUBE_DIR)/hcubature.so
	OPTC = -O3 -ffast-math -fpermissive -std=gnu++11 -pipe
endif

SHARED = -shared -fPIC
PROF = $(OPTC) -pg -g
DEBUG = -O0
IEGN = -I $(EIGEN_DIR)
ICUB = -I $(CUBE_DIR)/
TOTALINC = $(BUILD_DIR)/integrator.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/physics.o $(BUILD_DIR)/basis.o
TOTALINCS = $(BUILD_DIR)/integrator.so $(BUILD_DIR)/linalg.so $(BUILD_DIR)/physics.so $(BUILD_DIR)/basis.so
TOTALINCP = $(BUILD_DIR)/integratorP.o $(BUILD_DIR)/linalgP.o $(BUILD_DIR)/physicsP.o $(BUILD_DIR)/basisP.o

all: shared profile test test_linal test_MRI test_ROT

shared: $(TOTALINCS)
	$(CC) $(OPTC) $(SHARED) $(IEGN) $(ICUB) $(TOTALINCS) python.cpp -o $(BUILD_DIR)/turb.so $(LIBCUB)

profile: $(TOTALINCP)
	$(CC) $(OPTC) $(IEGN) $(ICUB) $(TOTALINCP) $(PROF) profile.cpp -o $(BUILD_DIR)/prof.e $(LIBCUB)

test: $(TOTALINC)
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(TOTALINC) test.cpp -o $(BUILD_DIR)/test.e $(LIBCUB)

test_linalg: $(BUILD_DIR)/linalg.o
	$(CC) $(OPTC) $(BUILD_DIR)/linalg.o test_linalg.cpp -o $(BUILD_DIR)/test_linalg.e

test_MRI: $(BUILD_DIR)/linalg.o $(BUILD_DIR)/physics.o $(BUILD_DIR)/basis.o
	$(CC) $(OPTC) $(IEGN) $(BUILD_DIR)/linalg.o $(BUILD_DIR)/physics.o $(BUILD_DIR)/basis.o test_MRI.cpp -o $(BUILD_DIR)/test_MRI.e

test_ROT: $(BUILD_DIR)/linalg.o $(BUILD_DIR)/physics.o $(BUILD_DIR)/basis.o
	$(CC) $(OPTC) $(IEGN) $(BUILD_DIR)/linalg.o $(BUILD_DIR)/physics.o $(BUILD_DIR)/basis.o test_rotational_shear.cpp -o $(BUILD_DIR)/test_ROT.e

$(BUILD_DIR)/basis.o: basis.cpp basis.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) -c basis.cpp -o $(BUILD_DIR)/basis.o

$(BUILD_DIR)/integrator.o: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(ICUB) -c integrator.cpp -o $(BUILD_DIR)/integrator.o

$(BUILD_DIR)/physics.o:
	$(CC) $(OPTC)  $(IEGN) -c physics.cpp -o $(BUILD_DIR)/physics.o	

$(BUILD_DIR)/linalg.o:
	$(CC) $(OPTC)  $(IEGN) -c linalg.cpp -o $(BUILD_DIR)/linalg.o
 
$(BUILD_DIR)/integrator.so: integrator.cpp integrator.hpp linalg.hpp physics.hpp
	$(CC) $(OPTC)  $(IEGN) $(SHARED) $(ICUB) -c integrator.cpp -o $(BUILD_DIR)/integrator.so

$(BUILD_DIR)/physics.so:
	$(CC) $(OPTC)  $(IEGN) $(SHARED) -c physics.cpp -o $(BUILD_DIR)/physics.so	

$(BUILD_DIR)/linalg.so:
	$(CC) $(OPTC)  $(IEGN) $(SHARED) -c linalg.cpp -o $(BUILD_DIR)/linalg.so
 
$(BUILD_DIR)/basis.so:
	$(CC) $(OPTC)  $(IEGN) $(SHARED) $(ICUB) -c basis.cpp -o $(BUILD_DIR)/basis.so

$(BUILD_DIR)/basisP.o:
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(PROF) -c basis.cpp -o $(BUILD_DIR)/basisP.o

$(BUILD_DIR)/integratorP.o:
	$(CC) $(OPTC)  $(IEGN) $(ICUB) $(PROF) -c integrator.cpp -o $(BUILD_DIR)/integratorP.o

$(BUILD_DIR)/physicsP.o:
	$(CC) $(OPTC)  $(IEGN) $(PROF) -c physics.cpp -o $(BUILD_DIR)/physicsP.o	

$(BUILD_DIR)/linalgP.o:
	$(CC) $(OPTC)  $(IEGN) $(PROF) -c linalg.cpp -o $(BUILD_DIR)/linalgP.o

clean:
	rm ../Build/*

