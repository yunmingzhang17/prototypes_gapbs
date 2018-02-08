# See LICENSE.txt for license details.
CXX=icpc -I /data/scratch/mengjiao/graphit/src/runtime_lib/ -DOPENMP

INCLUDE_DIR=./include/

ifeq ($(DEBUG), 1)
CXX_FLAGS += -std=c++11 -O0 -g -Wall -DDEBUG -march=native -mtune=native
PAR_FLAG = -fopenmp

else 

CXX_FLAGS += -std=c++11 -O3 -Wall -g
PAR_FLAG = -fopenmp

endif

ifeq ($(TIME_MSG), 1)
CXX_FLAGS += -DTIME_MSG
endif

ifeq ($(LOAD_MSG), 1)
CXX_FLAGS += -DLOAD_MSG
endif

ifeq ($(COUNT), 1)
CXX_FLAGS += -DCOUNT
endif

ifdef NUMA
CXX_FLAGS += -lnuma
endif

ifeq ($(PROFILE), 1)
CXX_FLAGS += -DPROFILE
endif

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -qopenmp
endif

ifneq (,$(findstring sunCC,$(CXX)))
	CXX_FLAGS = -std=c++11 -xO3 -m64 -xtarget=native
	PAR_FLAG = -xopenmp
endif

ifneq ($(SERIAL), 1)
	CXX_FLAGS += $(PAR_FLAG)
endif


CXX_FLAGS += -I ${INCLUDE_DIR}

KERNELS = bc bfs cc pr sssp tc tc_migra multi_hop_simple_rec numa_cc numa_pr_delta numa_cf pr_cache
SUITE = $(KERNELS) converter

.PHONY: all
all: $(SUITE)

% : src/%.cc src/*.h
	$(CXX) $(CXX_FLAGS) $< -o $@

# Testing
include test/test.mk

.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
