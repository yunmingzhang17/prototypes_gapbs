# See LICENSE.txt for license details.
CXX=

INCLUDE_DIR=./include/

ifeq ($(DEBUG), 1)
CXX_FLAGS += -std=c++11 -O0 -g -Wall -DDEBUG
PAR_FLAG = -fopenmp

else 

CXX_FLAGS += -std=c++11 -O3 -Wall -g
PAR_FLAG = -fopenmp

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

KERNELS = bc bfs cc pr sssp sssp_benchmark sssp_bucket_opt tc tc_migra multi_hop_simple_rec cc_sv kcore ppsp ppsp_benchmark
SUITE = $(KERNELS) converter

.PHONY: all
all: $(SUITE)

% : src/%.cc src/*.h
	$(CXX) $(CXX_FLAGS) $< -o $@

# Testing
include test/test.mk

# Benchmark Automation
include benchmark/bench.mk


.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
