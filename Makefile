# Copyright (C) 2022 Iver Håkonsen
#
# cardiac-fibers is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
#
# Authors: Iver Håkonsen <hakonseniver@yahoo.no


CFLAGS = -std=c++11
ifneq ($(DEBUG), YES)
CFLAGS += -O3 -Wall
else
CFLAGS += -O0 -DDEBUG -g3 -Wall
endif

MPI_CXX ?= mpic++
CXX = $(MPI_CXX)

# MFEM
ifeq ($(DEBUG), YES)
# Set the debug build to MFEM_ROOT if MFEM_DBG_ROOT is not set
ifdef MFEM_DBG_ROOT
MFEM_ROOT=$(MFEM_DBG_ROOT)
endif
endif

MFEM_INCDIR = $(MFEM_ROOT)/include
MFEM_LIBDIR = $(MFEM_ROOT)/lib

MPI_COMPILE_FLAGS = $(shell ${MPI_CXX} --showme:compile)
MPI_LINK_FLAGS = $(shell ${MPI_CXX} --showme:link)

# Setup the include and link flags
IFLAGS = -I$(MFEM_INCDIR) \
		 -I$(HYPRE_INCDIR) \
		 -I$(METIS_INCDIR) \
		 ${MPI_COMPILE_FLAGS}

LFLAGS = -L$(MFEM_LIBDIR) -lmfem \
		 -L$(HYPRE_LIBDIR) -lHYPRE \
		 -L$(METIS_LIBDIR) -lmetis \
		 ${MPI_LINK_FLAGS}

# openmp
ifeq ($(CARDIAC_FIBERS_HAS_OPENMP), YES)
CFLAGS += -fopenmp
LFLAGS += -lgomp
endif

# HIP
ifeq ($(CARDIAC_FIBERS_HAS_HIP), YES)
HIP_CXX ?= hipcc
CXX=$(HIP_CXX)
CFLAGS += $(shell hipconfig -C)
LFLAGS += -lhipsparse

ifeq ($(HIP_TRACE), YES)
CFLAGS += -DHIP_TRACE
IFLAGS += -I$(shell hipconfig --path)/../roctracer/include
LFLAGS += -lroctx64
endif
endif

# CUDA
ifeq ($(CARDIAC_FIBERS_HAS_CUDA), YES)
CUDA_CXX ?= nvcc
CXX=$(CUDA_CXX)
# Wrap the old CFLAGS into --compiler-options, and set the needed nvcc specific flags
CFLAGS := -ccbin=mpiCC -x=cu --expt-extended-lambda --compiler-options="$(CFLAGS)"
LFLAGS += -lcusparse -lrt
endif

SRC = cardiac-fibers.cpp util.cpp calculus.cpp fem.cpp
OBJ=$(SRC:.cpp=.o)

all: cardiac-fibers

cardiac-fibers: check-env $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LFLAGS)

tests: tests.o calculus.o util.o
	$(CXX) -o $@ $^ $(LFLAGS)

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(IFLAGS) -o $@ $^

# Check that all the needed environment variables are set
check-env:
ifndef MFEM_ROOT
	$(error MFEM_ROOT is not set!)
endif
ifeq ($(DEBUG), YES)
ifndef MFEM_DBG_ROOT
	$(warning Warning: MFEM_DBG_ROOT is not set, using MFEM_ROOT)
endif
endif
ifndef HYPRE_INCDIR
	$(error HYPRE_INCDIR is not set!)
endif
ifndef HYPRE_LIBDIR
	$(error HYPRE_LIBDIR is not set!)
endif
ifndef METIS_INCDIR
	$(error METIS_INCDIR is not set!)
endif
ifndef METIS_LIBDIR
	$(error METIS_LIBDIR is not set!)
endif

clean:
	$(RM) cardiac-fibers tests $(OBJ) tests.o

.PHONY: check-env clean
