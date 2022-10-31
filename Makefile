# Copyright (C) 2022 Iver Håkonsen
#
# ldrb-gpu is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ldrb-gpu is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ldrb-gpu.  If not, see <https://www.gnu.org/licenses/>.
#
# Authors: Iver Håkonsen <hakonseniver@yahoo.no


CFLAGS =
ifneq ($(DEBUG), YES)
CFLAGS += -O3
else
CFLAGS += -O0 -DDEBUG -g3 -Wall
endif

ifeq ($(GPU_CALCULUS), YES)
CFLAGS += -DGPU_CALCULUS
endif

CC=mpiCC

# MFEM
ifeq ($(DEBUG), YES)
MFEM_ROOT=$(MFEM_DBG_ROOT)
endif

MFEM_INCDIR = $(MFEM_ROOT)/include
MFEM_LIBDIR = $(MFEM_ROOT)/lib

# Setup the include and link flags
IFLAGS = -I$(MFEM_INCDIR) \
		 -I$(MPI_INCDIR) \
		 -I$(HYPRE_INCDIR) \
		 -I$(METIS_INCDIR)

LFLAGS = -L$(MFEM_LIBDIR) -lmfem \
		 -L$(MPI_LIBDIR) -lmpi -lmpi_cxx \
		 -L$(HYPRE_LIBDIR) -lHYPRE \
		 -L$(METIS_LIBDIR) -lmetis

# HIP
ifeq ($(LDRB_HAS_HIP), YES)
CC=hipcc
CFLAGS += $(shell hipconfig -C)
LFLAGS += -lhipsparse

ifeq ($(HIP_TRACE), YES)
CFLAGS += -DHIP_TRACE
IFLAGS += -I$(ROCm_ROOT)/roctracer/include
LFLAGS += -lroctx64
endif
endif

# CUDA
ifeq ($(LDRB_HAS_CUDA), YES)
CC=nvcc
CFLAGS += -ccbin=mpiCC -x=cu --expt-extended-lambda
LFLAGS += -lcusparse -lrt
endif

SRC = ldrb-gpup.cpp util.cpp
ifeq ($(GPU_CALCULUS), YES)
SRC += calculus_gpu.cpp
else
SRC += calculus.cpp
endif

OBJ=$(SRC:.cpp=.o)

all: ldrb-gpup


ldrb-gpup: .check-env
ldrb-gpup: $(OBJ)
	$(CC) -o $@ $(OBJ) $(LFLAGS)

tests: tests.o calculus.o util.o
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $^ $(LFLAGS)

tests_gpu: tests_gpu.o calculus_gpu.o util.o
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(IFLAGS) -o $@ $^

# Check that all the needed environment variables are set
.check-env:
ifeq ($(LDRB_HAS_CUDA), YES)
ifndef CUDA_INCDIR
	$(error CUDA_INCDIR is not set!)
endif
ifndef CUDA_LIBDIR
	$(error CUDA_LIBDIR is not set!)
endif
endif
ifndef MPI_INCDIR
	$(error MPI_INCDIR is not set!)
endif
ifndef MPI_LIBDIR
	$(error MPI_LIBDIR is not set!)
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
	$(RM) ldrb-gpup tests test_gpu calculus.o calculus_gpu.o ldrb-gpu.o ldrb-gpup.o util.o tests.o

.PHONY: .check-env clean
