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
CFLAGS += -O3 -Wall
else
CFLAGS += -O0 -DDEBUG -g3 -Wall
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

# openmp
ifeq ($(LDRB_HAS_OPENMP), YES)
CFLAGS += -fopenmp
LFLAGS += -lgomp
endif

# HIP
ifeq ($(LDRB_HAS_HIP), YES)
CC=hipcc
CFLAGS += $(shell hipconfig -C)
LFLAGS += -lhipsparse

ifeq ($(HIP_TRACE), YES)
CFLAGS += -DHIP_TRACE
IFLAGS += -I$(shell hipconfig --path)/../roctracer/include
LFLAGS += -lroctx64
endif
endif

# CUDA
ifeq ($(LDRB_HAS_CUDA), YES)
CC=nvcc
# Wrap the old CFLAGS into --compiler-options, and set the needed nvcc specific flags
CFLAGS := -ccbin=mpiCC -x=cu --expt-extended-lambda --compiler-options="$(CFLAGS)"
LFLAGS += -lcusparse -lrt
endif

SRC = ldrb-gpu.cpp util.cpp calculus.cpp
OBJ=$(SRC:.cpp=.o)

all: ldrb-gpu

ldrb-gpu: check-env $(OBJ)
	$(CC) -o $@ $(OBJ) $(LFLAGS)

tests: tests.o calculus.o util.o
	$(CC) -o $@ $^ $(LFLAGS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(IFLAGS) -o $@ $^

# Check that all the needed environment variables are set
check-env:
ifndef MFEM_ROOT
	$(error MPI_ROOT is not set!)
endif
ifeq ($(DEBUG), YES)
ifndef MFEM_DBG_ROOT
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
	$(RM) ldrb-gpu tests ldrb-gpu.o calculus.o util.o tests.o

.PHONY: check-env clean
