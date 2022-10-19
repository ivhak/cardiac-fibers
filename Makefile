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

MFEM_SERIAL_ROOT=${HOME}/packages/mfem-serial-4.4
MFEM_PARALLEL_ROOT=${HOME}/packages/mfem-4.4

MFEM_ROOT=
MFEM_INCDIR=$(MFEM_ROOT)/include
MFEM_LIBDIR=$(MFEM_ROOT)/lib

CFLAGS =
ifneq ($(DEBUG), YES)
CFLAGS += -march=native -O3
else
CFLAGS += -O0 -DDEBUG -g -Wall
MFEM_SERIAL_ROOT=${HOME}/packages/mfem-serial-dbg-4.4
MFEM_PARALLEL_ROOT=${HOME}/packages/mfem-dbg-4.4
endif

IFLAGS=
LFLAGS=

# HIP
CFLAGS += $(shell hipconfig -C) --gcc-toolchain=${GCC_ROOT}
HIP_LFLAGS = -lhipsparse

# MFEM
MFEM_IFLAGS = -I$(MFEM_INCDIR)
MFEM_LFLAGS = -L$(MFEM_LIBDIR) -lmfem

# MPI
MPI_IFLAGS = -I${MPI_HOME}/include
MPI_LFLAGS = -L${MPI_HOME}/lib -lmpi -lmpi_cxx

# HYPRE
HYPRE_IFLAGS = -I${HYPRE_INCDIR}
HYPRE_LFLAGS = -L${HYPRE_LIBDIR} -lHYPRE

# Metis
METIS_IFLAGS = -I${METIS_INCDIR}
METIS_LFLAGS = -L${METIS_LIBDIR} -lmetis

CC=hipcc

PARALLEL_IFLAGS = $(HIP_IFLAGS) $(MFEM_IFLAGS) $(MPI_IFLAGS) $(HYPRE_IFLAGS) $(METIS_IFLAGS)
PARALLEL_LFLAGS = $(HIP_LFLAGS) $(MFEM_LFLAGS) $(MPI_LFLAGS) $(HYPRE_LFLAGS) $(METIS_LFLAGS)

SERIAL_IFLAGS = $(HIP_IFLAGS) $(MFEM_IFLAGS)
SERIAL_LFLAGS = $(HIP_LFLAGS) $(MFEM_LFLAGS)


all: ldrb-gpup ldrb-gpu

ldrb-gpup: MFEM_ROOT=$(MFEM_PARALLEL_ROOT)
ldrb-gpup: IFLAGS=$(PARALLEL_IFLAGS)
ldrb-gpup: LFLAGS=$(PARALLEL_LFLAGS)
ldrb-gpup: ldrb-gpup.o calculus.o util.o
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $^ $(LFLAGS)

ldrb-gpu: MFEM_ROOT=$(MFEM_SERIAL_ROOT)
ldrb-gpu: IFLAGS=$(SERIAL_IFLAGS)
ldrb-gpu: LFLAGS=$(SERIAL_LFLAGS)
ldrb-gpu: ldrb-gpu.o calculus.o util.o
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $^ $(LFLAGS)

tests: MFEM_ROOT=$(MFEM_SERIAL_ROOT)
tests: IFLAGS=$(SERIAL_IFLAGS)
tests: LFLAGS=$(SERIAL_LFLAGS)
tests: tests.o calculus.o util.o
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(IFLAGS) -o $@ $^

.PHONY: clean
clean:
	$(RM) ldrb-gpu ldrb-gpup tests calculus.o ldrb-gpu.o util.o tests.o
