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
# along with hyprep.  If not, see <https://www.gnu.org/licenses/>.
#
# Authors: Iver Håkonsen <hakonseniver@yahoo.no

MFEM_ROOT=${HOME}/packages/mfem-4.4
MFEM_INCDIR=$(MFEM_ROOT)/include
MFEM_LIBDIR=$(MFEM_ROOT)/lib

CFLAGS =
#CFLAGS = -Wall -Wextra -march=native -O3
CFLAGS += $(shell hipconfig -C)

IFLAGS=
LFLAGS=

# HIP
CFLAGS += --gcc-toolchain=${GCC_ROOT}
LFLAGS += -lhipsparse

# MFEM
IFLAGS += -I$(MFEM_INCDIR)
LFLAGS += -L$(MFEM_LIBDIR) -lmfem

# MPI
IFLAGS += -I${MPI_HOME}/include
LFLAGS += -L${MPI_HOME}/lib -lmpi

# HYPRE
IFLAGS += -I${HYPRE_INCDIR}
LFLAGS += -L${HYPRE_LIBDIR} -lHYPRE

# Metis
IFLAGS += -I${METIS_INCDIR}
LFLAGS += -L${METIS_LIBDIR} -lmetis

CC=hipcc

ldrb-gpu: ldrb-gpu.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $^  $(LFLAGS)

debug: CFLAGS += -DDEBUG -g -Wall
debug: MFEM_ROOT=${HOME}/packages/mfem-dbg-4.4
debug: ldrb-gpu

.PHONY: clean
clean:
	$(RM) ldrb-gpu
