#!usr/bin/env sh

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

# Script for setting the environment when running on MacOS.
# Tested on MacOS Catalina 10.15.7, with mfem-4.5 built from source and
# installed into ~/.local/packages/mfem-4.5. MPI, hypre and metis all installed
# with homebrew:
#
#     $ brew install mpi
#     $ brew install hypre
#     $ brew install metis
#
# The installation of mfem used is built using the script in
# mfem_build_scripts/build_mfem_darwin.sh

PKG_ROOT="${HOME}/.local/packages"

CARDIAC_FIBERS_HAS_HIP=NO
CARDIAC_FIBERS_HAS_CUDA=NO
CARDIAC_FIBERS_HAS_OPENMP=NO

MFEM_ROOT=${PKG_ROOT}/mfem-4.5
HYPRE_ROOT=/usr/local
METIS_ROOT=/usr/local
MPI_ROOT=/usr/local

export MFEM_ROOT

export HYPRE_INCDIR=${HYPRE_ROOT}/include
export HYPRE_LIBDIR=${HYPRE_ROOT}/lib

export METIS_INCDIR=${METIS_ROOT}/include
export METIS_LIBDIR=${METIS_ROOT}/lib

export MPI_INCDIR=${MPI_ROOT}/include
export MPI_LIBDIR=${MPI_ROOT}/lib

export CARDIAC_FIBERS_HAS_HIP=NO
export CARDIAC_FIBERS_HAS_CUDA=NO
export CARDIAC_FIBERS_HAS_OPENMP=NO
