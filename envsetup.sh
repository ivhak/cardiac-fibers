#!usr/bin/env sh

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

# Script for setting the environment when running on the mi100q of eX³.

PKG_ROOT="/global/D1/homes/iverh/packages"

LDRB_HAS_HIP=NO

module purge

if [ "$1" = "mi100q" ] || [ -z "$1" ]; then
    # HIP build, default

    module use /global/D1/homes/james/ex3modules/mi100q/1.0.0/modulefiles

    module load amd/rocm/5.1.3
    module load hypre-32-2.24.0
    module load openmpi-4.1.4
    module load metis-5.1.0

    # The openmpi-4.1.4 module does not properly export MPI_INCDIR and MPI_LIBDIR
    export MPI_INCDIR=${MPI_HOME}/include
    export MPI_LIBDIR=${MPI_HOME}/lib

    export MFEM_ROOT="${PKG_ROOT}/mi100q/mfem-4.5"
    export MFEM_DBG_ROOT="${PKG_ROOT}/packages/mi100q/mfem-dbg-4.5"
    export LDRB_HAS_HIP=YES

elif [ "$1" = "defq" ]; then
    # CPU build

    module use /global/D1/homes/james/ex3modules/defq/1.0.0/modulefiles

    module load hypre-32-2.25.0
    module load openmpi-4.1.4
    module load metis-32-5.1.0

    # The openmpi-4.1.4 module does not properly export MPI_INCDIR and MPI_LIBDIR
    export MPI_INCDIR=${MPI_HOME}/include
    export MPI_LIBDIR=${MPI_HOME}/lib

    export MFEM_ROOT="${PKG_ROOT}/defq/mfem-4.5"
    export MFEM_DBG_ROOT="${PKG_ROOT}/defq/mfem-dbg-4.5"
fi

export LDRB_HAS_HIP

