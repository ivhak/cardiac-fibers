#!usr/bin/env sh

# Copyright (C) 2022 Iver Håkonsen
#
# cardiac-fibers is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with
# cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
#
# Authors: Iver Håkonsen <hakonseniver@yahoo.no

# Script for setting the environment when running on eX³.


# The builds of mfem are in the following directory
PKG_ROOT="/global/D1/homes/iverh/packages"

CARDIAC_FIBERS_HAS_HIP=NO
CARDIAC_FIBERS_HAS_CUDA=NO
CARDIAC_FIBERS_HAS_OPENMP=NO

module load slurm/20.02.7

if [ "$1" = "mi100q" ] || [ -z "$1" ]; then
    # HIP build, default
    echo "Setting up environment for mi100q (HIP)"

    module use /global/D1/homes/james/ex3modules/mi100q/1.0.0/modulefiles

    module load amd/rocm/5.1.3
    module load hypre-32-2.24.0
    module load openmpi-4.1.4
    module load metis-32-5.1.0

    export MFEM_ROOT="${PKG_ROOT}/mi100q/mfem-4.5"
    export MFEM_DBG_ROOT="${PKG_ROOT}/mi100q/mfem-dbg-4.5"

    CARDIAC_FIBERS_HAS_HIP=YES

elif [ "$1" = "mi210q" ]; then
    # HIP build
    echo "Setting up environment for mi210q (HIP)"

    module use /global/D1/homes/james/ex3modules/mi210q/20221215/modulefiles

    module load amd/rocm/5.1.3
    module load hypre-32-2.25.0-gfx90a
    module load openmpi-4.1.4
    module load metis-32-5.1.0

    export MFEM_ROOT="${PKG_ROOT}/mi210q/mfem-4.5"
    export MFEM_DBG_ROOT="${PKG_ROOT}/mi210q/mfem-dbg-4.5"

    CARDIAC_FIBERS_HAS_HIP=YES

elif [ "$1" = "dgx2q" ]; then
    # CUDA build
    echo "Setting up environment for dgx2q (CUDA)"

    module use /global/D1/homes/james/ex3modules/dgx2q/20221210/modulefiles

    module load cuda11.0/toolkit/11.0.3
    module load hypre-32-2.25.0
    module load openmpi-4.1.4
    module load metis-32-5.1.0

    export MFEM_ROOT="${PKG_ROOT}/dgx2q/mfem-4.5"
    export MFEM_DBG_ROOT="${PKG_ROOT}/dgx2q/mfem-dbg-4.5"

    CARDIAC_FIBERS_HAS_CUDA=YES

elif [ "$1" = "hgx2q" ]; then
    # CUDA build
    echo "Setting up environment for hgx2q (CUDA)"

    module use /global/D1/homes/james/ex3modules/hgx2q/2022-08-17/modulefiles

    module load cuda11.0/toolkit/11.0.3
    module load hypre-32-2.25.0
    module load openmpi-4.1.4
    module load metis-32-5.1.0

    export MFEM_ROOT="${PKG_ROOT}/hgx2q/mfem-4.5"
    export MFEM_DBG_ROOT="${PKG_ROOT}/hgx2q/mfem-dbg-4.5"

    CARDIAC_FIBERS_HAS_CUDA=YES

elif [ "$1" = "defq" ]; then
    # CPU build
    echo "Setting up environment for defq (cpu)"

    module use /global/D1/homes/james/ex3modules/defq/1.0.0/modulefiles

    module load hypre-32-2.25.0
    module load openmpi-4.1.4
    module load metis-32-5.1.0

    export MFEM_ROOT="${PKG_ROOT}/defq/mfem-4.5"
    export MFEM_DBG_ROOT="${PKG_ROOT}/defq/mfem-dbg-4.5"

    CARDIAC_FIBERS_HAS_OPENMP=YES
fi

export CARDIAC_FIBERS_HAS_HIP
export CARDIAC_FIBERS_HAS_CUDA
export CARDIAC_FIBERS_HAS_OPENMP
