#!/usr/bin/env bash
#SBATCH --job-name=build-mfem-dgx2q-dbg
#SBATCH --partition=dgx2q
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=03:00:00

set -e
module purge
module use /global/D1/homes/james/ex3modules/dgx2q/1.0.0/modulefiles

module load cuda10.1/toolkit/10.1.243
module load hypre-32-2.24.0
module load openmpi-4.1.4
module load metis-32-5.1.0

make BUILD_DIR=./build-dgx2q-dbg config \
    MFEM_DEBUG=YES \
    MFEM_USE_MPI=YES \
    MPICXX=mpiCC \
    MFEM_USE_CUDA=YES \
    MFEM_USE_METIS=YES \
    MFEM_USE_METIS_5=YES \
    HYPRE_LIB="-L${HYPRE_LIBDIR} -lHYPRE" \
    HYPRE_OPT="-I${HYPRE_INCDIR}" \
    METIS_LIB="-L${METIS_LIBDIR} -lmetis" \
    METIS_OPT="-I${METIS_INCDIR}"

make BUILD_DIR=./build-dgx2q-dbg -j 4

make BUILD_DIR=./build-dgx2q-dbg check

make BUILD_DIR=./build-dgx2q-dbg install PREFIX=/global/D1/homes/iverh/packages/dgx2q/mfem-dbg-4.5

