#!/usr/bin/env bash
#SBATCH --job-name=build-mfem-dgx2q-dbg
#SBATCH --partition=dgx2q
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=03:00:00

# Script for building a debug build of mfem-4.5 on the defq partition of eX3.
# Assumes that it is run inside the root of the mfem source directory.
#
# $ tar -xzvf mfem-4.5.tgz
# $ cd mfem-4.5
# $ # Copy this script in here
# $ sbatch build_mfem_defq_dbg.sbatch

MFEM_BUILD_DIR=./build-dgx2q-dbg
MFEM_INSTALL_DIR=/global/D1/homes/iverh/packages/dgx2q/mfem-dbg-4.5

set -e
module purge
module use /global/D1/homes/james/ex3modules/dgx2q/20221215/modulefiles

module load cuda11.0/toolkit/11.0.3
module load hypre-32-2.25.0
module load openmpi-4.1.4
module load metis-32-5.1.0

make BUILD_DIR=${MFEM_BUILD_DIR} config \
    MFEM_DEBUG=YES \
    MFEM_USE_MPI=YES \
    MPICXX=mpic++ \
    MFEM_USE_CUDA=YES \
    MFEM_USE_METIS=YES \
    MFEM_USE_METIS_5=YES \
    HYPRE_LIB="-L${HYPRE_LIBDIR} -lHYPRE" \
    HYPRE_OPT="-I${HYPRE_INCDIR}" \
    METIS_LIB="-L${METIS_LIBDIR} -lmetis" \
    METIS_OPT="-I${METIS_INCDIR}"
    CUDA_ARCH=sm_70

make BUILD_DIR=${MFEM_BUILD_DIR} -j 4

make BUILD_DIR=${MFEM_BUILD_DIR} check

make BUILD_DIR=${MFEM_BUILD_DIR} install PREFIX=${MFEM_INSTALL_DIR}

