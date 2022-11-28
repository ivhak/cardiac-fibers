#!/usr/bin/env bash
#SBATCH --job-name=build-mfem-mi100q
#SBATCH --partition=mi100q
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=03:00:00
#SBATCH --nodelist=n004

# Script for building mfem-4.5 on the mi100q partition of eX3.
# Assumes that it is run inside the root of the mfem source directory.
#
# $ tar -xzvf mfem-4.5.tgz
# $ cd mfem-4.5
# $ # Copy this script in here
# $ sbatch build_mfem_mi100q.sbatch

MFEM_BUILD_DIR=./build-mi100q
MFEM_INSTALL_DIR=/global/D1/homes/iverh/packages/mi100q/mfem-4.5

set -e
module purge
module use /global/D1/homes/james/ex3modules/mi100q/1.0.0/modulefiles

module load amd/rocm/5.1.3
module load hypre-32-2.24.0
module load openmpi-4.1.4
module load metis-32-5.1.0

make BUILD_DIR=${MFEM_BUILD_DIR} config\
    MFEM_USE_MPI=YES \
    MPICXX=mpic++ \
    MFEM_USE_HIP=YES \
    MFEM_USE_METIS=YES \
    MFEM_USE_METIS_5=YES \
    HYPRE_LIB="-L${HYPRE_LIBDIR} -lHYPRE" \
    HYPRE_OPT="-I${HYPRE_INCDIR}" \
    METIS_LIB="-L${METIS_LIBDIR} -lmetis" \
    METIS_OPT="-I${METIS_INCDIR}" \
    HIP_ARCH=gfx908

make BUILD_DIR=${MFEM_BUILD_DIR} -j 4

make BUILD_DIR=${MFEM_BUILD_DIR} check

make BUILD_DIR=${MFEM_BUILD_DIR} install PREFIX=${MFEM_INSTALL_DIR}

