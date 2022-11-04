#!/usr/bin/env bash
#SBATCH --job-name=build-mfem-defq
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=03:00:00

set -e
module purge
module use /global/D1/homes/james/ex3modules/defq/1.0.0/modulefiles

module load hypre-32-2.25.0
module load openmpi-4.1.4
module load metis-32-5.1.0

make BUILD_DIR=./build-defq config\
    MFEM_USE_MPI=YES \
    MPICXX=mpiCC \
    MFEM_USE_METIS=YES \
    MFEM_USE_METIS_5=YES \
    HYPRE_LIB="-L${HYPRE_LIBDIR} -lHYPRE" \
    HYPRE_OPT="-I${HYPRE_INCDIR}" \
    METIS_LIB="-L${METIS_LIBDIR} -lmetis" \
    METIS_OPT="-I${METIS_INCDIR}"

make BUILD_DIR=./build-defq -j 4

make BUILD_DIR=./build-defq check

make BUILD_DIR=./build-defq install PREFIX=/global/D1/homes/iverh/packages/defq/mfem-4.5

