#!/usr/bin/env bash
#SBATCH --job-name=build-mfem-mi100q
#SBATCH --partition=mi100q
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=03:00:00
#SBATCH --nodelist=n004

set -e
module purge
module use /global/D1/homes/james/ex3modules/mi100q/1.0.0/modulefiles

module load amd/rocm/5.1.3
module load hypre-32-2.24.0
module load openmpi-4.1.4
module load metis-32-5.1.0

make BUILD_DIR=./build-mi100q config\
    MFEM_USE_MPI=YES \
    MPICXX=mpiCC \
    MFEM_USE_HIP=YES \
    MFEM_USE_METIS=YES \
    MFEM_USE_METIS_5=YES \
    HYPRE_LIB="-L${HYPRE_LIBDIR} -lHYPRE" \
    HYPRE_OPT="-I${HYPRE_INCDIR}" \
    METIS_LIB="-L${METIS_LIBDIR} -lmetis" \
    METIS_OPT="-I${METIS_INCDIR}" \
    HIP_ARCH=gfx908

make BUILD_DIR=./build-mi100q -j 4

make BUILD_DIR=./build-mi100q check

make BUILD_DIR=./build-mi100q install PREFIX=/global/D1/homes/iverh/packages/mi100q/mfem-4.5

