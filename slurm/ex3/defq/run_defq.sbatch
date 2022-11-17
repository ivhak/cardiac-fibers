#!/usr/bin/env bash
#SBATCH --job-name=ldrb-gpu
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:00
#SBATCH --output=jobs/%x-%j.err
#SBATCH --error=jobs/%x-%j.err

LDRB_DIR="${HOME}/master/ldrb-gpu"

DATA_DIR="/global/D1/homes/iverh/data/meshes/martinez-navarro-etal/mesh/gmsh"

JOB_DIR="./jobs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}"

BUILD_DIR=


function run {
    mpirun --bind-to core -np ${SLURM_NTASKS} ${BUILD_DIR}/ldrb-gpu \
        --mesh ${DATA_DIR}/"$1".msh \
        --out ${JOB_DIR}/"$1" \
        --device omp \
        --time-to-file \
        --verbose 3 \
        --apex '345 1232 170' \
        --no-save-paraview \
        --no-save-mfem
}

function main {
    set -o errexit

    BUILD_DIR=/global/D1/homes/iverh/ldrb-gpu-${SLURM_JOB_ID}

    cp -r ${LDRB_DIR} ${BUILD_DIR}

    . ${BUILD_DIR}/envsetup-ex3.sh defq


    # Print some diagnostics
    echo "SLURM_JOB_NAME=${SLURM_JOB_NAME}"
    echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
    echo "SLURM_SUBMIT_HOST=${SLURM_SUBMIT_HOST}"
    echo "SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}"
    echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST}"
    echo "SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}"
    echo "BUILD_DIR=${BUILD_DIR}"
    set -x
    uname -a
    lscpu
    numactl -H
    git -C "${BUILD_DIR}" describe --always --abbrev=40 --dirty
    git -C "${BUILD_DIR}" log --oneline --decorate=short | head -n 10

    make -C ${BUILD_DIR} clean ldrb-gpu

    mkdir -p ${JOB_DIR}

    run heart01
    run heart02
    run heart03
    run heart04
    run heart05
    run heart06
    run heart07
    run heart08
    run heart09
    run heart10
    run heart11
    run heart12
    run heart13
}

main
