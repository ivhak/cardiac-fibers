#!/usr/bin/env bash
#SBATCH --job-name=ldrb-gpu
#SBATCH --partition=dgx2q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=jobs/%x-dgx2q-%j-out.txt
#SBATCH --error=jobs/%x-dgx2q-%j-err.txt

LDRB_DIR="${HOME}/master/ldrb-gpu"

DATA_DIR="/global/D1/homes/iverh/data/meshes/martinez-navarro-etal/mesh/gmsh"

JOB_DIR="./jobs/${SLURM_JOB_NAME}-dgx2q-${SLURM_JOB_ID}"

BUILD_DIR=


function run {
    mpirun -np ${SLURM_NTASKS} ${BUILD_DIR}/ldrb-gpu \
        --mesh ${DATA_DIR}/"$1".msh \
        --out ${JOB_DIR}/"$1" \
        --time-to-file \
        --verbose 3 \
        --apex '345 1232 170' \
        --no-save-paraview \
        --no-save-mfem
}

function main {
    set -o errexit

    TMPDIR=$(mktemp -d)


    BUILD_DIR=${TMPDIR}/ldrb-gpu
    cp -r ${LDRB_DIR} ${BUILD_DIR}

    . ${BUILD_DIR}/envsetup-ex3.sh dgx2q


    # Print some diagnostics
    echo "SLURM_JOB_NAME=${SLURM_JOB_NAME}"
    echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
    echo "SLURM_SUBMIT_HOST=${SLURM_SUBMIT_HOST}"
    echo "SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}"
    echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST}"
    echo "SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}"
    echo "TMPDIR=${TMPDIR}"
    echo "BUILD_DIR=${BUILD_DIR}"
    set -x
    uname -a
    lscpu
    numactl -H
    nvidia-smi
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
}

main
