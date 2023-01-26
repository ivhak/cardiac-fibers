#!/usr/bin/env bash
#SBATCH --job-name=cardiac-fibers
#SBATCH --partition=hgx2q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4-00:00:00
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --output=jobs/%x-%j.out
#SBATCH --error=jobs/%x-%j.err

LDRB_DIR="${HOME}/master/cardiac-fibers"

DATA_DIR="/global/D1/homes/iverh/data/meshes/martinez-navarro-etal/mesh/mfem"

JOB_DIR="./jobs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}"

BUILD_DIR=

function run {
    srun --cpu-bind=rank ./gpu-wrapper.sh ${BUILD_DIR}/cardiac-fibers/cardiac-fibers \
    --par-mesh \
    --mesh ${DATA_DIR}/"$1"/${SLURM_NTASKS}/"$1".mesh \
    --out ${JOB_DIR}/"$1" \
    --time-to-file \
    --verbose 4 \
    --apex '345 1232 170' \
    --no-save-mfem ${@:2}
}

function main {
    set -o errexit

    BUILD_DIR=/global/D1/homes/iverh/work/cardiac-fibers-${SLURM_JOB_ID}
    mkdir -p ${BUILD_DIR}
    pushd ${BUILD_DIR}
    git clone git@github.com:ivhak/cardiac-fibers.git
    pushd cardiac-fibers
    git checkout ${CARDIAC_FIBERS_COMMIT}
    . envsetup-ex3.sh hgx2q
    make TIMING_BARRIERS=YES clean cardiac-fibers
    popd
    popd

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
    nvidia-smi
    git -C "${BUILD_DIR}/cardiac-fibers" describe --always --abbrev=40 --dirty
    git -C "${BUILD_DIR}/cardiac-fibers" log --oneline --decorate=short | head -n 10

    make -C "${BUILD_DIR}/cardiac-fibers" clean cardiac-fibers

    mkdir -p ${JOB_DIR}

    run heart01 --no-save-paraview
    run heart02 --save-paraview
    run heart03 --no-save-paraview
    run heart04 --no-save-paraview
    run heart05 --no-save-paraview
    run heart06 --no-save-paraview
    run heart07 --no-save-paraview
    run heart08 --no-save-paraview
    run heart09 --no-save-paraview
}

main
