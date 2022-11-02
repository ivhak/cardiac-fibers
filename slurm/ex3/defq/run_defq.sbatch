#!/usr/bin/env bash
#SBATCH --job-name=ldrb-gpu
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --output=jobs/%x-defq-%j-out.txt
#SBATCH --error=jobs/%x-defq-%j-err.txt

LDRB_DIR="${HOME}/master/ldrb-gpu"

DATA_DIR="/global/D1/homes/iverh/data/meshes/martinez-navarro-etal/mesh/gmsh"

JOB_DIR="${SLURM_SUBMIT_DIR}/jobs/${SLURM_JOB_NAME}-defq-${SLURM_JOB_ID}"
OUT_DIR="${JOB_DIR}/out"


function run {
    mpirun -np ${SLURM_NTASKS} ${LDRB_DIR}/ldrb-gpu \
        --mesh ${DATA_DIR}/"$1".msh \
        --out "${OUT_DIR}/$1" \
        --verbose 3 \
        --apex '345 1232 170' \
        --paraview \
        --alpha-endo 60.0 \
        --alpha-epi -60.0 \
        --beta-endo -60.0 \
        --beta-epi   60.0

}

function main {
    set -o errexit


    pushd "${LDRB_DIR}"
    . envsetup.sh defq

    set -x
    make clean ldrb-gpu
    popd

    mkdir -p ${JOB_DIR}
    mkdir -p ${OUT_DIR}

    run heart01
    run heart02
    run heart03
    run heart04
}

main
