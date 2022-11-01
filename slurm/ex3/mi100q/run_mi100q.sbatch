#!/usr/bin/env bash
#SBATCH --job-name=ldrb-gpu
#SBATCH --partition=mi100q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --nodelist=n004
#SBATCH --output=jobs/%x-mi100q-%j-out.txt
#SBATCH --error=jobs/%x-mi100q-%j-err.txt

LDRB_DIR="${HOME}/master/ldrb-gpu"

DATA_DIR="/global/D1/homes/iverh/data/meshes/martinez-navarro-etal/mesh/gmsh"

JOB_DIR="${SLURM_SUBMIT_DIR}/jobs/${SLURM_JOB_NAME}-mi100q-${SLURM_JOB_ID}"
OUT_DIR="${JOB_DIR}/out"


function run {
    mpirun -np ${SLURM_NTASKS} ${LDRB_DIR}/ldrb-gpup \
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
    . envsetup.sh mi100q

    set -x
    make GPU_CALCULUS=YES clean ldrb-gpup
    popd

    mkdir -p ${JOB_DIR}
    mkdir -p ${OUT_DIR}

    run heart01
    run heart02
    run heart03
    run heart04
}

main
