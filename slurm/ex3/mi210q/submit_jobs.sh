#!/usr/bin/env bash

mkdir -p jobs

COMMIT=HEAD

sbatch --job-name cardiac-fibers-mi210q-N1-n8-g1 -N 1 --gres=rsmi:1 -n 8 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n4-g1 -N 1 --gres=rsmi:1 -n 4 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n2-g1 -N 1 --gres=rsmi:1 -n 2 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n1-g1 -N 1 --gres=rsmi:1 -n 1 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch

sbatch --job-name cardiac-fibers-mi210q-N1-n16-g2 -N 1 --gres=rsmi:2 -n 16 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n8-g2  -N 1 --gres=rsmi:2 -n 8  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n4-g2  -N 1 --gres=rsmi:2 -n 4  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n2-g2  -N 1 --gres=rsmi:2 -n 2  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch

sbatch --job-name cardiac-fibers-mi210q-N2-n32-g4 -N 2 --gres=rsmi:2 -n 32 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N2-n16-g4 -N 2 --gres=rsmi:2 -n 16 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N2-n4-g4  -N 2 --gres=rsmi:2 -n 4  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N2-n8-g4  -N 2 --gres=rsmi:2 -n 8  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch

