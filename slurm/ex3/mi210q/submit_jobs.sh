#!/usr/bin/env sh

mkdir -p jobs

COMMIT=HEAD

sbatch --job-name cardiac-fibers-mi210q-N1-n8-g1 --nodes 1 --gres=rsmi:1 --ntasks=8 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n4-g1 --nodes 1 --gres=rsmi:1 --ntasks=4 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n2-g1 --nodes 1 --gres=rsmi:1 --ntasks=2 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n1-g1 --nodes 1 --gres=rsmi:1 --ntasks=1 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_mi210q.sbatch
