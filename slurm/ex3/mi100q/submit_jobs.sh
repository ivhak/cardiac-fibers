#!/usr/bin/env sh

mkdir -p jobs

sbatch --job-name cardiac-fibers-mi100q-N1-n1-g1 --gres=rsmi:1 --nodes 1 --ntasks=1  run_mi100q.sbatch
