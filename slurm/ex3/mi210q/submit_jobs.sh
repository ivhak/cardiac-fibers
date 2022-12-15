#!/usr/bin/env sh

mkdir -p jobs

sbatch --job-name cardiac-fibers-mi210q-N1-n1-g1 --nodes 1 --gres=rsmi:1 --ntasks=1 run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n2-g2 --nodes 2 --gres=rsmi:2 --ntasks=2 run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N2-n4-g4 --nodes 2 --gres=rsmi:2 --ntasks=4 run_mi210q.sbatch
