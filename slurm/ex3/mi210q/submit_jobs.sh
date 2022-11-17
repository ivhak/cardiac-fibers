#!/usr/bin/env sh

mkdir -p jobs 

sbatch --job-name cardiac-fibers-mi210q-N1-n1  --nodes 1 --gres=rsmi:1 --ntasks=1  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n2  --nodes 1 --gres=rsmi:1 --ntasks=2  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n4  --nodes 1 --gres=rsmi:1 --ntasks=4  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n8  --nodes 1 --gres=rsmi:1 --ntasks=8  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n16 --nodes 1 --gres=rsmi:1 --ntasks=16 run_mi210q.sbatch

sbatch --job-name cardiac-fibers-mi210q-N2-n2  --nodes 2 --gres=rsmi:2 --ntasks=2  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N2-n4  --nodes 2 --gres=rsmi:2 --ntasks=4  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N2-n8  --nodes 2 --gres=rsmi:2 --ntasks=8  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N2-n16 --nodes 2 --gres=rsmi:2 --ntasks=16 run_mi210q.sbatch
