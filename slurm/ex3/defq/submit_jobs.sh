#!/usr/bin/env sh

mkdir -p jobs

sbatch --job-name cardiac-fibers-N1-n128 --nodes 1 --ntasks=128 run_defq.sbatch
sbatch --job-name cardiac-fibers-N1-n64  --nodes 1 --ntasks=64  run_defq.sbatch
sbatch --job-name cardiac-fibers-N1-n32  --nodes 1 --ntasks=32  run_defq.sbatch
sbatch --job-name cardiac-fibers-N1-n16  --nodes 1 --ntasks=16  run_defq.sbatch
sbatch --job-name cardiac-fibers-N1-n8   --nodes 1 --ntasks=8   run_defq.sbatch
sbatch --job-name cardiac-fibers-N1-n4   --nodes 1 --ntasks=4   run_defq.sbatch
sbatch --job-name cardiac-fibers-N1-n2   --nodes 1 --ntasks=2   run_defq.sbatch
sbatch --job-name cardiac-fibers-N1-n1   --nodes 1 --ntasks=1   run_defq.sbatch

