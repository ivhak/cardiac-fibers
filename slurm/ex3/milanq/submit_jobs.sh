#!/usr/bin/env sh

mkdir -p jobs

sbatch --job-name cardiac-fibers-milanq-N4-n512 --nodes 4 --ntasks=512 run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N2-n256 --nodes 2 --ntasks=256 run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n128 --nodes 1 --ntasks=128 run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n64  --nodes 1 --ntasks=64  run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n32  --nodes 1 --ntasks=32  run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n16  --nodes 1 --ntasks=16  run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n8   --nodes 1 --ntasks=8   run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n4   --nodes 1 --ntasks=4   run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n2   --nodes 1 --ntasks=2   run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n1   --nodes 1 --ntasks=1   run_milanq.sbatch

