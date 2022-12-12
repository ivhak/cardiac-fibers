#!/usr/bin/env sh

mkdir -p jobs 

sbatch --job-name cardiac-fibers-mi210q-N1-n128-g2 --nodes 1 --gres=rsmi:2 --ntasks=128 run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n64-g2  --nodes 1 --gres=rsmi:2 --ntasks=64  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n32-g2  --nodes 1 --gres=rsmi:2 --ntasks=32  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n16-g2  --nodes 1 --gres=rsmi:2 --ntasks=16  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n8-g2   --nodes 1 --gres=rsmi:2 --ntasks=8   run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n4-g2   --nodes 1 --gres=rsmi:2 --ntasks=4   run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n2-g2   --nodes 1 --gres=rsmi:2 --ntasks=2   run_mi210q.sbatch
sbatch --job-name cardiac-fibers-mi210q-N1-n1-g1   --nodes 1 --gres=rsmi:1 --ntasks=1   run_mi210q.sbatch

sbatch --job-name cardiac-fibers-N2-n128-g4 --nodes 2 --gres=rsmi:2 --ntasks=128 run_mi210q.sbatch
sbatch --job-name cardiac-fibers-N2-n64-g4  --nodes 2 --gres=rsmi:2 --ntasks=64  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-N2-n32-g4  --nodes 2 --gres=rsmi:2 --ntasks=32  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-N2-n16-g4  --nodes 2 --gres=rsmi:2 --ntasks=16  run_mi210q.sbatch
sbatch --job-name cardiac-fibers-N2-n8-g4   --nodes 2 --gres=rsmi:2 --ntasks=8   run_mi210q.sbatch
sbatch --job-name cardiac-fibers-N2-n4-g4   --nodes 2 --gres=rsmi:2 --ntasks=4   run_mi210q.sbatch
