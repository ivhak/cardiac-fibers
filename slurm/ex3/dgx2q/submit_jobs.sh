#!/usr/bin/env sh

mkdir -p jobs

sbatch --job-name cardiac-fibers-dgx2q-N1-n1  --nodes 1 --ntasks=1  --gres=gpu:1 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n2  --nodes 1 --ntasks=2  --gres=gpu:1 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n4  --nodes 1 --ntasks=4  --gres=gpu:1 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n8  --nodes 1 --ntasks=8  --gres=gpu:1 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n16 --nodes 1 --ntasks=16 --gres=gpu:1 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n32 --nodes 1 --ntasks=32 --gres=gpu:1 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n64 --nodes 1 --ntasks=64 --gres=gpu:1 run_dgx2q.sbatch

sbatch --job-name cardiac-fibers-dgx2q-N1-n8-gpu2  --nodes 1 --ntasks=8  --gres=gpu:2 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n16-gpu2 --nodes 1 --ntasks=16 --gres=gpu:2 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n32-gpu2 --nodes 1 --ntasks=32 --gres=gpu:2 run_dgx2q.sbatch
sbatch --job-name cardiac-fibers-dgx2q-N1-n64-gpu2 --nodes 1 --ntasks=64 --gres=gpu:2 run_dgx2q.sbatch

