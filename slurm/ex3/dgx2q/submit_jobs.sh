#!/usr/bin/env sh

mkdir -p jobs

sbatch --job-name ldrb-gpu-dgx2q-N1-n1  --nodes 1 --ntasks=1  run_dgx2q.sbatch
sbatch --job-name ldrb-gpu-dgx2q-N1-n2  --nodes 1 --ntasks=2  run_dgx2q.sbatch
sbatch --job-name ldrb-gpu-dgx2q-N1-n4  --nodes 1 --ntasks=4  run_dgx2q.sbatch
sbatch --job-name ldrb-gpu-dgx2q-N1-n8  --nodes 1 --ntasks=8  run_dgx2q.sbatch
sbatch --job-name ldrb-gpu-dgx2q-N1-n16 --nodes 1 --ntasks=16 run_dgx2q.sbatch
