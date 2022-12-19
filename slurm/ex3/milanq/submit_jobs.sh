#!/usr/bin/env sh

mkdir -p jobs

COMMIT=HEAD

sbatch --job-name cardiac-fibers-milanq-N4-n512 --nodes 4 --ntasks=512 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N2-n256 --nodes 2 --ntasks=256 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n128 --nodes 1 --ntasks=128 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n64  --nodes 1 --ntasks=64  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n32  --nodes 1 --ntasks=32  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n16  --nodes 1 --ntasks=16  --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n8   --nodes 1 --ntasks=8   --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n4   --nodes 1 --ntasks=4   --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n2   --nodes 1 --ntasks=2   --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch
sbatch --job-name cardiac-fibers-milanq-N1-n1   --nodes 1 --ntasks=1   --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_milanq.sbatch

