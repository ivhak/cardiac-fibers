#!/usr/bin/env bash

mkdir -p jobs

COMMIT=45c15c9

sbatch --job-name cardiac-fibers-hgx2q-N1-n8-g8 -n 8 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_hgx2q.sbatch
sbatch --job-name cardiac-fibers-hgx2q-N1-n4-g4 -n 4 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_hgx2q.sbatch
sbatch --job-name cardiac-fibers-hgx2q-N1-n2-g2 -n 2 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_hgx2q.sbatch
sbatch --job-name cardiac-fibers-hgx2q-N1-n1-g1 -n 1 --export=ALL,CARDIAC_FIBERS_COMMIT=${COMMIT} run_hgx2q.sbatch
