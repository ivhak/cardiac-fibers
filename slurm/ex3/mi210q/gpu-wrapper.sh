#!/usr/bin/env sh

# Each node on the mi210q has two GPUs. Given a job with n tasks per node,
# assign the first n/2 ranks to GPU 0, and the remaining n/2 ranks to GPU 1.

JOBS_PER_NODE="$(echo "$SLURM_NTASKS/$SLURM_JOB_NUM_NODES" | bc)"
GPUS_PER_NODE=2
GPU_ID="$(echo "$SLURM_LOCALID >= $JOBS_PER_NODE/$GPUS_PER_NODE" | bc)"

echo "Job $SLURM_LOCALID on node $SLURM_NODEID: GPU $GPU_ID"

export ROCR_VISIBLE_DEVICES="$GPU_ID"

"$@"
