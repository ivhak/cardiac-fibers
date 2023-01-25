#!/usr/bin/env sh

export CUDA_VISIBLE_DEVICES="$SLURM_LOCALID"

"$@"
