#!/usr/bin/env sh

mkdir -p jobs

sbatch --job-name ldrb-gpu-mi100q-N1-n1  --fres=rsmi:1 --nodes 1 --ntasks=1  run_mi100q.sbatch
sbatch --job-name ldrb-gpu-mi100q-N1-n2  --fres=rsmi:1 --nodes 1 --ntasks=2  run_mi100q.sbatch
sbatch --job-name ldrb-gpu-mi100q-N1-n4  --fres=rsmi:1 --nodes 1 --ntasks=4  run_mi100q.sbatch
sbatch --job-name ldrb-gpu-mi100q-N1-n8  --fres=rsmi:1 --nodes 1 --ntasks=8  run_mi100q.sbatch
sbatch --job-name ldrb-gpu-mi100q-N1-n16 --fres=rsmi:1 --nodes 1 --ntasks=16 run_mi100q.sbatch

sbatch --job-name ldrb-gpu-mi100q-N2-n2  --gres=rsmi:2 --nodes 1 --ntasks=2  run_mi100q.sbatch
sbatch --job-name ldrb-gpu-mi100q-N2-n4  --gres=rsmi:2 --nodes 1 --ntasks=4  run_mi100q.sbatch
sbatch --job-name ldrb-gpu-mi100q-N2-n8  --gres=rsmi:2 --nodes 1 --ntasks=8  run_mi100q.sbatch
sbatch --job-name ldrb-gpu-mi100q-N2-n16 --gres=rsmi:2 --nodes 1 --ntasks=16 run_mi100q.sbatch
