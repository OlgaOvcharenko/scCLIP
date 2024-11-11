#!/bin/bash

#SBATCH -o logs_cpu/log-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=30:00:00
#SBATCH --gpus=1
#SBATCH --gres=gpumem:16G
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-type=END,FAIL

mkdir -p logs

conda activate myenv

# python3 metrics.py --data simulated
python3 metrics.py --data human_cite
# python3 metrics.py --data human_multiome

# python3 qr_scclip.py --data simulated
# python3 qr_scclip.py --data human_cite
# python3 qr_scclip.py --data human_multiome