#!/bin/bash

#SBATCH -o logs/log-%j-multimodal.out
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=10:00:00
#SBATCH --gres=gpumem:21G
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=END,FAIL

mkdir -p logs

#module load gcc/8.2.0 python_gpu/3.9.9
# module load stack/.2024-04-silent stack/2024-04
# module load gcc/8.5.0
# module --ignore_cache load python/3.9.18
nvidia-smi

conda init
conda activate myenv

# python3 train_clip.py --data_dir data/Fetal/fetal.h5mu --logit_scale 1
# python3 train_clip.py --data_dir data/simulated/mudata_simulated_full.h5mu --logit_scale 1
python3 train_clip.py --data_dir data/simulated/mudata_simulated_train.h5mu --logit_scale 1
