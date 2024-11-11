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
# python3 train_clip.py --data_dir data/human_cite/mudata_human_cite_full.h5mu --logit_scale 1
# python3 train_clip.py --data_dir data/human_multiome/mudata_human_multiome_full.h5mu --logit_scale 1

# python3 train_clip.py --data_dir simulated_train --logit_scale 1
# python3 train_clip.py --data_dir simulated_test --logit_scale 1 --checkpoint results/simulated_train/1.0_False_30000_0.00015_/lightning_logs/checkpoints/epoch=85-step=30000.ckpt

# python3 train_clip.py --data_dir data/human_cite/mudata_human_cite_train.h5mu --logit_scale 1
# python3 train_clip.py --data_dir data/human_cite/mudata_human_cite_test.h5mu --logit_scale 1 --checkpoint results/data/human_cite/mudata_human_cite_train.h5mu/1.0_False_30000_0.00015_/lightning_logs/checkpoints/epoch=115-step=30000.ckpt

# python3 train_clip.py --data_dir data/human_multiome/mudata_human_multiome_train.h5mu --logit_scale 1
python3 train_clip.py --data_dir data/human_multiome/mudata_human_multiome_test.h5mu --logit_scale 1 --checkpoint results/data/human_multiome/mudata_human_multiome_train.h5mu/1.0_False_30000_0.00015_/lightning_logs/checkpoints/epoch=180-step=30000.ckpt

