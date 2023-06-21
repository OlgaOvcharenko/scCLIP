# Multi-modal Single-cell Integration using Transformers with Contrastive Learning

## Download data 
Download [multiome](https://www.dropbox.com/sh/70caiyjydx3jnq1/AAB51h6PCX9IGgi8jyT5KMhaa?dl=0) data folder under data/  
- dataset:
  - Brain
    - train: AD
    - test: human_brain_3k
  - Fetal
    - train: fetal
    - test: 
## Run 
```
python train_clip.py --data_dir AD --logit_scale 1
```
## Load Pretrain model
The model should saved under `path`/lightning_logs/checkpoints/epochxx-stepxx.ckpt
```
python train_clip.py --data_dir AD --checkpoint model.ckpt
```
The results will be saved under the same `path` 
