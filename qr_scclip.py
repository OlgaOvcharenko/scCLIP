import argparse
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import scanpy as sc
import anndata as ad
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import pandas as pd

import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, f1_score


def get_args():
    parser = argparse.ArgumentParser(description='scCLIP.')

    parser.add_argument('--data', type=str, required=True,
                        help='Adata path')

    args = parser.parse_args()
    return args


args = get_args()
data = args.data
repeat = 0

data_train = f"{data}_train"
data_test = f"{data}_test"


only_files_train = [f'results/{data_train}/1.0_False_30000_0.00015_/{data_train}/{mod}.h5ad' for mod in ["rna", "atac"]]
only_files_test = [f'results/{data_train}/1.0_False_30000_0.00015_/{data_test}/{mod}.h5ad' for mod in ["rna", "atac"]]

# results/simulated_train/1.0_False_30000_0.00015_/simulated_test/rna.h5ad
# results/simulated_test/1.0_False_30000_0.00015_/simulated_test/rna.h5ad

#  concat
if data == "simulated":
    adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_RNA_full.h5ad')
    adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_Protein_full.h5ad')

elif data == "human_multiome":
    adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_GEX_multiome_full.h5ad')
    adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_ATAC_multiome_full.h5ad')

elif data == "human_cite":
    adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_GEX_full.h5ad')
    adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_ADT_full.h5ad')

adata_merged = ad.concat([adata_RNA, adata_Protein], axis=1, merge="same")

adata_RNA_train = sc.read_h5ad(only_files_train[0]) 
adata_Protein_train = sc.read_h5ad(only_files_train[1]) 
adata_train = ad.concat([adata_RNA_train, adata_Protein_train], axis=1, merge="same")

adata_RNA_test = sc.read_h5ad(only_files_test[0]) 
adata_Protein_test = sc.read_h5ad(only_files_test[1]) 
adata_test= ad.concat([adata_RNA_test, adata_Protein_test], axis=1, merge="same")

knn = KNeighborsClassifier(n_neighbors=5)
knn.fit(adata_train, adata_merged.obs["cell_type_l1"].tolist())
cat_preds = knn.predict(adata_test)

cell_types_list = pd.unique(adata_RNA_test.obs['cell_type_l1']).tolist()
acc = accuracy_score(adata_RNA_test.obs['cell_type_l1'].to_list(), cat_preds)
f1 = f1_score(adata_RNA_test.obs['cell_type_l1'].to_list(), cat_preds, labels=cell_types_list, average=None)
f1_weighted = f1_score(adata_RNA_test.obs['cell_type_l1'].to_list(), cat_preds, labels=cell_types_list, average='weighted')
f1_macro = f1_score(adata_RNA_test.obs['cell_type_l1'].to_list(), cat_preds, labels=cell_types_list, average='macro')
f1_median = np.median(f1)

print(f"Per class {cell_types_list} F1 {f1}")
print('Accuracy {:.3f}, F1 median {:.3f}, F1 macro {:.3f}, F1 weighted {:.3f} '.format(acc, f1_median, f1_macro, f1_weighted),)
