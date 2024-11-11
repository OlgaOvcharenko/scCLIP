import argparse
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import scanpy as sc
import anndata as ad
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import pandas as pd


_BIO_METRICS = BioConservation(isolated_labels=True, 
                               nmi_ari_cluster_labels_leiden=True, 
                               nmi_ari_cluster_labels_kmeans=False, 
                               silhouette_label=True, 
                               clisi_knn=True
                               )
_BATCH_METRICS = BatchCorrection(graph_connectivity=True, 
                                 kbet_per_label=True, 
                                 ilisi_knn=True, 
                                 pcr_comparison=True, 
                                 silhouette_batch=True
                                 )


def get_args():
    parser = argparse.ArgumentParser(description='CONCERTO Batch Correction.')

    parser.add_argument('--data', type=str, required=True,
                        help='Adata path')

    args = parser.parse_args()
    return args


def evaluate_model(adata, batch_key="batch", cell_type_label="cell_type_l1"):
    names_obs = ['scclip']
    print(names_obs)
    bm = Benchmarker(
                adata,
                batch_key=batch_key,
                label_key=cell_type_label,
                embedding_obsm_keys=names_obs,
                bio_conservation_metrics=_BIO_METRICS,
                batch_correction_metrics=_BATCH_METRICS,
                n_jobs=4,
            )
    bm.benchmark()
    a = bm.get_results(False, True)
    results = a.round(decimals=4)
    return results

def scale_result(df):
    cols = ['Isolated labels', 'Leiden NMI', 'Leiden ARI', 'Silhouette label', 'cLISI', 'Silhouette batch', 'iLISI', 'KBET', 'Graph connectivity', 'PCR comparison']
    df[cols] = MinMaxScaler().fit_transform(df[cols])
    scaled= pd.DataFrame(df, columns=df.columns, index=df.index)
    
    biometrics = [i for i in ['Isolated labels', 'Leiden NMI', 'Leiden ARI', 'Silhouette label', 'cLISI']]
    batchmetrics = [i for i in ['Silhouette batch', 'iLISI', 'KBET', 'Graph connectivity', 'PCR comparison']]
    scaled[f"Batch correction final"] = scaled[batchmetrics].mean(1)
    scaled[f"Bio conservation final"] = scaled[biometrics].mean(1)
    scaled[f"Total final"] = 0.6 * scaled[f"Bio conservation final"] + 0.4 * scaled[f"Batch correction final"]
    df[f"Bio conservation final"] = scaled[f"Bio conservation final"].copy()
    df[f"Batch correction final"] = scaled[f"Batch correction final"].copy()
    df[f"Total final"] = scaled[f"Total final"].copy()
    return df

args = get_args()
data = args.data
repeat = 0


only_files = [f'results/data/{data}/mudata_{data}_full.h5mu/1.0_False_30000_0.00015_/data/{data}/mudata_{data}_full.h5mu/' + f for f in listdir(f'results/data/{data}/mudata_{data}_full.h5mu/1.0_False_30000_0.00015_/data/{data}/mudata_{data}_full.h5mu/') if isfile(join(f'results/data/{data}/mudata_{data}_full.h5mu/1.0_False_30000_0.00015_/data/{data}/mudata_{data}_full.h5mu/', f)) if f.endswith(f'.h5ad')]

final_df = pd.DataFrame(columns=['combine_omics', 'model_type', 'batch_size', 'epoch', 'lr', 'drop_rate', 'heads', 'Embedding', 'Isolated labels', 'Leiden NMI', 'Leiden ARI', 'Silhouette label', 'cLISI', 'Silhouette batch', 'iLISI', 'KBET', 'Graph connectivity', 'PCR comparison', 'Batch correction', 'Bio conservation', 'Total'])

combine_omics = 0
model_type = 0
batch_size = 256
epoch = 30000
lr = 0.00015
drop_rate = 0.1
heads = 256

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

adata_RNA_new = sc.read_h5ad(only_files[0]) 
adata_Protein_new = sc.read_h5ad(only_files[1]) 
adata = ad.concat([adata_RNA_new, adata_Protein_new], axis=1, merge="same")
adata_merged.obsm["scclip"] = adata.X

df = evaluate_model(adata=adata_merged, cell_type_label="cell_type_l1")
df = df.assign(**{"combine_omics": combine_omics, "model_type": model_type, "batch_size": batch_size, "epoch": epoch, "lr": lr, "drop_rate": drop_rate, "heads": heads})
print(df.columns)

final_df = pd.concat([df, final_df], ignore_index=True)

final_df.to_csv(f'results_downstream/{data}/scclip_new_metrics_unscaled.csv')
