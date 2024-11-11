import anndata as ad
import scanpy as sc
import mudata as md
from mudata import MuData 
import numpy as np

def create_mudata(data, adata_RNA, adata_Protein, save_path):
    mdata = MuData({"rna": adata_RNA, "atac": adata_Protein})
    mdata.write(f"{save_path}")


# Simulated
# data = 'simulated'
# save_path = f'data/{data}/mudata_simulated_full.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_RNA_full.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_Protein_full.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)

# save_path = f'data/{data}/mudata_simulated_train.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_RNA_train.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_Protein_train.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)

# save_path = f'data/{data}/mudata_simulated_test.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_RNA_test.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_Protein_test.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)


# Multiome
data = 'human_multiome'
# save_path = f'data/{data}/mudata_human_multiome_full.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_GEX_multiome_full.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_ATAC_multiome_full.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)

# save_path = f'data/{data}/mudata_human_multiome_train.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_GEX_multiome_train.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_ATAC_multiome_train.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)

save_path = f'data/{data}/mudata_human_multiome_test.h5mu'
adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_GEX_multiome_test.h5ad')
adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_ATAC_multiome_test.h5ad')
create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)


# Cite
# data = 'human_cite'
# save_path = f'data/{data}/mudata_human_cite_full.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_GEX_full.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_neurips_ADT_full.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)

# save_path = f'data/{data}/mudata_human_cite_train.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_GEX_train.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_ADT_train.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)

# save_path = f'data/{data}/mudata_human_cite_test.h5mu'
# adata_RNA = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_GEX_test.h5ad')
# adata_Protein = sc.read_h5ad('../concerto-reproducibility/Multimodal_pretraining/adata_ADT_test.h5ad')
# create_mudata(data=data, adata_RNA=adata_RNA, adata_Protein=adata_Protein, save_path=save_path)

