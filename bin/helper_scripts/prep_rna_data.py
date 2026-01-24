from pathlib import Path
import scanpy as sc
import tangram as tg
import pandas as pd
import os

sp_ad = sc.read_h5ad('../data/shareseq/c_sp_ad.h5ad')
sc_ad = sc.read_h5ad('../data/shareseq/c_sc_ad.h5ad')

tg.pp_adatas(sc_ad, sp_ad, genes=sc_ad.var_names)

# Run Tangram mapping (basic example, adjust mode/labels as needed)
tg_map = tg.map_cells_to_space(
   sc_ad,
    sp_ad,
    mode = 'cells', 
    num_epochs=500,
    device="cpu"
)

sc_ad.obs.index = sc_ad.obs_names
common_cells = tg_map.obs_names.intersection(sc_ad.obs_names)
common_cells = sorted(list(common_cells))  # Ensure same order
shared_index = pd.Index(common_cells, dtype="str", name=None)
print(f"Number of common cells for gene projection: {len(shared_index)}")
ad_map_sub = tg_map[shared_index].copy()
ad_sc_sub = sc_ad[shared_index].copy()  # type: ignore 

ad_map_sub.obs_names = shared_index
ad_sc_sub.obs_names = shared_index  # type: ignore 
ad_map_sub.obs.index = shared_index
ad_sc_sub.obs.index = shared_index

ad_ge = tg.project_genes(ad_map_sub, ad_sc_sub)

tg.project_cell_annotations(ad_map_sub, ad_ge, "celltype")
ad_ge.obsm['spatial'] = sp_ad.obsm['spatial']


# Create path object for saving
path = Path('../data/tangram_output')
ad_ge.write_h5ad(path / 'sp_rna.h5ad')