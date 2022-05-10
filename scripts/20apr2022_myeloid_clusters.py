#!/usr/bin/env python
# coding: utf-8

# # wrapping up medullo for 440
# ## what we know so far
# 1. SPP1 should be added to our assay and there's a really solid Ab for it
# 2. Batch effect normalization is worth paying attention, maybe do analysis w/ and w/o
# 3. myeloid vs lymphocyte classifications are valid, the rest are not
# 
# ## what to do
# 1. recluster all myeloid with harmony and without
# 2. check for representation across tumor subgroups
# 3. identify interesting clusters
# 4. try to classify what those clusters represent
# 5. use different samples to get statistical significance

# In[2]:


import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from scipy import stats
print(np.__version__)
print(pd.__version__)
print(sc.__version__)
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
results_file = 'write/440.h5ad'


# ## stat testing
# bleh

# In[3]:


tcell_counts = pd.read_csv('compmyeloid_cluster_counts.csv', index_col=0)
tcell_counts


# In[4]:


adata = sc.read(
    "medullo_immune_cells/matrix.csv", 
    cache=True)
non_normalized = [data for data in adata.var.index if 'ALRA_' not in data]
nn_adata = adata[:, non_normalized]

meta_df = pd.read_csv("meta.tsv", delimiter="\t", index_col = 0)
getem = ['cell_types','clusters','subgroup']
w_meta = pd.merge(nn_adata.obs, meta_df[getem], right_index=True, left_index=True)
nn_adata.obs = w_meta
nn_adata.obs['batch'] = meta_df['UPN'].astype("category")

sc.pp.log1p(nn_adata)
sc.pp.highly_variable_genes(nn_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)#, batch_key='batch')

sc.pp.filter_cells(nn_adata, min_genes=200)
sc.pp.filter_genes(nn_adata, min_cells=3)
nn_adata.var['mt'] = nn_adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(nn_adata, 
                           qc_vars=['mt'], 
                           percent_top=None, 
                           log1p=False, 
                           inplace=True)


# In[5]:


nn_adata = nn_adata[nn_adata.obs.n_genes_by_counts < 5000, :]
nn_adata = nn_adata[nn_adata.obs.pct_counts_mt < 5.5, :]
sc.pp.normalize_total(nn_adata, target_sum=1e4)
nn_adata.raw = nn_adata
nn_adata = nn_adata[:, nn_adata.var.highly_variable]
sc.pp.regress_out(nn_adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(nn_adata, max_value=10)


# In[6]:


sc.tl.pca(nn_adata)#, svd_solver='arpack')
sc.pl.pca(nn_adata, color='cell_types', save=".png")
sce.pp.harmony_integrate(nn_adata, 'batch')


# In[60]:


smol_df = nn_adata.obs[['subgroup','batch']].sort_values('batch')
set([(smol_df.loc[row][0], smol_df.loc[row][1]) for row in smol_df.index])
    


# In[7]:


sc.pp.neighbors(nn_adata, use_rep='X_pca', n_pcs = 12)#, n_neighbors = 20)
sc.tl.umap(nn_adata)
sc.set_figure_params(figsize=(10,6))
sc.tl.leiden(nn_adata)
# sc.pl.umap(nn_adata, color='leiden')#, save="_clustered.png")


# In[65]:


# sc.pl.umap(nn_adata, color='cell_types')#, save="_clustered.png")


# In[198]:


clusters = list(set(nn_adata.obs["cell_types"]))
batch_to_sg = {'925': 'GP3', '945': 'GP3', '1028': 'GP3', '1130': 'GP3',
               '1167': 'GP3', '1355': 'GP3', '1433': 'GP3', '753': 'GP4',
               '934': 'GP4', '966': 'GP4', '996': 'GP4', '1066': 'GP4',
               '1070': 'GP4', '1125': 'GP4', '1155': 'GP4', '1177': 'GP4',
               '1195': 'GP4', '1238': 'GP4', '801': 'SHH', '831': 'SHH',
               '877': 'SHH', '898': 'SHH', '1224': 'SHH', '1235': 'SHH',
               '1325': 'SHH', '1397': 'SHH', '1416': 'SHH'
              }
batch_dict = {}
# what order are clusters in? 
for b in set(nn_adata.obs["batch"]):
    batch_df = nn_adata[nn_adata.obs["batch"] == b]
    batch_dict[b] = [len(batch_df[batch_df.obs['cell_types'] == c])/len(batch_df) for c in clusters]

statstest = pd.DataFrame(index=["GP3/GP4 T val","GP3/GP4 p val", "GP4/SHH T val","GP4/SHH p val", "GP3/SHH T val","GP3/SHH p val"])
for index in range(len(clusters)):
    sg_dict = {"GP3":[],"GP4":[],"SHH":[]}
    col = []
    for b in batch_dict:
        sg_dict[batch_to_sg[str(b)]].append(batch_dict[b][index])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["GP4"])[0])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["GP4"])[1])
    col.append(stats.ttest_ind(sg_dict["GP4"], sg_dict["SHH"])[0])
    col.append(stats.ttest_ind(sg_dict["GP4"], sg_dict["SHH"])[1])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["SHH"])[0])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["SHH"])[1])
    statstest[clusters[index]] = col
    if clusters[index] == '10':
        print(sg_dict)
statstest = statstest.sort_index(axis=1)
statstest.to_csv("tvals_T_cells.csv")
print(batch_dict.keys())
statstest


# In[8]:


myeloid_labels = ['Complement myeloid (Complement-M)', 
                  'Dendritic cell-like myeloid (DC-M)', 
                  'M2-activated myeloid  (M2-M)',
                  'Chemokine myeloid (Chemokine-M)',
                  'Neutrophil (Nt)',
                  'Non-activated microglia (NA-Microglia)'
                 ]

myeloid_bools = [True if x in myeloid_labels else False for x in nn_adata.obs['cell_types']]
all_myeloid = nn_adata[myeloid_bools]

sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(all_myeloid, )#, n_pcs = 12, n_neighbors = 20)
sc.tl.leiden(all_myeloid)#, 1.1)
sc.tl.umap(all_myeloid)
sc.pl.umap(all_myeloid, color='leiden', save="_USE_all_myeloid_leiden.png", title="All Myeloid cells clustered")


# In[9]:


sc.pl.umap(all_myeloid, color='SPP1', save="all_myeloid_leiden.png", title="All Myeloid cells clustered")


# In[8]:


sc.pl.umap(all_myeloid, color='CD163', save="all_myeloid_leiden.png", title="All Myeloid cells clustered")


# In[110]:


sc.pl.umap(all_myeloid, color='cell_types', title="All myeloid cells by cell type")


# In[103]:


tcell_columns = list(set(all_myeloid.obs['leiden']))
tcell_rows = list(set(all_myeloid.obs['subgroup']))
tcell_counts = pd.DataFrame(columns=tcell_columns, index=tcell_rows)
for ct in tcell_columns:
    for sg in tcell_rows:
        new_data = all_myeloid[all_myeloid.obs['subgroup'] == sg]
        new_data = new_data[new_data.obs['leiden'] == ct]
        tcell_counts.loc[sg,ct] = len(new_data)
tcell_counts.to_csv('allmyeloid_cluster_counts.csv') 

sums = []
for c in tcell_counts.index:
    sums.append(sum(tcell_counts.loc[c]))
    
tcell_counts['sums'] = sums
for c in tcell_counts.index:
    for r in tcell_counts.columns[:-1]:
        tcell_counts.loc[c, r] = tcell_counts.loc[c, r]/tcell_counts.loc[c,'sums']
   
tcell_counts = tcell_counts.sort_index(axis=1)
# tcell_counts.to_csv('normalized_allmyeloid_cluster_counts.csv') 
tcell_counts


# In[163]:


clusters = list(set(all_myeloid.obs["leiden"]))
batch_to_sg = {'925': 'GP3', '945': 'GP3', '1028': 'GP3', '1130': 'GP3',
               '1167': 'GP3', '1355': 'GP3', '1433': 'GP3', '753': 'GP4',
               '934': 'GP4', '966': 'GP4', '996': 'GP4', '1066': 'GP4',
               '1070': 'GP4', '1125': 'GP4', '1155': 'GP4', '1177': 'GP4',
               '1195': 'GP4', '1238': 'GP4', '801': 'SHH', '831': 'SHH',
               '877': 'SHH', '898': 'SHH', '1224': 'SHH', '1235': 'SHH',
               '1325': 'SHH', '1397': 'SHH', '1416': 'SHH'
              }
batch_dict = {}
# what order are clusters in? 
for b in set(all_myeloid.obs["batch"]):
    batch_df = all_myeloid[all_myeloid.obs["batch"] == b]
    batch_dict[b] = [len(batch_df[batch_df.obs['leiden'] == c])/len(batch_df) for c in clusters]

statstest = pd.DataFrame(index=["GP3/GP4 T val","GP3/GP4 p val", "GP4/SHH T val","GP4/SHH p val", "GP3/SHH T val","GP3/SHH p val"])
for index in range(len(clusters)):
    sg_dict = {"GP3":[],"GP4":[],"SHH":[]}
    col = []
    for b in batch_dict:
        sg_dict[batch_to_sg[str(b)]].append(batch_dict[b][index])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["GP4"])[0])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["GP4"])[1])
    col.append(stats.ttest_ind(sg_dict["GP4"], sg_dict["SHH"])[0])
    col.append(stats.ttest_ind(sg_dict["GP4"], sg_dict["SHH"])[1])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["SHH"])[0])
    col.append(stats.ttest_ind(sg_dict["GP3"], sg_dict["SHH"])[1])
    statstest[clusters[index]] = col
    if clusters[index] == '10':
        print(sg_dict)
statstest = statstest.sort_index(axis=1)
# statstest.to_csv("CORRECTED_statstestdf.csv")
print(batch_dict.keys())
statstest


# ### depleted in SHH
# groups 10, 2 and 6

# ## which are statistically significant
# for each batch, 

# In[147]:


sc.tl.rank_genes_groups(all_myeloid, "leiden", method='t-test')
sc.pl.rank_genes_groups(all_myeloid, n_genes=15, sharey=False, ncols=4, fontsize=12, save="all_myeloid_top_genes.png")


# In[112]:


cluster_cat = [int(l) for l in all_myeloid.obs['leiden']]
set(cluster_cat)
all_myeloid.obs["cluster cat"] = cluster_cat
new_data = all_myeloid[all_myeloid.obs['cluster cat'] == 10]
print(set(new_data.obs['clusters']))
sc.pl.umap(new_data, color='subgroup', size=95, 
           title="Cluster 10 by subgroup", 
           save="_myeloid_clus10.png")


# In[166]:


new_data = all_myeloid[all_myeloid.obs['leiden'] == '10']
sc.pl.umap(new_data, color='batch', size=95, 
           title="Cluster 10 by batch") 
#            save="_myeloid_clus10.png")


# ### cluster 3

# In[113]:


cluster_cat = [int(l) for l in all_myeloid.obs['leiden']]
set(cluster_cat)
all_myeloid.obs["cluster cat"] = cluster_cat
new_data = all_myeloid[all_myeloid.obs['cluster cat'] == 3]
print(set(new_data.obs['clusters']))
sc.pl.umap(new_data, color='subgroup', size=95, title="Cluster 3 by subgroup", save="_myeloid_clus3.png") #save="_all_myeloid_leiden.png", title="All Myeloid cells clustered")


# In[167]:


new_data = all_myeloid[all_myeloid.obs['leiden'] == '3']
sc.pl.umap(new_data, color='batch', size=95, 
           title="Cluster 3 by batch", 
           save="_myeloid_clus3_batch.png")


# In[159]:


new_data = all_myeloid[all_myeloid.obs['leiden'] == '5']
sc.pl.umap(new_data, color='subgroup', size=95, title="Cluster 5 by subgroup", save="_myeloid_clus5.png") #save="_all_myeloid_leiden.png", title="All Myeloid cells clustered")


# In[160]:


new_data = all_myeloid[all_myeloid.obs['leiden'] == '1']
sc.pl.umap(new_data, color='subgroup', size=95, title="Cluster 1 by subgroup", save="_myeloid_clus1.png") #save="_all_myeloid_leiden.png", title="All Myeloid cells clustered")


# ## HARMONy

# In[118]:


B_cell = nn_adata[nn_adata.obs["cell_types"] == 'B cell (B)']
B_cell


# In[125]:


sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(B_cell)
sc.tl.leiden(B_cell)
sc.tl.umap(B_cell)
sc.pl.umap(B_cell, color='leiden', save="bcells_leiden.png", title="B cells clustered", size=300)


# In[126]:


sc.pl.umap(B_cell, color='subgroup', save="bcells_sg.png", title="B cells by subgroup", size=300)


# In[124]:


sc.tl.rank_genes_groups(B_cell, "leiden", method='t-test')
sc.pl.rank_genes_groups(B_cell, n_genes=15, sharey=False, ncols=3, fontsize=12, save="all_bcell_top_genes.png")


# In[127]:


prolif = nn_adata[nn_adata.obs["cell_types"] == 'Proliferative (Prolif.)']
prolif


# In[128]:


sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(prolif)
sc.tl.leiden(prolif)
sc.tl.umap(prolif)
sc.pl.umap(prolif, color='leiden', save="_prolif_leiden.png", title="Proliferative cells clustered", size=300)


# In[131]:


sc.pl.umap(prolif, color='subgroup', save="prolif_sg.png", title="Proliferative cells by subgroup", size=300)


# In[134]:


treg = nn_adata[nn_adata.obs["cell_types"] == 'Treg cell (Treg)']
sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(treg)
sc.tl.leiden(treg)
sc.tl.umap(treg)
sc.pl.umap(treg, color='leiden', save="_treg_leiden.png", title="Treg cells clustered", size=300)


# In[136]:


sc.pl.umap(prolif, color='subgroup', save="_treg_sg.png", title="Proliferative cells by subgroup", size=300)


# In[137]:


nk = nn_adata[nn_adata.obs["cell_types"] == 'NK cell (NK)']
sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(nk)
sc.tl.leiden(nk)
sc.tl.umap(nk)
sc.pl.umap(nk, color='leiden', save="_nk_leiden.png", title="NK cells clustered", size=300)


# In[146]:


sc.pl.umap(nk, color='subgroup', save="_nk_sg.png", title="NK cells by subgroup", size=300)


# In[188]:


dcm = nn_adata[nn_adata.obs["cell_types"] == 'Dendritic cell-like myeloid (DC-M)']
sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(dcm)
sc.tl.leiden(dcm, 0.5)
sc.tl.umap(dcm)
sc.pl.umap(dcm, color='leiden', save="_dcm_leiden.png", title="DC-M cells clustered", size=300)


# In[189]:


sc.pl.umap(dcm, color='subgroup', save="_dcm_sg.png", title="DC-M cells by subgroup", size=300)


# In[196]:


sc.pl.umap(dcm, color='batch', size=300, cmap="viridis",
           save="_dcm_batch.png", 
           title="DC-M cells by subgroup"
          )
dcm_gp3 = dcm[dcm.obs["subgroup"] == 'GP3']
sc.pl.umap(dcm_gp3, color='batch', size=300, cmap="viridis",
           save="_dcm_gp3_batch.png", 
           title="DC-M cells by subgroup"
          )


# In[145]:


sc.tl.rank_genes_groups(dcm, "leiden", method='t-test')
sc.pl.rank_genes_groups(dcm, n_genes=15, sharey=False, ncols=3, fontsize=12, save="all_bcell_top_genes.png")


# In[171]:


dcm = nn_adata[nn_adata.obs["cell_types"] == 'Non-activated microglia (NA-Microglia)']
sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(dcm)
sc.tl.leiden(dcm, 0.5)
sc.tl.umap(dcm)
sc.pl.umap(dcm, color='leiden', size=300,
           save="_namicroglia_leiden.png", 
           title="Non-activated microglia cells clustered"
          )


# In[172]:


sc.pl.umap(dcm, color='subgroup', size=300,
           save="_namicroglia_sg.png", 
           title="Non-activated microglia cells by subgroup"
          )


# In[175]:


sc.pl.umap(dcm, color='batch', size=300,
#            save="_namicroglia_sg.png", 
           title="Non-activated microglia cells by batch"
          )


# In[174]:


sc.tl.rank_genes_groups(dcm, "leiden", method='t-test')
sc.pl.rank_genes_groups(dcm, n_genes=15, sharey=False, 
                        ncols=2, fontsize=12, 
                        save="all_namicroglia_top_genes.png")


# In[176]:


m2m = nn_adata[nn_adata.obs["cell_types"] == 'M2-activated myeloid  (M2-M)']
sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(m2m)
sc.tl.leiden(m2m, 0.5)
sc.tl.umap(m2m)
sc.pl.umap(m2m, color='leiden', size=300,
           save="_m2m_leiden.png", 
           title="M2-activated myeloid cells clustered"
          )


# In[177]:


sc.pl.umap(m2m, color='subgroup', size=300,
           save="_m2m_sg.png", 
           title="M2-activated myeloid cells by subgroup"
          )


# In[178]:


sc.pl.umap(m2m, color='batch', size=300,
#            save="_m2m_leiden.png", 
           title="M2-activated myeloid cells by batch"
          )


# In[184]:


nt = nn_adata[nn_adata.obs["cell_types"] == 'Neutrophil (Nt)']
sc.set_figure_params(figsize=(8,6))
sc.pp.neighbors(nt)
sc.tl.leiden(nt)
sc.tl.umap(nt)
sc.pl.umap(nt, color='leiden', size=300,
           save="_nt_leiden.png", 
           title="Neutrophil cells clustered"
          )


# In[185]:


sc.pl.umap(nt, color='subgroup', size=300,
           save="_nt_sg.png", 
           title="Neutrophil cells by subgroup"
          )


# In[187]:


sc.pl.umap(nt, color='batch', size=300,
           save="_nt_batch.png", 
           title="Neutrophil cells by batch"
          )


# In[ ]:




