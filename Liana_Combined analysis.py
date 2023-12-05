#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import liana
import liana as li
# needed for visualization and toy data
import scanpy as sc
import anndata
# import liana's rank_aggregate
from liana.mt import rank_aggregate
# import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean


# In[2]:


file_path = "/home/guest/Disha/Liana code/ColonscRNAseq_pro_combined.h5ad"
adata = anndata.read_h5ad(file_path)


# In[3]:


adata


# In[4]:


sc.pl.umap(adata, color=['Cell_type'])


# In[5]:


adata.raw = adata.copy()


# In[6]:


adata.raw.X


# In[7]:


# Run rank_aggregate
li.mt.rank_aggregate(adata, groupby='Cell_type', expr_prop=0, verbose=True)


# In[ ]:


Aggregate_combined_low_highIP_plot = li.pl.dotplot(adata = adata,
                        colour='magnitude_rank',
                        inverse_colour=True,
                        size='specificity_rank',
                        inverse_size=True,
                        source_labels=['B cells', 'LowIP epithelial cells', 'HighIP epithelial cells', 'T cells', 'Myeloids', 'Stromal cells'],
                        target_labels=['B cells', 'LowIP epithelial cells', 'HighIP epithelial cells', 'T cells', 'Myeloids', 'Stromal cells'],
                        filterby='specificity_rank',
                        filter_lambda=lambda x: x <= 0.01,
                        figure_size=(15,180),
                       )
Aggregate_combined_low_highIP_plot


# In[8]:


import pandas as pd
adata_liana_res = adata.uns['liana_res']


# In[9]:


columns_to_keep = ['source', 'target', 'specificity_rank', 'magnitude_rank', 'ligand_complex', 'receptor_complex']


# In[10]:


filtered_adata_Liana_res = pd.DataFrame(adata_liana_res, columns=columns_to_keep)


# In[11]:


filtered_adata_Liana_res


# In[12]:


filtered_adata_Liana_res = filtered_adata_Liana_res[filtered_adata_Liana_res['specificity_rank'] <= 0.05]


# In[13]:


filtered_adata_Liana_res = pd.DataFrame(filtered_adata_Liana_res)


# In[14]:


filtered_adata_Liana_res['ligand_receptor_pair'] = filtered_adata_Liana_res['ligand_complex'] + '-' + filtered_adata_Liana_res['receptor_complex']


# In[15]:


filtered_adata_Liana_res


# In[16]:


import numpy as np
filtered_adata_Liana_res['magnitude_rank'] = -np.log10(filtered_adata_Liana_res['magnitude_rank'])


# In[17]:


filtered_adata_Liana_res


# In[110]:


ST_cells_filtered_adata_Liana_res = filtered_adata_Liana_res[filtered_adata_Liana_res['source'] == 'T cells']


# In[111]:


import matplotlib.pyplot as plt
import seaborn as sns
heatmap_data = ST_cells_filtered_adata_Liana_res.pivot_table(index='ligand_receptor_pair', columns='target', values='magnitude_rank')

# Create the tile plot
plt.figure(figsize=(12, 24))
sns.set(font_scale=1)
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Magnitude Rank'})

plt.title('Ligand-Receptor Interaction Tile Plot with T cells as source')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.savefig('heatmap_plot.pdf', format='pdf')
plt.show()


# In[21]:


Epi_cells_filtered_Liana_res = filtered_adata_Liana_res[filtered_adata_Liana_res['source'] == 'Epithelial cells']
Epi_cells_filtered_Liana_res = Epi_cells_filtered_Liana_res[Epi_cells_filtered_Liana_res['target'] == 'T cells']
Epi_cells_filtered_Liana_res = Epi_cells_filtered_Liana_res[Epi_cells_filtered_Liana_res['magnitude_rank'] >= 1.5]
Epi_cells_filtered_Liana_res


# In[25]:


import matplotlib.pyplot as plt
import seaborn as sns
heatmap_data = Epi_cells_filtered_Liana_res.pivot_table(index='ligand_receptor_pair', columns='target', values='magnitude_rank')

# Create the tile plot
plt.figure(figsize=(6, 9))
sns.set(font_scale=1)
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Magnitude Rank'})

plt.title('Ligand-Receptor Interaction Tile Plot with Epithelial cells as source')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.savefig('heatmap_plot.pdf', format='pdf')
plt.show()


# In[60]:


SLEpi_cells_filtered_Liana_res = filtered_adata_Liana_res[filtered_adata_Liana_res['source'] == 'LowIP epithelial cells']
SLEpi_cells_filtered_Liana_res = SLEpi_cells_filtered_Liana_res[SLEpi_cells_filtered_Liana_res['target'] == 'T cells']
SLEpi_cells_filtered_Liana_res


# In[29]:


heatmap_data = SLEpi_cells_filtered_Liana_res.pivot_table(index='ligand_receptor_pair', columns='target', values='magnitude_rank')

# Create the tile plot
plt.figure(figsize=(15, 54))
sns.set(font_scale=1)
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Magnitude Rank'})

plt.title('Ligand-Receptor Interaction Tile Plot with LowIP Epithelial cells as source')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.savefig('heatmap_plot.pdf', format='pdf')
plt.show()


# In[61]:


combined_data = pd.concat([SLEpi_cells_filtered_Liana_res, SHEpi_cells_filtered_Liana_res])
combined_data


# In[65]:


heatmap_data = combined_data.pivot(index='ligand_receptor_pair', columns='source', values='magnitude_rank')

# Create the heatmap
plt.figure(figsize=(12, 45))
sns.heatmap(heatmap_data, cmap="YlGnBu", annot=True, fmt=".2f", cbar_kws={'label': 'Magnitude'})
plt.title('Combined Ligand-Receptor Interaction Heatmap')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.show()

