#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#This file consists of code to run Liana using the LowIP and highIP groups
#Consensus method (aggregate) of liana is first used to perform the cell-cell interaction analysis
#Dot plots can  be conputed after getting the Liana output
#This file also contains codes to compute tileplots using the Liana output and by using Tcells or Epithelial cells as the source


# In[205]:


pip install liana


# In[206]:


# import liana
import liana as li
# needed for visualization and toy data
import scanpy as sc
# needed for AnnData analysis
import anndata


# In[207]:


# HighIP data analysis
file_path = "/home/guest/Disha/Liana code/ColonscRNAseq_pro_HighIP.h5ad"
adata_HighIP = anndata.read_h5ad(file_path)


# In[208]:


adata_HighIP


# In[209]:


sc.pl.umap(adata_HighIP, color='Cell_type', title='', frameon=False)


# In[210]:


adata_HighIP.raw = adata_HighIP.copy()


# In[211]:


adata_HighIP.raw.X


# In[212]:


# import liana's rank_aggregate
from liana.mt import rank_aggregate


# In[213]:


# import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean


# In[214]:


# Run rank_aggregate
li.mt.rank_aggregate(adata_HighIP, groupby='Cell_type', expr_prop=0, verbose=True)


# In[215]:


adata_HighIP.uns['liana_res'].head()


# In[53]:


Aggregate_HighIP_plot = li.pl.dotplot(adata = adata_HighIP,
                        colour='magnitude_rank',
                        inverse_colour=True,
                        size='specificity_rank',
                        inverse_size=True,
                        source_labels=['B cells', 'Epithelial cells', 'T cells', 'Myeloids', 'Stromal cells'],
                        target_labels=['B cells', 'Epithelial cells', 'T cells', 'Myeloids', 'Stromal cells'],
                        filterby='specificity_rank',
                        filter_lambda=lambda x: x <= 0.01,
                        figure_size=(12,160),
                       )
Aggregate_HighIP_plot


# In[172]:


# LowIP data analysis
file_path = "/home/guest/Disha/Liana code/ColonscRNAseq_pro_LowIP.h5ad"
adata_LowIP = anndata.read_h5ad(file_path)


# In[216]:


adata_LowIP


# In[217]:


sc.pl.umap(adata_LowIP, color='Cell_type', title='UMAP_LowIP', frameon=False)


# In[218]:


adata_LowIP.raw = adata_LowIP.copy()


# In[219]:


adata_LowIP.raw.X


# In[66]:


# import liana's rank_aggregate
from liana.mt import rank_aggregate


# In[67]:


# import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean


# In[220]:


# Run rank_aggregate
li.mt.rank_aggregate(adata_LowIP, groupby='Cell_type', expr_prop=0, verbose=True)


# In[69]:


Aggregate_LowIP_plot = li.pl.dotplot(adata = adata_LowIP,
                        colour='magnitude_rank',
                        inverse_colour=True,
                        size='specificity_rank',
                        inverse_size=True,
                        source_labels=['B cells', 'Epithelial cells', 'T cells', 'Myeloids', 'Stromal cells'],
                        target_labels=['B cells', 'Epithelial cells', 'T cells', 'Myeloids', 'Stromal cells'],
                        filterby='specificity_rank',
                        filter_lambda=lambda x: x <= 0.01,
                        figure_size=(12,160),
                       )
Aggregate_LowIP_plot


# In[105]:


get_ipython().run_line_magic('pinfo', 'li.pl.dotplot')


# In[221]:


adata_HighIP.uns['liana_res']


# In[222]:


import pandas as pd
HighIP_liana_res = adata_HighIP.uns['liana_res']


# In[223]:


columns_to_keep = ['source', 'target', 'specificity_rank', 'magnitude_rank', 'ligand_complex', 'receptor_complex']


# In[224]:


filtered_HighIP_Liana_res = pd.DataFrame(HighIP_liana_res, columns=columns_to_keep)


# In[225]:


filtered_HighIP_Liana_res


# In[227]:


filtered_HighIP_Liana_res = filtered_HighIP_Liana_res[filtered_HighIP_Liana_res['specificity_rank'] <= 0.01]


# In[228]:


filtered_HighIP_Liana_res


# In[229]:


filtered_HighIP_Liana_res = pd.DataFrame(filtered_HighIP_Liana_res)


# In[230]:


filtered_HighIP_Liana_res['ligand_receptor_pair'] = filtered_HighIP_Liana_res['ligand_complex'] + '-' + filtered_HighIP_Liana_res['receptor_complex']


# In[231]:


filtered_HighIP_Liana_res


# In[232]:


import numpy as np
filtered_HighIP_Liana_res['magnitude_rank'] = -np.log10(filtered_HighIP_Liana_res['magnitude_rank'])


# In[233]:


filtered_HighIP_Liana_res


# In[234]:


ST_cells_filtered_HighIP_Liana_res = filtered_HighIP_Liana_res[filtered_HighIP_Liana_res['source'] == 'T cells']


# In[235]:


ST_cells_filtered_HighIP_Liana_res


# In[249]:


import matplotlib.pyplot as plt
import seaborn as sns
heatmap_data = ST_cells_filtered_HighIP_Liana_res.pivot_table(index='ligand_receptor_pair', columns='target', values='magnitude_rank')

# Create the tile plot
plt.figure(figsize=(12, 24))
sns.set(font_scale=1)
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Magnitude Rank'})

plt.title('Ligand-Receptor Interaction Tile Plot with T cells as source_HighIP')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.savefig('heatmap_plot.pdf', format='pdf')
plt.show()


# In[237]:


LowIP_liana_res = adata_LowIP.uns['liana_res']


# In[238]:


columns_to_keep = ['source', 'target', 'specificity_rank', 'magnitude_rank', 'ligand_complex', 'receptor_complex']


# In[239]:


filtered_LowIP_Liana_res = pd.DataFrame(LowIP_liana_res, columns=columns_to_keep)
filtered_LowIP_Liana_res


# In[240]:


filtered_LowIP_Liana_res = filtered_LowIP_Liana_res[filtered_LowIP_Liana_res['specificity_rank'] <= 0.01]
filtered_LowIP_Liana_res


# In[241]:


filtered_LowIP_Liana_res = pd.DataFrame(filtered_LowIP_Liana_res)
filtered_LowIP_Liana_res['ligand_receptor_pair'] = filtered_LowIP_Liana_res['ligand_complex'] + '-' + filtered_LowIP_Liana_res['receptor_complex']


# In[242]:


import numpy as np
filtered_LowIP_Liana_res['magnitude_rank'] = -np.log10(filtered_LowIP_Liana_res['magnitude_rank'])
filtered_LowIP_Liana_res


# In[243]:


ST_cells_filtered_LowIP_Liana_res = filtered_LowIP_Liana_res[filtered_LowIP_Liana_res['source'] == 'T cells']


# In[250]:


heatmap_data = ST_cells_filtered_LowIP_Liana_res.pivot_table(index='ligand_receptor_pair', columns='target', values='magnitude_rank')

# Create the tile plot
plt.figure(figsize=(12, 24))
sns.set(font_scale=1)
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Magnitude Rank'})

plt.title('Ligand-Receptor Interaction Tile Plot with T cells as source_LowIP')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.savefig('heatmap_plot.pdf', format='pdf')
plt.show()


# In[245]:


SEpi_cells_filtered_HighIP_Liana_res = filtered_HighIP_Liana_res[filtered_HighIP_Liana_res['source'] == 'Epithelial cells']
SEpi_cells_filtered_HighIP_Liana_res


# In[251]:


heatmap_data = SEpi_cells_filtered_HighIP_Liana_res.pivot_table(index='ligand_receptor_pair', columns='target', values='magnitude_rank')

# Create the tile plot
plt.figure(figsize=(15, 54))
sns.set(font_scale=1)
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Magnitude Rank'})

plt.title('Ligand-Receptor Interaction Tile Plot with Epithelial cells as source_HighIP')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.savefig('heatmap_plot.pdf', format='pdf')
plt.show()


# In[247]:


SEpi_cells_filtered_LowIP_Liana_res = filtered_LowIP_Liana_res[filtered_LowIP_Liana_res['source'] == 'Epithelial cells']
SEpi_cells_filtered_LowIP_Liana_res


# In[252]:


heatmap_data = SEpi_cells_filtered_LowIP_Liana_res.pivot_table(index='ligand_receptor_pair', columns='target', values='magnitude_rank')

# Create the tile plot
plt.figure(figsize=(15, 54))
sns.set(font_scale=1)
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Magnitude Rank'})

plt.title('Ligand-Receptor Interaction Tile Plot with Epithelial cells as source_LowIP')
plt.xlabel('Target Cells')
plt.ylabel('Ligand-Receptor Pair')
plt.savefig('heatmap_plot.pdf', format='pdf')
plt.show()

