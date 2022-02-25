# %%

import pymaid
import contools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# allows text to be editable in Illustrator
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# font settings
plt.rcParams['font.size'] = 6
plt.rcParams['font.family'] = 'arial'

from pymaid_creds import url, name, password, token
rm = pymaid.CatmaidInstance(url, token, name, password)

PNs = contools.Celltype_Analyzer.get_skids_from_meta_annotation('mw brain PNs')

left_KC = list(np.intersect1d(pymaid.get_skids_by_annotation('mw KC'), pymaid.get_skids_by_annotation('mw left')))
KC_adj = pymaid.adjacency_matrix(pymaid.get_skids_by_annotation('mw KC'), pymaid.get_skids_by_annotation('mw KC'))

KC_no_connections = np.intersect1d(KC_adj.index[KC_adj.sum(axis=0)==0], KC_adj.index[KC_adj.sum(axis=0)==0]) # neurons with no inputs/outputs
KC_adj.drop(index=KC_no_connections, columns=KC_no_connections, inplace=True)

# %%
# cosine similarity of KCs

mat = KC_adj

cols = mat.columns.values
sim_mat = pd.DataFrame(np.zeros(shape=[len(cols), len(cols)]), index=cols, columns=cols)

for i in cols:
    for j in cols:
        a = mat.loc[:, i].values
        b = mat.loc[:, j].values
        dot = np.dot(a, b)
        norma = np.linalg.norm(a)
        normb = np.linalg.norm(b)
        if(norma!=0 and normb!=0):
            calculation = dot / (norma * normb)
        else:
            calculation = 0

        sim_mat.loc[i, j] = calculation

cluster = sns.clustermap(sim_mat)

# %%
# cut dendrogram and collect clusters
### These clusters aren't what we want ####

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import fcluster
from scipy.cluster import hierarchy

cluster_inds = cluster.dendrogram_row.reordered_ind
cluster_dendrogram = cluster.dendrogram_row.dendrogram
cluster_linkage = cluster.dendrogram_row.linkage
D = np.array(cluster_dendrogram['dcoord'])
I = np.array(cluster_dendrogram['icoord'])

max_d = 2.35 # splits into three groups on left and three groups on right
  
Z = cluster_linkage
plt.figure()
plt.title("Dendrograms")
dendrogram = hierarchy.dendrogram(Z)
  
# Cutting the dendrogram at max_d
plt.axhline(y=max_d, c='k')

# determine clusters
cluster_ids = fcluster(cluster_linkage, max_d, criterion='distance')

clusters = []
for i in np.arange(min(cluster_ids), max(cluster_ids), 1):
    cluster_identify = cluster_ids==i
    cluster = list(KC_adj.index[cluster_identify])
    clusters.append(cluster)
    
# doesn't give desired results, but will keep for future reference

# %%
# using Louvain clusters from Kathi's paper
### These annotations aren't what we want ####

group1 = pymaid.get_skids_by_annotation('KC group 1')
group2 = pymaid.get_skids_by_annotation('KC group 2')
group3 = list(np.setdiff1d(pymaid.get_skids_by_annotation('mw KC'), group1+group2))

KC_all = group1 + group2
left_KC = list(np.intersect1d(KC_all, pymaid.get_skids_by_annotation('mw left')))
right_KC = list(np.intersect1d(KC_all, pymaid.get_skids_by_annotation('mw right')))
KC_all = left_KC + right_KC

KC_all_adj = pymaid.adjacency_matrix(KC_all, KC_all)
sns.heatmap(KC_all_adj, vmax = 20)

# doesn't quite add up or make sense; I think I will default to Marta's suggestion

# %%
# sort connectivity matrix as such: KC 1-claw, 2-claw, 3-claw, etc., young

KC1 = pymaid.get_skids_by_annotation('mw KC_subclass_1claw')
KC2 = pymaid.get_skids_by_annotation('mw KC_subclass_2claw')
KC3 = pymaid.get_skids_by_annotation('mw KC_subclass_3claw')
KC4 = pymaid.get_skids_by_annotation('mw KC_subclass_4claw')
KC5 = pymaid.get_skids_by_annotation('mw KC_subclass_5claw')
KC6 = pymaid.get_skids_by_annotation('mw KC_subclass_6claw')

KC = pymaid.get_skids_by_annotation('mw KC')
young_KC = np.setdiff1d(KC, KC1+KC2+KC3+KC4+KC5+KC6)

KC_types = [KC1, KC2, KC3, KC4, KC5, KC6, young_KC]

left = pymaid.get_skids_by_annotation('mw left')
right = pymaid.get_skids_by_annotation('mw right')

KC_types_left = [list(np.intersect1d(skids, left)) for skids in KC_types]
KC_types_right = [list(np.intersect1d(skids, right)) for skids in KC_types]
KC_all_sorted = [x for sublist in KC_types_right for x in sublist] + [x for sublist in KC_types_left for x in sublist]

KC_adj = pymaid.adjacency_matrix(KC_all_sorted, KC_all_sorted)
KC_adj.index.name=''
KC_adj.columns.name=''

# add coloring for labels on side
color1 = '#7c1824'
color2 = '#d3676f'
color3 = '#e5aeb5'

coloring = [[color1]*len(KC_types_right[i]) for i in [0, 1]] + [[color2]*len(KC_types_right[i]) for i in np.arange(2,6)] + [[color3]*len(KC_types_right[6])]
coloring = coloring + [[color1]*len(KC_types_left[i]) for i in [0, 1]] + [[color2]*len(KC_types_left[i]) for i in np.arange(2,6)] + [[color3]*len(KC_types_left[6])]
coloring = [x for sublist in coloring for x in sublist]

# using clustermap simply because it supports row_colors and col_colors
g = sns.clustermap(
    KC_adj, vmax = 10, cmap='Blues', 
    row_colors=coloring, col_colors=coloring, 
    row_cluster=False, col_cluster=False, 
    yticklabels=False, xticklabels=False
)

# Drawing the frame
ax = g.ax_heatmap
ax.axhline(y = 0, color='k',linewidth = 1)
ax.axhline(y = len(coloring), color = 'k',
            linewidth = 1)
  
ax.axvline(x = 0, color = 'k',
            linewidth = 1)
  
ax.axvline(x = len(coloring), 
            color = 'k', linewidth = 1)

# save PDF
plt.savefig('plots/KC-KC_adj.pdf', format='pdf', bbox_inches='tight')

# %%
# plot morphology of these KC types: 1) 1claw + 2claw, 2) 3+ claw, 3) young KC

from contools import Celltype, Celltype_Analyzer
import navis 

KC12 = Celltype(name='KCs 1/2 claws', skids=KC1+KC2, color=color1)
KC3p = Celltype(name='KCs 3+ claws', skids=KC3+KC4+KC5+KC6, color=color2)
KCy = Celltype(name='KCs young', skids=young_KC, color=color3)

volume = 'PS_Neuropil_manual'
neuropil = pymaid.get_volume(volume)
vol_color = (250, 250, 250, .05)
neuropil.color = vol_color

neurons = pymaid.get_neurons(KC12.skids)

color=color1
alpha = 0.5
linewidth=1.5
connectors=False

fig, ax = navis.plot2d([neurons, neuropil], method='3d_complex', color=color, linewidth=linewidth, connectors=connectors, cn_size=2, alpha=alpha)
ax.azim=-90
ax.elev=-90
ax.dist=3
ax.xlim3d=(-4500, 110000)
ax.ylim3d=(-4500, 110000)

# used manually generated images in the end anyways

# %%
