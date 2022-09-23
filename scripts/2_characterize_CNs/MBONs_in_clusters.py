# %%
# determine which clusters have aversive and appetitive MBONs
import pymaid
import contools
from contools import Promat, Adjacency_matrix, Celltype_Analyzer, Celltype
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pymaid_creds import url, name, password, token
from data_settings import pairs_path, data_date
rm = pymaid.CatmaidInstance(url, token, name, password)

#pairs = Promat.get_pairs(pairs_path)

MBON = Celltype('MBON', pymaid.get_skids_by_annotation('mw MBON'))
print(MBON)
annots = pymaid.get_annotated('mw brain clusters level 7').name
clusters = [Celltype(annot, pymaid.get_skids_by_annotation(annot)) for annot in annots]

cta = Celltype_Analyzer(clusters)
cta.set_known_types([MBON])
cta.memberships()
# %%
from matplotlib.gridspec import GridSpec

def plot_marginal_cell_type_cluster(size, particular_cell_type, particular_color, cluster_level, path, all_celltypes=None):

    # all cell types plot data
    if(all_celltypes==None):
        _, all_celltypes = Celltype_Analyzer.default_celltypes()
        
    clusters = [Celltype(annot, pymaid.get_skids_by_annotation(annot)) for annot in annots]

    #all_clusters = [Celltype(lvl.cluster_df.cluster[i], lvl.cluster_df.skids[i]) for i in range(0, len(lvl.clusters))]
    cluster_analyze = Celltype_Analyzer(clusters)

    cluster_analyze.set_known_types(all_celltypes)
    celltype_colors = [x.get_color() for x in cluster_analyze.get_known_types()]
    all_memberships = cluster_analyze.memberships()
    all_memberships = all_memberships.iloc[[0,1,2,3,4,5,6,7,8,9,10,11,12,17,13,14,15,16], :] # switching order so unknown is not above outputs and RGNs before pre-outputs
    celltype_colors = [celltype_colors[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11,12,17,13,14,15,16]] # switching order so unknown is not above outputs and RGNs before pre-outputs
    
    # particular cell type data
    cluster_analyze.set_known_types([particular_cell_type])
    membership = cluster_analyze.memberships()

    # plot
    fig = plt.figure(figsize=size) 
    fig.subplots_adjust(hspace=0.1)
    gs = GridSpec(4, 1)

    ax = fig.add_subplot(gs[0:3, 0])
    ind = np.arange(0, len(cluster_analyze.Celltypes))
    ax.bar(ind, membership.iloc[0, :], color=particular_color)
    ax.set(xlim = (-1, len(ind)), ylim=(0,1), xticks=([]), yticks=([]), title=particular_cell_type.get_name())

    ax = fig.add_subplot(gs[3, 0])
    ind = np.arange(0, len(cluster_analyze.Celltypes))
    ax.bar(ind, all_memberships.iloc[0, :], color=celltype_colors[0])
    bottom = all_memberships.iloc[0, :]
    for i in range(1, len(all_memberships.index)):
        plt.bar(ind, all_memberships.iloc[i, :], bottom = bottom, color=celltype_colors[i])
        bottom = bottom + all_memberships.iloc[i, :]
    ax.set(xlim = (-1, len(ind)), ylim=(0,1), xticks=([]), yticks=([]))
    ax.axis('off')
    ax.axis('off')

    plt.savefig(path, format='pdf', bbox_inches='tight')

aversiveMBON = Celltype('MBON-av', pymaid.get_skids_by_annotation('mw MBON subclass_aversive'))
appetitiveMBON = Celltype('MBON-app', pymaid.get_skids_by_annotation('mw MBON subclass_appetitive'))
plot_marginal_cell_type_cluster(size=(4,2), particular_cell_type=MBON, particular_color='gray', cluster_level=7, path='plots/MBON-in-clusters.pdf')
plot_marginal_cell_type_cluster(size=(4,2), particular_cell_type=appetitiveMBON, particular_color='blue', cluster_level=7, path='plots/MBONapp-in-clusters.pdf')
plot_marginal_cell_type_cluster(size=(4,2), particular_cell_type=aversiveMBON, particular_color='red', cluster_level=7, path='plots/MBONav-in-clusters.pdf')

app_val_MBON = Celltype('MBON-app-val', pymaid.get_skids_by_annotation('mw MBON excit Appet') + pymaid.get_skids_by_annotation('mw MBON inhib Avers'))
avers_val_MBON = Celltype('MBON-av-val', pymaid.get_skids_by_annotation('mw MBON excit Avers') + pymaid.get_skids_by_annotation('mw MBON inhib Appet'))
plot_marginal_cell_type_cluster(size=(4,2), particular_cell_type=app_val_MBON, particular_color='blue', cluster_level=7, path='plots/MBON-app-val_in-clusters.pdf')
plot_marginal_cell_type_cluster(size=(4,2), particular_cell_type=avers_val_MBON, particular_color='red', cluster_level=7, path='plots/MBON-av-val_in-clusters.pdf')

# %%
MBONs = pymaid.get_skids_by_annotation('mw MBON')
MBINs = pymaid.get_skids_by_annotation('mw MBIN')
FBNs = pymaid.get_skids_by_annotation('mw FBN')

clusters_MBINs = [x for x in clusters if len(np.intersect1d(x.skids, MBINs))!=0]
clusters_MBONs = [x for x in clusters if len(np.intersect1d(x.skids, MBONs))!=0]
clusters_FBNs = [x for x in clusters if len(np.intersect1d(x.skids, FBNs))!=0]

[pymaid.add_annotations(clusters_MBINs[i].skids, f'mw MBIN-cluster {clusters_MBINs[i].name}') for i in range(len(clusters_MBINs))]
[pymaid.add_annotations(clusters_MBONs[i].skids, f'mw MBON-cluster {clusters_MBONs[i].name}') for i in range(len(clusters_MBONs))]
[pymaid.add_annotations(clusters_FBNs[i].skids, f'mw FBN-cluster {clusters_FBNs[i].name}') for i in range(len(clusters_FBNs))]

# %%
meta = pymaid.get_annotated('MB nomenclature')
pymaid.get_annotations(10682660)

import re

def convert_skids_to_MBnames(ct):
    left = pymaid.get_skids_by_annotation('mw left')
    regex = re.compile('\d+')

    skids = ct.skids
    skid_names = []
    for skid in skids:
        if(skid in left):
            MBnames = np.intersect1d(pymaid.get_annotations(skid)[f'{skid}'], meta.name)
            if(len(MBnames)>0):
                skid_name = ', '.join(MBnames)
            else: 
                skid_name = pymaid.get_neurons(skid).name
                skid_name = skid_name.replace('left', '')

            cluster_name = regex.match(ct.name).group()
            cluster_name = 'cluster_' + cluster_name

            skid_names.append([cluster_name, skid, skid_name])

    return(skid_names)

df = pd.DataFrame([x for sublist in [convert_skids_to_MBnames(ct) for ct in clusters_FBNs] for x in sublist], df = pd.DataFrame([x for sublist in [convert_skids_to_MBnames(ct) for ct in clusters_FBNs] for x in sublist], columns=['cluster', 'skid', 'name'])

# %%
