# %%

import pymaid
import contools
from contools import Promat, Adjacency_matrix, Celltype_Analyzer, Celltype
import numpy as np
import pandas as pd

from pymaid_creds import url, name, password, token
from data_settings import data_date, pairs_path
rm = pymaid.CatmaidInstance(url, token, name, password)

pairs = Promat.get_pairs(pairs_path)

# %%
# prep data

# load input counts
inputs = pd.read_csv(f'data/adj/inputs_{data_date}.csv', index_col=0)
outputs = pd.read_csv(f'data/adj/outputs_{data_date}.csv', index_col=0)

# load appropriate sensory and ascending types
modalities = 'mw brain sensory modalities'
input_types, input_names = Celltype_Analyzer.get_skids_from_meta_meta_annotation(modalities, split=True)
input_names = [annot.replace('mw ', '') for annot in input_names]
input_types = pd.DataFrame(zip(input_names, input_types), columns = ['type', 'source'])

# %%
# identify all 2nd/3rd/4th-order neurons

threshold = 0.01
edges_ad = pd.read_csv(f'data/edges_threshold/ad_pairwise-input-threshold-{threshold}_all-edges_{data_date}.csv', index_col=0)
edges_ad_pairs = pd.read_csv(f'data/edges_threshold/ad_pairwise-input-threshold-{threshold}_paired-edges_{data_date}.csv', index_col=0)

modalities = 'mw brain sensory modalities'
brain_inputs = Celltype_Analyzer.get_skids_from_meta_meta_annotation(modalities, split=False)
brain_inputs = brain_inputs + pymaid.get_skids_by_annotation('mw A1 ascending unknown')
brain = pymaid.get_skids_by_annotation('mw brain neurons') + brain_inputs

pdiff = pymaid.get_skids_by_annotation('mw partially differentiated')
SEZ_motorneurons = pymaid.get_skids_by_annotation('mw motor')
exclude = brain_inputs+pdiff+SEZ_motorneurons

order2 = [Promat.downstream_multihop(edges=edges_ad, sources=skids, hops=1, exclude=exclude, exclude_unpaired=True, pairs=pairs)[0] for skids in input_types.source]
input_types['order2'] = order2 

all_order2 = list(np.unique([x for sublist in input_types.order2 for x in sublist]))
order3 = [Promat.downstream_multihop(edges=edges_ad, sources=skids, hops=1, exclude=exclude+all_order2, exclude_unpaired=True, pairs=pairs)[0] for skids in input_types.order2]
input_types['order3'] = order3

all_order3 = list(np.unique([x for sublist in input_types.order3 for x in sublist]))
order4 = [Promat.downstream_multihop(edges=edges_ad, sources=skids, hops=1, exclude=exclude+all_order3+all_order2, exclude_unpaired=True, pairs=pairs)[0] for skids in input_types.order3]
input_types['order4'] = order4

all_order4 = list(np.unique([x for sublist in input_types.order4 for x in sublist]))
order5 = [Promat.downstream_multihop(edges=edges_ad, sources=skids, hops=1, exclude=exclude+all_order4+all_order3+all_order2, exclude_unpaired=True, pairs=pairs)[0] for skids in input_types.order4]
input_types['order5'] = order5

all_order5 = list(np.unique([x for sublist in input_types.order5 for x in sublist]))
order6 = [Promat.downstream_multihop(edges=edges_ad, sources=skids, hops=1, exclude=exclude+all_order5+all_order4+all_order3+all_order2, exclude_unpaired=True, pairs=pairs)[0] for skids in input_types.order5]
input_types['order6'] = order6

all_order6 = list(np.unique([x for sublist in input_types.order6 for x in sublist]))

# %%
# export IDs for modality layers

[pymaid.add_annotations(input_types.order2.loc[index], f'mw {input_types.type.loc[index]} 2nd_order') for index in input_types.index]
pymaid.add_meta_annotations([f'mw {input_types.type.loc[index]} 2nd_order' for index in input_types.index], 'mw brain inputs 2nd_order')
[pymaid.add_annotations(input_types.order3.loc[index], f'mw {input_types.type.loc[index]} 3rd_order') for index in input_types.index]
pymaid.add_meta_annotations([f'mw {input_types.type.loc[index]} 3rd_order' for index in input_types.index], 'mw brain inputs 3rd_order')
[pymaid.add_annotations(input_types.order4.loc[index], f'mw {input_types.type.loc[index]} 4th_order') for index in input_types.index]
pymaid.add_meta_annotations([f'mw {input_types.type.loc[index]} 4th_order' for index in input_types.index], 'mw brain inputs 4th_order')
[pymaid.add_annotations(input_types.order5.loc[index], f'mw {input_types.type.loc[index]} 5th_order') for index in input_types.index if len(input_types.order5.loc[index])!=0]
pymaid.add_meta_annotations([f'mw {input_types.type.loc[index]} 5th_order' for index in input_types.index if len(input_types.order5.loc[index])!=0], 'mw brain inputs 5th_order')

input_types = input_types.set_index('type') # for future chunks

# %%
# identify LNs in each layer

pymaid.clear_cache()

# load all-all (summed) adjacency matrix and axo-axonic adjacency matrix
subgraph = ['mw brain paper clustered neurons', 'mw brain accessory neurons']
summed_adj = Promat.pull_adj(type_adj='all-all', date=data_date, subgraph=subgraph)
aa_adj = Promat.pull_adj(type_adj='aa', date=data_date, subgraph=subgraph)

order = ['olfactory', 'gustatory-external', 'gustatory-pharyngeal', 'enteric', 'thermo-warm', 'thermo-cold', 'visual', 'noci', 'mechano-Ch', 'mechano-II/III', 'proprio', 'respiratory']

sens = [Celltype_Analyzer.get_skids_from_meta_annotation(f'mw {celltype}') for celltype in order]
order2_ct = [Celltype(f'2nd-order {celltype}', pymaid.get_skids_by_annotation(f'mw {celltype} 2nd_order')) for celltype in order]
order3_ct = [Celltype(f'3rd-order {celltype}', pymaid.get_skids_by_annotation(f'mw {celltype} 3rd_order')) for celltype in order]
order4_ct = [Celltype(f'4th-order {celltype}', pymaid.get_skids_by_annotation(f'mw {celltype} 4th_order')) for celltype in order]
exclude = [pymaid.get_skids_by_annotation(x) for x in ['mw MBON', 'mw MBIN', 'mw RGN', 'mw dVNC', 'mw dSEZ', 'mw KC', 'mw motor']] # 'mw motor' refers to motorneurons in SEZ
exclude = [x for sublist in exclude for x in sublist]
exclude = exclude + brain_inputs

# identify LNs
# use 0.5 output fraction within group threshold
threshold = 0.5
LNs_2nd = [celltype.identify_LNs(threshold, summed_adj, aa_adj, sens[i], outputs, exclude=exclude, pairs_path=pairs_path)[0] for i, celltype in enumerate(order2_ct)]
LNs_3rd = [celltype.identify_LNs(threshold, summed_adj, aa_adj, order2_ct[i].get_skids(), outputs, exclude=exclude, pairs_path=pairs_path)[0] for i, celltype in enumerate(order3_ct)]
LNs_4th = [celltype.identify_LNs(threshold, summed_adj, aa_adj, order3_ct[i].get_skids(), outputs, exclude=exclude, pairs_path=pairs_path)[0] for i, celltype in enumerate(order4_ct)]

# export LNs
[pymaid.add_annotations(LNs_2nd[i], f'mw brain 2nd_order LN {name}') for i, name in enumerate(order) if len(LNs_2nd[i])>0]
pymaid.add_meta_annotations([f'mw brain 2nd_order LN {name}' for i, name in enumerate(order) if len(LNs_2nd[i])>0], 'mw brain inputs 2nd_order LN')
[pymaid.add_annotations(LNs_3rd[i], f'mw brain 3rd_order LN {name}') for i, name in enumerate(order) if len(LNs_3rd[i])>0]
pymaid.add_meta_annotations([f'mw brain 3rd_order LN {name}' for i, name in enumerate(order) if len(LNs_3rd[i])>0], 'mw brain inputs 3rd_order LN')
[pymaid.add_annotations(LNs_4th[i], f'mw brain 4th_order LN {name}') for i, name in enumerate(order) if len(LNs_4th[i])>0]
pymaid.add_meta_annotations([f'mw brain 4th_order LN {name}' for i, name in enumerate(order) if len(LNs_4th[i])>0], 'mw brain inputs 4th_order LN')

# add special case for olfactory/gustatory 2nd-order because it's so interconnected
pymaid.clear_cache()
ct_skids = pymaid.get_skids_by_annotation('mw olfactory 2nd_order') + pymaid.get_skids_by_annotation('mw gustatory-external 2nd_order') + pymaid.get_skids_by_annotation('mw gustatory-pharyngeal 2nd_order')
input_skids = Celltype_Analyzer.get_skids_from_meta_annotation('mw olfactory') + Celltype_Analyzer.get_skids_from_meta_annotation('mw gustatory-external') + Celltype_Analyzer.get_skids_from_meta_annotation('mw gustatory-pharyngeal')
olf_gust_order2 = Celltype('2nd-order olfactory-gustatory', ct_skids)
olf_gust_LN_order2 = olf_gust_order2.identify_LNs(threshold, summed_adj, aa_adj, input_skids, outputs, exclude=exclude, pairs_path=pairs_path)[0] 
pymaid.add_annotations(olf_gust_LN_order2, 'mw brain 2nd_order LN olfactory-gustatory')
pymaid.clear_cache()
pymaid.add_meta_annotations('mw brain 2nd_order LN olfactory-gustatory', 'mw brain inputs 2nd_order LN')
# APL not picked up because it doesn't receive input from 2nd-order olfactory neurons

# add input/output type of LN
# use 0.5 output fraction within group threshold
threshold = 0.5
LNs_io_2nd = [celltype.identify_in_out_LNs(threshold, summed_adj, outputs, inputs, exclude=exclude, pairs_path=pairs_path)[0] for i, celltype in enumerate(order2_ct)]
LNs_io_3rd = [celltype.identify_in_out_LNs(threshold, summed_adj, outputs, inputs, exclude=exclude, pairs_path=pairs_path)[0] for i, celltype in enumerate(order3_ct)]
LNs_io_4th = [celltype.identify_in_out_LNs(threshold, summed_adj, outputs, inputs, exclude=exclude, pairs_path=pairs_path)[0] for i, celltype in enumerate(order4_ct)]

# export LNs
[pymaid.add_annotations(LNs_io_2nd[i], f'mw brain 2nd_order LN_io {name}') for i, name in enumerate(order) if len(LNs_io_2nd[i])>0]
pymaid.add_meta_annotations([f'mw brain 2nd_order LN_io {name}' for i, name in enumerate(order) if len(LNs_io_2nd[i])>0], 'mw brain inputs 2nd_order LN_io')
[pymaid.add_annotations(LNs_io_3rd[i], f'mw brain 3rd_order LN_io {name}') for i, name in enumerate(order) if len(LNs_io_3rd[i])>0]
pymaid.add_meta_annotations([f'mw brain 3rd_order LN_io {name}' for i, name in enumerate(order) if len(LNs_io_3rd[i])>0], 'mw brain inputs 3rd_order LN_io')
[pymaid.add_annotations(LNs_io_4th[i], f'mw brain 4th_order LN_io {name}') for i, name in enumerate(order) if len(LNs_io_4th[i])>0]
pymaid.add_meta_annotations([f'mw brain 4th_order LN_io {name}' for i, name in enumerate(order) if len(LNs_io_4th[i])>0], 'mw brain inputs 4th_order LN_io')

pymaid.clear_cache()
LNs_o = [Celltype_Analyzer.get_skids_from_meta_annotation(f'mw brain inputs {order} LN') for order in ['2nd_order', '3rd_order', '4th_order']]
LNs_io = [Celltype_Analyzer.get_skids_from_meta_annotation(f'mw brain inputs {order} LN_io') for order in ['2nd_order', '3rd_order']]
LNs_o = [x for sublist in LNs_o for x in sublist]
LNs_io = [x for sublist in LNs_io for x in sublist]

pymaid.add_annotations(LNs_o, 'mw LNs_cohort')
pymaid.add_annotations(LNs_io, 'mw LNs_noncohort')

# %%
# identify PNs by exclusion (not LN and not output neuron[dVNC, dSEZ, RGN])

# pull skids for 2nd-order modalities
order = ['olfactory', 'gustatory-external', 'gustatory-pharyngeal', 'enteric', 'thermo-warm', 'thermo-cold', 'visual', 'noci', 'mechano-Ch', 'mechano-II/III', 'proprio', 'respiratory']
order2_ct = [Celltype(f'{celltype}', pymaid.get_skids_by_annotation(f'mw {celltype} 2nd_order')) for celltype in order]

LNs = Celltype_Analyzer.get_skids_from_meta_annotation('mw brain LNs')
accessory_nonPN = pymaid.get_skids_by_annotation('mw brain accessory non-PN') # SEZ neurons that also mostly output to SEZ; so not really dSEZ but also not brain PN
outputs = pymaid.get_skids_by_annotation(['mw dVNC', 'mw dSEZ', 'mw RGN'])

PNs_2nd = [Celltype(f'mw {celltype.name} 2nd_order PN', np.setdiff1d(celltype.skids, LNs + accessory_nonPN + outputs)) for celltype in order2_ct]

[pymaid.add_annotations(celltype.skids, celltype.name) for celltype in PNs_2nd]
pymaid.add_meta_annotations([celltype.name for celltype in PNs_2nd], 'mw brain inputs 2nd_order PN')

# %%
