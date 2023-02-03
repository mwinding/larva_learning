# %%
from pymaid_creds import url, name, password, token
from data_settings import data_date, pairs_path
import pymaid

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx
from contools import Promat, Prograph

rm = pymaid.CatmaidInstance(url, token, name, password)
FBN = pymaid.get_skids_by_annotation('mw FBN')
DAN = pymaid.get_skids_by_annotation('mw MBIN subclass_DAN')
MBON = pymaid.get_skids_by_annotation('mw MBON')

ad_edges_frac = Promat.pull_edges(type_edges='ad', threshold=0.01, data_date=data_date, pairs_combined=False)
ad_edges_frac = ad_edges_frac.set_index('downstream_skid', drop=False)

adj_ad = Promat.pull_adj(type_adj='ad', data_date=data_date)
inputs = pd.read_csv(f'data/adj/inputs_{data_date}.csv', index_col=0)
pairs = Promat.get_pairs(pairs_path=pairs_path)

# %%
# upstream of DAN

DAN_edges = ad_edges_frac.loc[DAN]
DAN_edges = DAN_edges.drop(['upstream_side', 'downstream_side', 'type', 'upstream_status', 'downstream_status'], axis=1)

all_syn = []
for i in range(len(DAN_edges.index)):
    upstream_skid = DAN_edges.iloc[i, :].upstream_skid
    downstream_skid = DAN_edges.iloc[i, :].downstream_skid
    syn = adj_ad.loc[upstream_skid, downstream_skid]
    all_syn.append(int(syn))

DAN_edges['edge_weight_syn'] = all_syn
DAN_edges.index = range(len(DAN_edges.index))

# change formatting and save
DAN_edges.columns = ['upstream_skid', 'downstream_skid', 'edge_weight_fraction', 'edge_weight']
DAN_edges = DAN_edges.loc[:, ['upstream_skid', 'downstream_skid', 'edge_weight', 'edge_weight_fraction']]
DAN_edges.to_csv('plots/all-to-DAN_edges.csv')

# %%
# from MBON to FBN

MBON_FBN_edges = ad_edges_frac.loc[FBN]
MBON_FBN_edges.index = MBON_FBN_edges.upstream_skid
MBON_FBN_edges = MBON_FBN_edges.loc[np.intersect1d(MBON, MBON_FBN_edges.index)]
MBON_FBN_edges = MBON_FBN_edges.drop(['upstream_side', 'downstream_side', 'type', 'upstream_status', 'downstream_status'], axis=1)

all_syn = []
for i in range(len(MBON_FBN_edges.index)):
    upstream_skid = MBON_FBN_edges.iloc[i, :].upstream_skid
    downstream_skid = MBON_FBN_edges.iloc[i, :].downstream_skid
    syn = adj_ad.loc[upstream_skid, downstream_skid]
    all_syn.append(int(syn))

MBON_FBN_edges['edge_weight_syn'] = all_syn
MBON_FBN_edges.index = range(len(MBON_FBN_edges.index))

# change formatting and save
MBON_FBN_edges.columns = ['upstream_skid', 'downstream_skid', 'edge_weight_fraction', 'edge_weight']
MBON_FBN_edges = MBON_FBN_edges.loc[:, ['upstream_skid', 'downstream_skid', 'edge_weight', 'edge_weight_fraction']]

MBON_FBN_edges.to_csv('plots/MBON-to-FBN_edges.csv')
# %%
# all edges to DANs, with both hemispheres combined

ad_edges_frac = Promat.pull_edges(type_edges='ad', threshold=0.01, data_date=data_date, pairs_combined=True)
ad_edges_frac = ad_edges_frac.set_index('downstream_pair_id', drop=False)

DAN_edges = ad_edges_frac.loc[np.intersect1d(DAN, ad_edges_frac.index)]

left_syns=[]
right_syns=[]
left_totals=[]
right_totals=[]
for i in range(len(DAN_edges)):
    row = DAN_edges.iloc[i, :]
    if(row.type=='ipsilateral'):
        left_syn = adj_ad.loc[row.upstream_pair_id, row.downstream_pair_id]
        right_syn = adj_ad.loc[Promat.get_paired_skids(row.upstream_pair_id, pairs)[1], Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]]
        left_total = inputs.loc[row.downstream_pair_id].dendrite_input
        right_total = inputs.loc[Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]].dendrite_input

    if(row.type=='contralateral'):
        left_syn = adj_ad.loc[Promat.get_paired_skids(row.upstream_pair_id, pairs)[1], row.downstream_pair_id]
        right_syn = adj_ad.loc[row.upstream_pair_id, Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]]
        left_total = inputs.loc[row.downstream_pair_id].dendrite_input
        right_total = inputs.loc[Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]].dendrite_input

    left_syns.append(int(left_syn))
    right_syns.append(int(right_syn))
    left_totals.append(int(left_total))
    right_totals.append(int(right_total))

DAN_edges['left_syn'] = left_syns
DAN_edges['right_syn'] = right_syns
DAN_edges['left_total'] = left_totals
DAN_edges['right_total'] = right_totals

DAN_edges = DAN_edges.loc[:, ['upstream_pair_id', 'downstream_pair_id', 'left_syn', 'right_syn', 'left_total', 'right_total', 'left', 'right']]
DAN_edges.columns = ['upstream_pair_id', 'downstream_pair_id', 'left_input', 'right_input', 'left_total_input', 'right_total_input', 'left_fraction', 'right_fraction']

DAN_edges.index = pd.MultiIndex.from_frame(DAN_edges.loc[:, ['upstream_pair_id', 'downstream_pair_id']])
rows = []
for skid in np.unique(DAN_edges.downstream_pair_id):
    upstream_pair_id = -1
    downstream_pair_id = skid
    left_total_input = DAN_edges.loc[(slice(None), skid), :].left_total_input.iloc[0]
    right_total_input = DAN_edges.loc[(slice(None), skid), :].right_total_input.iloc[0]
    left_input = left_total_input - sum(DAN_edges.loc[(slice(None), skid), :].left_input)
    right_input = right_total_input - sum(DAN_edges.loc[(slice(None), skid), :].left_input)
    left_fraction = left_input/left_total_input
    right_fraction = right_input/right_total_input
    rows.append([upstream_pair_id, downstream_pair_id, int(left_input), int(right_input), int(left_total_input), int(right_total_input), left_fraction, right_fraction])

DAN_edges = pd.concat([DAN_edges, pd.DataFrame(rows, columns=DAN_edges.columns)], axis=0)
DAN_edges.index = range(len(DAN_edges.index))
DAN_edges = DAN_edges.sort_values(by=['downstream_pair_id', 'upstream_pair_id'])
DAN_edges.index = range(len(DAN_edges.index))

DAN_edges.to_csv('plots/all-to-DAN_edges.csv')

# %%
# all edges to FBNs, with both hemispheres combined


ad_edges_frac = Promat.pull_edges(type_edges='ad', threshold=0.01, data_date=data_date, pairs_combined=True)
ad_edges_frac = ad_edges_frac.set_index('downstream_pair_id', drop=False)

FBN_edges = ad_edges_frac.loc[np.intersect1d(FBN, ad_edges_frac.index)]

left_syns=[]
right_syns=[]
left_totals=[]
right_totals=[]
for i in range(len(FBN_edges)):
    row = FBN_edges.iloc[i, :]
    # some nonpaired neurons were upstream of FBNs, excluded those
    if((row.type=='ipsilateral') & (row.upstream_status=='paired') & (row.downstream_status=='paired')):
        left_syn = adj_ad.loc[row.upstream_pair_id, row.downstream_pair_id]
        right_syn = adj_ad.loc[Promat.get_paired_skids(row.upstream_pair_id, pairs)[1], Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]]
        left_total = inputs.loc[row.downstream_pair_id].dendrite_input
        right_total = inputs.loc[Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]].dendrite_input

    if((row.type=='contralateral') & (row.upstream_status=='paired') & (row.downstream_status=='paired')):
        left_syn = adj_ad.loc[Promat.get_paired_skids(row.upstream_pair_id, pairs)[1], row.downstream_pair_id]
        right_syn = adj_ad.loc[row.upstream_pair_id, Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]]
        left_total = inputs.loc[row.downstream_pair_id].dendrite_input
        right_total = inputs.loc[Promat.get_paired_skids(row.downstream_pair_id, pairs)[1]].dendrite_input

    left_syns.append(int(left_syn))
    right_syns.append(int(right_syn))
    left_totals.append(int(left_total))
    right_totals.append(int(right_total))

FBN_edges['left_syn'] = left_syns
FBN_edges['right_syn'] = right_syns
FBN_edges['left_total'] = left_totals
FBN_edges['right_total'] = right_totals

FBN_edges = FBN_edges.loc[:, ['upstream_pair_id', 'downstream_pair_id', 'left_syn', 'right_syn', 'left_total', 'right_total', 'left', 'right']]
FBN_edges.columns = ['upstream_pair_id', 'downstream_pair_id', 'left_input', 'right_input', 'left_total_input', 'right_total_input', 'left_fraction', 'right_fraction']

FBN_edges.index = pd.MultiIndex.from_frame(FBN_edges.loc[:, ['upstream_pair_id', 'downstream_pair_id']])
rows = []
for skid in np.unique(FBN_edges.downstream_pair_id):
    upstream_pair_id = -1
    downstream_pair_id = skid
    left_total_input = FBN_edges.loc[(slice(None), skid), :].left_total_input.iloc[0]
    right_total_input = FBN_edges.loc[(slice(None), skid), :].right_total_input.iloc[0]
    left_input = left_total_input - sum(FBN_edges.loc[(slice(None), skid), :].left_input)
    right_input = right_total_input - sum(FBN_edges.loc[(slice(None), skid), :].left_input)
    left_fraction = left_input/left_total_input
    right_fraction = right_input/right_total_input
    rows.append([upstream_pair_id, downstream_pair_id, int(left_input), int(right_input), int(left_total_input), int(right_total_input), left_fraction, right_fraction])

FBN_edges = pd.concat([FBN_edges, pd.DataFrame(rows, columns=FBN_edges.columns)], axis=0)
FBN_edges.index = range(len(FBN_edges.index))
FBN_edges = FBN_edges.sort_values(by=['downstream_pair_id', 'upstream_pair_id'])
FBN_edges.index = range(len(FBN_edges.index))

FBN_edges.to_csv('plots/all-to-FBN_edges.csv')

# %%
