# %%

import pymaid
import contools
from contools import generate_adjs

from pymaid_creds import url, name, password, token
rm = pymaid.CatmaidInstance(url, token, name, password)

# %%

all_neurons = pymaid.get_skids_by_annotation(['mw brain paper clustered neurons', 'mw brain accessory neurons'])
split_tag = 'mw axon split'
special_split_tags = ['mw axon start', 'mw axon end']
not_split_skids = pymaid.get_skids_by_annotation(['mw unsplittable', 'mw partially differentiated', 'mw brain incomplete'])

generate_adjs.adj_split_axons_dendrites(all_neurons, split_tag, special_split_tags, not_split_skids)

# %%

# generate edge list with minimum pairwise threshold = 1
pairs_path = 'data/pairs/pairs-2021-04-06.csv'
pairs = contools.Promat.get_pairs(pairs_path=pairs_path)
generate_adjs.edge_thresholds(path='data/adj', threshold=1, left_annot='mw left', right_annot='mw right', pairs = pairs, fraction_input=False, date='2022-02-03')

# %%

# generate edge list with average pairwise threshold = 3
pairs_path = 'data/pairs/pairs-2021-04-06.csv'
pairs = contools.Promat.get_pairs(pairs_path=pairs_path)
generate_adjs.edge_thresholds(path='data/adj', threshold=3, left_annot='mw left', right_annot='mw right', pairs = pairs, fraction_input=False, date='2022-02-03')

# %%

# generate edge list with %input threshold = 0.01
pairs_path = 'data/pairs/pairs-2021-04-06.csv'
pairs = contools.Promat.get_pairs(pairs_path=pairs_path)
generate_adjs.edge_thresholds(path='data/adj', threshold=0.01, left_annot='mw left', right_annot='mw right', pairs = pairs, fraction_input=True, date='2022-02-03')

# %%
