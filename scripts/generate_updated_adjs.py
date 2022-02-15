# %%

import pymaid
import contools
from contools import generate_adjs
import numpy as np

from pymaid_creds import url, name, password, token
rm = pymaid.CatmaidInstance(url, token, name, password)

# %%

all_neurons = pymaid.get_skids_by_annotation(['mw brain paper clustered neurons', 'mw brain accessory neurons'])
remove_neurons = pymaid.get_skids_by_annotation('mw brain very incomplete')
all_neurons = list(np.setdiff1d(all_neurons, remove_neurons)) # remove neurons that are so incomplete, they have no split point

split_tag = 'mw axon split'
special_split_tags = ['mw axon start', 'mw axon end']
not_split_skids = pymaid.get_skids_by_annotation(['mw unsplittable'])

generate_adjs.adj_split_axons_dendrites(all_neurons, split_tag, special_split_tags, not_split_skids)

# %%

# generate edge list with average pairwise threshold = 3
threshold = 3
pairs_path = 'data/pairs/pairs-2022-02-14.csv'
pairs = contools.Promat.get_pairs(pairs_path=pairs_path)
date = '2022-02-15'
generate_adjs.edge_thresholds(path='data/adj', threshold=threshold, left_annot='mw left', right_annot='mw right', pairs = pairs, fraction_input=False, date=date)

# %%

# generate edge list with %input threshold = 0.01
threshold = 0.01
pairs_path = 'data/pairs/pairs-2022-02-14.csv'
pairs = contools.Promat.get_pairs(pairs_path=pairs_path)
date = '2022-02-10'
generate_adjs.edge_thresholds(path='data/adj', threshold=threshold, left_annot='mw left', right_annot='mw right', pairs = pairs, fraction_input=True, date=date)

# %%
