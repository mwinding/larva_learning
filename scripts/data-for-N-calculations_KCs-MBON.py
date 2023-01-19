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

brain = pymaid.get_skids_by_annotation('mw brain neurons')
pairs = Promat.get_pairs(pairs_path=pairs_path)

inputs = pd.read_csv(f'data/adj/inputs_{data_date}.csv', index_col=0)


# %%
# 

KC = pymaid.get_skids_by_annotation('mw KC')
ad_adj = Promat.pull_adj('ad', data_date=data_date)

MBONi1 = Promat.load_pairs_from_annotation('MBON-i1', pairList=pairs)
MBONj1 = Promat.load_pairs_from_annotation('MBON-j1', pairList=pairs)
MBONi1.index = ['MBON-i1']
MBONj1.index = ['MBON-j1']

MBON_data = pd.concat([MBONi1, MBONj1])

new_rows = []

for i in range(len(MBON_data.index)):
    leftid = MBON_data.loc[MBON_data.index[i], 'leftid']
    rightid = MBON_data.loc[MBON_data.index[i], 'rightid']

    left_input = int(inputs.loc[leftid, 'dendrite_input'])
    right_input = int(inputs.loc[rightid, 'dendrite_input'])

    left_KC_input_ad = int(ad_adj.loc[np.intersect1d(KC, ad_adj.index), leftid].sum())
    right_KC_input_ad = int(ad_adj.loc[np.intersect1d(KC, ad_adj.index), rightid].sum())

    left_fraction = left_KC_input_ad/left_input
    right_fraction = right_KC_input_ad/right_input

    new_rows.append([left_KC_input_ad, right_KC_input_ad, left_input, right_input, left_fraction, right_fraction])

df = pd.DataFrame(new_rows, index=['MBON-i1', 'MBON-j1'])
df.columns = ['left_KC_input', 'right_KC_input', 'left_total_input', 'right_total_input', 'left_fraction', 'right_fraction']

df.to_csv('plots/synapses-from-KCs-to-MBONs.csv')
# %%
