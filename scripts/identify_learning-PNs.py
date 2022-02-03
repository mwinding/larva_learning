# %%

import pymaid
import contools
from contools import generate_adjs

from pymaid_creds import url, name, password, token
rm = pymaid.CatmaidInstance(url, token, name, password)
pairs_path = 'data/pairs/pairs-2021-04-06.csv'

pairs = contools.Promat.get_pairs(pairs_path)
# %%
