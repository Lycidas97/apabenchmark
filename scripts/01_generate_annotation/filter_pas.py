#!/usr/bin/env python3
"""Script export of filter_pas.ipynb."""


# %%
import pandas as pd
import os

os.makedirs("../../data/int_data/annotations", exist_ok=True)

# %%
mouse_pas = pd.read_csv("../../data/int_data/annotations/mouse_integrated_pas.bed", sep='\t', header=None).drop_duplicates(subset=[3])
human_pas = pd.read_csv("../../data/int_data/annotations/human_integrated_pas.bed", sep='\t', header=None).drop_duplicates(subset=[3])

# %%
mouse_pas_forward = mouse_pas[mouse_pas[5] == '+'].copy().sort_values(by=[0, 1, 2])
mouse_pas_reverse = mouse_pas[mouse_pas[5] == '-'].copy().sort_values(by=[0, 1, 2])

mouse_pas_forward["gap_downstream"] = mouse_pas_forward[2].diff()
mouse_pas_reverse["gap_downstream"] = mouse_pas_reverse[1].diff()
mouse_pas_forward["gap_upstream"] = -mouse_pas_forward[2].diff(periods=-1)
mouse_pas_reverse["gap_upstream"] = -mouse_pas_reverse[1].diff(periods=-1)

mouse_pas_merged = pd.concat([mouse_pas_forward, mouse_pas_reverse]).sort_values(by=[0, 1, 2])
mouse_pas_merged["gap_downstream"] = mouse_pas_merged["gap_downstream"].fillna(10000)
mouse_pas_merged["gap_upstream"] = mouse_pas_merged["gap_upstream"].fillna(10000)

# %%
human_pas_forward = human_pas[human_pas[5] == '+'].copy().sort_values(by=[0, 1, 2])
human_pas_reverse = human_pas[human_pas[5] == '-'].copy().sort_values(by=[0, 1, 2])

human_pas_forward["gap_downstream"] = human_pas_forward[2].diff()
human_pas_reverse["gap_downstream"] = human_pas_reverse[1].diff()
human_pas_forward["gap_upstream"] = -human_pas_forward[2].diff(periods=-1)
human_pas_reverse["gap_upstream"] = -human_pas_reverse[1].diff(periods=-1)

human_pas_merged = pd.concat([human_pas_forward, human_pas_reverse]).sort_values(by=[0, 1, 2])
human_pas_merged["gap_downstream"] = human_pas_merged["gap_downstream"].fillna(10000)
human_pas_merged["gap_upstream"] = human_pas_merged["gap_upstream"].fillna(10000)

# %%
mouse_pas_500 = mouse_pas_merged[(mouse_pas_merged["gap_upstream"] > 500) & (mouse_pas_merged["gap_downstream"] > 500)].copy()
mouse_pas_500 = mouse_pas_500.iloc[:, 0:10]

mouse_pas_1000 = mouse_pas_merged[(mouse_pas_merged["gap_upstream"] > 1000) & (mouse_pas_merged["gap_downstream"] > 1000)].copy()
mouse_pas_1000 = mouse_pas_1000.iloc[:, 0:10]

pas_num = mouse_pas.groupby(6).count()[0].reset_index()
single_pas_gene = pas_num[pas_num[0] == 1][6].tolist()
mouse_pas_single = mouse_pas[mouse_pas[6].isin(single_pas_gene)].copy()

mouse_pas_500.to_csv('../../data/int_data/annotations/mouse_pas_500.bed', sep='\t', header=None, index=None)
mouse_pas_1000.to_csv('../../data/int_data/annotations/mouse_pas_1000.bed', sep='\t', header=None, index=None)
mouse_pas_single.to_csv('../../data/int_data/annotations/mouse_pas_single.bed', sep='\t', header=None, index=None)

# %%
human_pas_500 = human_pas_merged[(human_pas_merged["gap_upstream"] > 500) & (human_pas_merged["gap_downstream"] > 500)].copy()
human_pas_500 = human_pas_500.iloc[:, 0:10]

human_pas_1000 = human_pas_merged[(human_pas_merged["gap_upstream"] > 1000) & (human_pas_merged["gap_downstream"] > 1000)].copy()
human_pas_1000 = human_pas_1000.iloc[:, 0:10]

pas_num = human_pas.groupby(6).count()[0].reset_index()
single_pas_gene = pas_num[pas_num[0] == 1][6].tolist()
human_pas_single = human_pas[human_pas[6].isin(single_pas_gene)].copy()

human_pas_500.to_csv('../../data/int_data/annotations/human_pas_500.bed', sep='\t', header=None, index=None)
human_pas_1000.to_csv('../../data/int_data/annotations/human_pas_1000.bed', sep='\t', header=None, index=None)
human_pas_single.to_csv('../../data/int_data/annotations/human_pas_single.bed', sep='\t', header=None, index=None)

# %%
print(f"mouse_pas_500: {len(mouse_pas_500)}")
print(f"mouse_pas_1000: {len(mouse_pas_1000)}")
print(f"mouse_pas_single: {len(mouse_pas_single)}")
print(f"mouse_pas: {len(mouse_pas)}")
