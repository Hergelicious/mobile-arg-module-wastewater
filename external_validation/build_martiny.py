"""Build the Martiny sewage abundance matrix for the frozen Tier A tests."""
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np, json
U='/mnt/user-data/uploads/'
use=['gene','complete_name','genepid','fragmentCountAln_adj_kb','fragmentCountAln_adj']
chunks=[]
for ch in pd.read_csv(U+'panres_counts.csv',usecols=use,chunksize=400_000):
    chunks.append(ch)
d=pd.concat(chunks,ignore_index=True); print('rows',len(d),'| samples',d.complete_name.nunique(),'| genes',d.gene.nunique())
# one row per sample x gene (sum in case of duplicates), length-adjusted fragments
m=d.groupby(['complete_name','gene'],observed=True)['fragmentCountAln_adj_kb'].sum().unstack(fill_value=0.0)
print('matrix',m.shape)
m.to_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx')
print('metadata rows',len(meta),'| matched to matrix:',len(set(m.index)&set(meta.complete_name)))
