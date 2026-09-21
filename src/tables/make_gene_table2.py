"""Adds abundance, variability and detection columns to gene_table.csv -> gene_table2.csv."""
import numpy as np, pandas as pd
s=pd.read_excel('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table.csv',index_col=0); core=list(g.index)
G=s[core]; L=np.log10(G+G.replace(0,np.nan).min()/2)
g['mean_log']=np.log10(G.mean()); g['sd_log']=L.groupby(s.WWTPID).mean().std(); g['zero']=(G==0).mean()
g.to_csv('gene_table2.csv')
