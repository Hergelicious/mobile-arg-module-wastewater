"""Per-gene within-country Spearman correlation with intI1 (same computation as run_frozen.py test B) -> intI1_rho.csv"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
from scipy.stats import rankdata
s=pd.read_excel('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
mge=pd.read_excel('/mnt/user-data/uploads/41467_2025_59019_MOESM5_ESM.xlsx',header=2)
rows=mge[mge['MGE_Gene_name'].astype(str).str.lower()=='inti1']
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
intI=rows[[c for c in rows.columns if c in s.index]].astype(float).sum(0).reindex(s.index)
G=s[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2); li=np.log10(intI+intI[intI>0].min()/2)
LP=L.groupby(s.WWTPID).mean(); IP=li.groupby(s.WWTPID).mean().loc[LP.index]
ctry=s.groupby('WWTPID').CountryRegion_full.first().loc[LP.index]
Ld=LP-LP.groupby(ctry.values).transform('mean'); Id=IP-IP.groupby(ctry.values).transform('mean')
RL=np.apply_along_axis(rankdata,0,Ld.values); ri=rankdata(Id.values); ri=(ri-ri.mean())/ri.std()
r=((RL-RL.mean(0))/RL.std(0)*ri[:,None]).mean(0)
out=pd.Series(r,index=core,name='rho_intI1_within_country')
FR=json.load(open('frozen_results.json'))
assert abs(out[g.M_prior].mean()-FR['B_primary']['mean_rho_mobile'])<1e-12, 'mismatch with frozen run'
out.to_csv('intI1_rho.csv'); print('saved; matches frozen run. intI1 ORFs summed:',len(rows),'| MGE names:',sorted(rows['MGE_Gene_name'].unique()))
