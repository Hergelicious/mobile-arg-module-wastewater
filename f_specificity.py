"""EXPLORATORY (added after test F): is the consumption association specific to the module?"""
import warnings; warnings.filterwarnings("ignore")
import json, re, numpy as np, pandas as pd
from scipy.stats import spearmanr
U='/mnt/user-data/uploads/'; S=json.load(open('frozen_spec3.json'))
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq=set(sd2[sd2.database=='resfinder'].gene); fg=set(sd2[sd2.database!='resfinder'].gene)
MOD=sorted({p for v in S['gene_mapping']['module_mapped'].values() for p in v if p in m.columns})
prev=(m>0).mean(); core=set(prev[prev>=0.5].index)
sets={'module':[g for g in MOD if g in core],'acquired_not_module':sorted((acq&core)-set(MOD)),'functional_latent':sorted((fg-acq)&core)}
X=m.values.copy()
for i in range(X.shape[0]):
    row=X[i]; pos=row[row>0]
    if len(pos): row[row==0]=0.65*pos.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m.index,columns=m.columns)
F=pd.read_csv('F_country_table.csv')
out={}
for k,gs in sets.items():
    sc=LD[gs].mean(1).groupby(meta.country.values).mean()
    v=F.country.map(sc); r,p=spearmanr(v,F.ddd); out[k]={'n_genes':len(gs),'rho':float(r),'P':float(p)}
json.dump(out,open('f_specificity.json','w'),indent=1); print(json.dumps(out,indent=1))
