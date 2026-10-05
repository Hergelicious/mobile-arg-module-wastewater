"""Transferability with imputation fitted inside each training fold (no leakage)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, cross_val_predict
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s=pd.read_excel(x,'Fig. S2b and Fig. S4',index_col=0)
L=pd.read_csv('sample_layers.csv',index_col=0).loc[s.index]
mcols=[c for c in s.columns[:s.columns.get_loc('ARG_per_cell')] if pd.api.types.is_numeric_dtype(s[c])
       and s[c].isna().mean()<0.3 and c not in ['GC','rrn.copy.mean','comm.mean.copy']]
X=s[mcols].values   # raw, with NaNs: imputation happens inside the pipeline
targets={'total':L.total,'intrinsic_layer':L.intrinsic_layer,'mobile_layer':L.mobile_layer,'mobile_fraction':L.mobile_fraction}
lev={'plant':s.WWTPID,'city':s.City,'country':s.CountryRegion_full}
out={}
for tn,y in targets.items():
    for ln,grp in lev.items():
        pipe=make_pipeline(SimpleImputer(strategy='median'),RandomForestRegressor(500,min_samples_leaf=3,random_state=0,n_jobs=-1))
        p=cross_val_predict(pipe,X,y.values,groups=grp.values,cv=GroupKFold(10))
        out[f'{tn}|{ln}']=float(1-((y.values-p)**2).sum()/((y.values-y.values.mean())**2).sum())
old=json.load(open('results.json'))['transfer_R2']
K=json.load(open('canonical.json')); K['transfer_R2']=out; K['n_predictors']=len(mcols)
json.dump(K,open('canonical.json','w'),indent=1)
for k in out: print(f'{k:28s} leaky {old[k]:+.3f}   in-fold {out[k]:+.3f}')
