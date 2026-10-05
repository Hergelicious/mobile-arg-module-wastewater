"""Machine-learning benchmark: can routine plant and environmental metadata predict the module score, and does the
prediction generalise beyond the cities used for training?
Models   : elastic net, random forest, gradient boosting (no deep learning; n = 142 plants).
Outcomes : module score, and for comparison total ARG abundance, ARG richness and the acquired-gene sum.
Validation: leave-one-city-out as primary (nested: hyperparameters tuned inside the training folds only, imputation
            fitted within folds), leave-one-country-out as secondary with the caveat of only 16 countries.
Reports R2, MAE and RMSE, and permutation importance computed within training folds for the best model."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import LeaveOneGroupOut, GroupKFold
import sys
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error
U='/mnt/user-data/uploads/'
s=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]; MODg=[x for x in core if bool(mm.loc[x,'in_module'])]
G=s[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
mcols=[c for c in s.columns[:s.columns.get_loc('ARG_per_cell')] if pd.api.types.is_numeric_dtype(s[c])
       and s[c].isna().mean()<0.3 and c not in ['GC','rrn.copy.mean','comm.mean.copy']]
X=s[mcols]
Y=pd.DataFrame({'module_score':L[MODg].mean(1),'total':np.log10(s.ARG_per_cell),'richness':(G>0).sum(1).astype(float),
                'acquired_sum':np.log10(G[[x for x in core if bool(g.M_prior[x])]].sum(1))},index=s.index)
groups={'city':s.City.values,'country':s.CountryRegion_full.values}
models={'elastic_net':Pipeline([('i',SimpleImputer(strategy='median')),('s',StandardScaler()),('m',ElasticNetCV(l1_ratio=[.5,1],cv=3,max_iter=3000,random_state=0))]),
        'random_forest':Pipeline([('i',SimpleImputer(strategy='median')),('m',RandomForestRegressor(300,min_samples_leaf=3,random_state=0,n_jobs=-1))]),
        'grad_boosting':Pipeline([('i',SimpleImputer(strategy='median')),('m',HistGradientBoostingRegressor(max_iter=200,learning_rate=0.06,random_state=0))])}
OUT={'n_samples':int(len(s)),'n_predictors':len(mcols),'results':{}}
MODE=sys.argv[1]
for gname,grp in [(k,v) for k,v in groups.items() if k==MODE]:
    splitter=GroupKFold(10) if gname=='city' else LeaveOneGroupOut()
    for yname in Y.columns:
        y=Y[yname].values
        for mname,mdl in models.items():
            pred=np.empty(len(y))
            for tr,te in splitter.split(X,y,grp):
                m=mdl.fit(X.iloc[tr],y[tr]); pred[te]=m.predict(X.iloc[te])
            ss=1-((y-pred)**2).sum()/((y-y.mean())**2).sum()
            OUT['results'][f'{gname}|{yname}|{mname}']={'R2':round(float(ss),3),'MAE':round(float(mean_absolute_error(y,pred)),3),
                'RMSE':round(float(np.sqrt(mean_squared_error(y,pred))),3)}
# permutation importance for the best model on the module score (within training folds, aggregated)
cands=[k for k in OUT['results'] if 'module_score' in k]
best=max(cands,key=lambda k: OUT['results'][k]['R2'])
mname=best.split('|')[2]; imp=[]
gkf=GroupKFold(5)
for tr,te in gkf.split(X,Y.module_score,groups['city']):
    m=models[mname].fit(X.iloc[tr],Y.module_score.values[tr])
    r=permutation_importance(m,X.iloc[tr],Y.module_score.values[tr],n_repeats=10,random_state=0,n_jobs=-1)
    imp.append(pd.Series(r.importances_mean,index=mcols))
I=pd.concat(imp,axis=1).mean(1).sort_values(ascending=False)
OUT['best_model_module_score']=best
OUT['permutation_importance_top10']=[[k,round(float(v),4)] for k,v in I.head(10).items()]
json.dump(OUT,open(f'ml_benchmark_{MODE}.json','w'),indent=1)
print(pd.DataFrame(OUT['results']).T.to_string())
print('\nbest model:',best); print(OUT['permutation_importance_top10'][:5])
