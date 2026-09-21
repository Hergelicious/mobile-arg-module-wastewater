"""One consistent decoding analysis: paired draws, same folds, permutation null."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
rng=np.random.default_rng(31415)
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s=pd.read_excel(x,'Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table2.csv',index_col=0); core=list(g.index); M=g.M_prior.values
G=s[core].astype(float); LP=np.log10(G+G.replace(0,np.nan).min()/2).groupby(s.WWTPID).mean()
pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
multi=['Australia','Brazil','China','Germany','South Africa','Taiwan','USA']; k=pm.country.isin(multi)
X=LP[k]; y=pm.country[k].values; grp=pm.city[k].values
mob=list(g.index[M]); oth=list(g.index[~M])
def dec(genes):
    clf=make_pipeline(StandardScaler(),LogisticRegression(C=0.5,max_iter=3000,class_weight='balanced'))
    return balanced_accuracy_score(y,cross_val_predict(clf,X[list(genes)].values,y,groups=grp,cv=LeaveOneGroupOut()))
NDRAW=50
pa=np.array([dec(rng.choice(mob,30,replace=False)) for _ in range(NDRAW)])
oa=np.array([dec(rng.choice(oth,30,replace=False)) for _ in range(NDRAW)])
diff=pa.mean()-oa.mean()
nd=np.array([ (lambda p: dec(p[:30])-dec(p[30:60]))(rng.permutation(core)) for _ in range(300)])
cpl=g.M_coupling.values; hi=list(g.index[cpl>=np.quantile(cpl,2/3)]); lo=list(g.index[cpl<=np.quantile(cpl,1/3)])
rand=np.array([dec(rng.choice(core,30,replace=False)) for _ in range(200)])
D={'n_draws':NDRAW,'mobile':float(pa.mean()),'mobile_sd':float(pa.std()),'other':float(oa.mean()),'other_sd':float(oa.std()),
   'diff':float(diff),'p':float((np.sum(nd>=diff)+1)/301),'null_sd':float(nd.std()),
   'hi_tertile':float(dec(hi)),'lo_tertile':float(dec(lo)),'n_hi':len(hi),'n_lo':len(lo),
   'random_mean':float(rand.mean()),'random_q95':float(np.quantile(rand,.95)),
   'chance':1/len(multi),'n_plants':int(k.sum()),'n_countries':len(multi)}
np.save('rand_acc.npy',rand); json.dump(D,open('decode.json','w'),indent=1); print(json.dumps(D,indent=1))
